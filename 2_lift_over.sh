#!/usr/bin/env bash

#########################
## eQTLs
#########################
# get rsIDs & pull hg38 coordinates
mkdir data/intermediary
zcat data/eqtls/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz | awk -F"\t" 'NR > 1 {print $2}' | sort -u > data/intermediary/eqtl_rsids.txt 
zgrep -wFf data/intermediary/eqtl_rsids.txt data/dbSNP/GCF_000001405.40.gz > data/intermediary/eQTLGen_matched_snps_hg38.vcf
awk 'BEGIN {OFS="\t"; print "CHR", "POS", "ID"} {print $1, $2, $3}' data/intermediary/eQTLGen_matched_snps_hg38.vcf > data/intermediary/eQTLs_hg38.txt

# Lift over original file
python add_coords.py data/intermediary/eQTLs_hg38.txt \
                     data/dbSNP/GCF_000001405.40_GRCh38.p14_assembly_report.txt \
                     data/eqtls/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt \
                     eQTLGen

#########################
## bQTLs
#########################
# get rsIDs & pull hg38 coordinates
python get_baalnf_snps.py
zgrep -wFf data/intermediary/baal-nf-high-quality-snps.txt data/dbSNP/GCF_000001405.40.gz > data/intermediary/baalnf_matched_snps_hg38.vcf
awk 'BEGIN {OFS="\t"; print "CHR", "POS", "ID"} {print $1, $2, $3}' data/intermediary//baalnf_matched_snps_hg38.vcf > data/intermediary/baalnf_snps_hg38.txt

# Lift over original file
python add_coords.py data/intermediary/baalnf_snps_hg38.txt \
                     data/dbSNP/GCF_000001405.40_GRCh38.p14_assembly_report.txt \
                     /exports/igmm/eddie/ponting-lab/breeshey/data/baal-nf-all-tfs/compiled \
                     baal-nf

#########################
## transactors
#########################
# get rsIDs & pull hg19 coordinates
awk -F"\t" '{print $1}' data/additional_transactors/*/*_transactors_filtered_hg38.txt | sort -u > data/intermediary/transactor_rsids.txt
zgrep -v "^#" data/dbSNP/GCF_000001405.25.gz | awk '
    BEGIN {
        while ((getline < "data/intermediary/transactor_rsids.txt") > 0) 
            rsids[$1] = 1
        close("data/intermediary/transactor_rsids.txt")
    }
    $3 in rsids { print }
' > data/intermediary/transactor_matched_snps_hg19.vcf
awk '
BEGIN {
    OFS = "\t"
    print "ID", "CHR", "POS", "MAJOR", "MINOR"
}
{
    # Skip header lines
    if ($0 ~ /^#/) next

    # Extract fields from the VCF
    chr = $1
    pos = $2
    id = $3
    ref = $4
    alts = $5
    info = $8

    # Combine REF and ALT alleles into one array
    split(alts, alt_array, ",")
    alleles[1] = ref
    for (i = 1; i <= length(alt_array); i++) {
        alleles[i + 1] = alt_array[i]
    }

    # Extract the FREQ=1000Genomes data from the INFO column
    freq_data = ""
    if (match(info, /FREQ=1000Genomes:([^|;]+)/, match_array)) {
        freq_data = match_array[1] # Extract the frequencies for 1000Genomes
        split(freq_data, freq_array, ",") # Split into an array
    } else {
        next # Skip this record if no 1000Genomes data is found
    }

    # Map allele frequencies
    for (i = 1; i <= length(alleles); i++) {
        allele_freqs[i] = freq_array[i] + 0 # Convert string to number
    }

    # Find the two alleles with the highest frequencies
    max_freq1 = -1
    max_freq2 = -1
    major = ""
    minor = ""

    for (i = 1; i <= length(alleles); i++) {
        freq = allele_freqs[i]

        if (freq > max_freq1) {
            max_freq2 = max_freq1
            minor = major

            max_freq1 = freq
            major = alleles[i]
        } else if (freq > max_freq2) {
            max_freq2 = freq
            minor = alleles[i]
        }
    }

    # Print the SNP data with major and minor alleles
    if (major != "" && minor != "") {
        print id, chr, pos, major, minor
    }
}' data/intermediary/transactor_matched_snps_hg19.vcf > data/intermediary/transactors_hg19_with_major_minor.txt

# Lift over original files & add allele information
python add_coords.py data/intermediary/transactors_hg19_with_major_minor.txt \
                     data/dbSNP/GCF_000001405.25_GRCh37.p13_assembly_report.txt \
                     data/additional_transactors \
                     transactors
