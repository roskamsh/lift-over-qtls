import pandas as pd
import os
import sys

def add_genome_coordinates(snp_positions_file, reference_file, original_file, datatype):
    snp_positions=pd.read_csv(snp_positions_file, delim_whitespace=True)
    reference=pd.read_csv(reference_file, comment="#", sep="\t", header=None)
    reference.columns=["Name","Role","Molecule","Molecule_type","GenBank_Accn","Relationship",
                    "RefSeq_Accn","Assembly-unit","Sequence-length","UCSC-name"]
    reference = reference[reference["Role"]=="assembled-molecule"]
    merged = snp_positions.merge(reference, left_on = "CHR", right_on = "RefSeq_Accn")

    ## Now join back with original file
    if (datatype=="eQTLGen"):
        hg38_coords = merged[["ID","Molecule","POS"]].copy()
        hg38_coords.columns = ["ID","CHR_hg38","POS_hg38"]
        outname=os.path.basename(original_file).replace(".txt","") + '_hg38.txt'
        original=pd.read_csv(original_file,sep="\t")
        original = original.drop(columns=["SNPChr","SNPPos"])
        hg38_coords.columns = ["SNP","SNPChr","SNPPos"]
        original_annot = original.merge(hg38_coords, on="SNP")
        os.makedirs(os.path.join("data","lifted","eqtls"), exist_ok=True)
        original_annot.to_csv(os.path.join("data","lifted","eqtls",outname), index = False, sep="\t")
    elif (datatype=="baal-nf"):
        hg38_coords = merged[["ID","Molecule","POS"]].copy()
        hg38_coords.columns = ["ID","CHR_hg38","POS_hg38"]
        directory=original_file 
        for filename in os.listdir(directory):
            if filename.endswith(".csv"):  # Ensure it's a CSV file
                print("reading file", filename)
                filepath = os.path.join(directory, filename)
                outname=filename.replace(".csv","") + ".highQuality.hg38.csv"
                df = pd.read_csv(filepath)
                
                if "ASB_quality" in df.columns and "ID" in df.columns:
                    df["ASB_quality"] = df["ASB_quality"].astype(str)
                    filtered_df = df[df["ASB_quality"] == "High"]
                    if not filtered_df.empty:
                        filtered_df = filtered_df.drop(columns=["CHROM","POS"])
                        hg38_coords.columns = ["ID","CHROM","POS"]
                        updated = filtered_df.merge(hg38_coords, on = "ID")
                        updated["CHROM"] = 'chr' + updated['CHROM'].astype(str) # to exactly match baal-nf output
                        os.makedirs(os.path.join("data","lifted","baal-nf"), exist_ok=True)
                        updated.to_csv(os.path.join("data","lifted","baal-nf",outname), index = False)
                else:
                    print(f"Skipping {filename}: Missing required columns")
    elif (datatype=="transactors"):
        hg19_coords = merged[["ID","Molecule","POS","MAJOR","MINOR"]].copy()
        hg19_coords.columns = ["ID","CHR_hg19","POS_hg19","MAJOR","MINOR"]
        directory=original_file
        for tf in os.listdir(directory):
            file=f"{directory}/{tf}/{tf}_transactors_filtered_hg38.txt"
            print("reading file", file)
            df = pd.read_csv(file, sep = "\t")
            updated = df.merge(hg19_coords, left_on = "SNPS", right_on = "ID")
            # Format AssessedAllele and OtherAllele information
            updated["AssessedAllele"] = updated["RISK.ALLELE"]
            # Add the OtherAllele column, which is always equal to the other allele
            updated["OtherAllele"] = updated.apply(
                lambda row: row["MAJOR"] if row["AssessedAllele"] == row["MINOR"] else row["MINOR"],
                axis=1
            )
            updated["CHR_POS"] = updated["CHR_POS"].astype(int)
            updated["Pvalue"] = 0 # Set Pvalue column to zero for prioritization in ld-clump pipeline
            updated["GeneSymbol"] = tf
            # write hg19 table with alleles, formatted for ld-clump pipeline
            hg19_out = updated[['ID', 'CHR_hg19','POS_hg19', 'MAJOR', 'MINOR', 'AssessedAllele', 'OtherAllele','GeneSymbol',
                                'Pvalue','STRONGEST.SNP.RISK.ALLELE', 'RISK.ALLELE.FREQUENCY', 'DISEASE.TRAIT', 'MAPPED_TRAIT',
                                'MAPPED_TRAIT_URI', 'P.VALUE', 'STUDY', 'STUDY.ACCESSION', 'INITIAL.SAMPLE.SIZE', 
                                'DATE.ADDED.TO.CATALOG', 'OR.or.BETA', 'PVALUE_MLOG', 'RISK.ALLELE', 'SAMPLE.SIZE']].copy()
            hg19_out.columns = ['SNPChr','SNPChr','SNPPos','MajorAllele','MinorAllele','AssessedAllele', 'OtherAllele','GeneSymbol',
                                'Pvalue','STRONGEST.SNP.RISK.ALLELE', 'RISK.ALLELE.FREQUENCY', 'DISEASE.TRAIT', 'MAPPED_TRAIT',
                                'MAPPED_TRAIT_URI', 'P.VALUE', 'STUDY', 'STUDY.ACCESSION', 'INITIAL.SAMPLE.SIZE', 
                                'DATE.ADDED.TO.CATALOG', 'OR.or.BETA', 'PVALUE_MLOG', 'RISK.ALLELE', 'SAMPLE.SIZE']
            outname_hg19=os.path.basename(file).replace("_hg38.txt","_hg19_with_alleles.txt")
            os.makedirs(os.path.join("data","lifted","transactors",tf), exist_ok=True)
            hg19_out.to_csv(os.path.join("data","lifted","transactors",tf,outname_hg19), index=False,sep="\t")
            
            # write hg38 table with alleles, formatted for ld-clump pipeline
            hg38_out = updated[['ID', 'CHR_ID','CHR_POS', 'MAJOR', 'MINOR', 'AssessedAllele', 'OtherAllele','GeneSymbol',
                                'Pvalue','STRONGEST.SNP.RISK.ALLELE', 'RISK.ALLELE.FREQUENCY', 'DISEASE.TRAIT', 'MAPPED_TRAIT',
                                'MAPPED_TRAIT_URI', 'P.VALUE', 'STUDY', 'STUDY.ACCESSION', 'INITIAL.SAMPLE.SIZE', 
                                'DATE.ADDED.TO.CATALOG', 'OR.or.BETA', 'PVALUE_MLOG', 'RISK.ALLELE', 'SAMPLE.SIZE']].copy()
            hg38_out.columns = ['SNPChr','SNPChr','SNPPos','MajorAllele','MinorAllele','AssessedAllele', 'OtherAllele','GeneSymbol',
                                'Pvalue','STRONGEST.SNP.RISK.ALLELE', 'RISK.ALLELE.FREQUENCY', 'DISEASE.TRAIT', 'MAPPED_TRAIT',
                                'MAPPED_TRAIT_URI', 'P.VALUE', 'STUDY', 'STUDY.ACCESSION', 'INITIAL.SAMPLE.SIZE', 
                                'DATE.ADDED.TO.CATALOG', 'OR.or.BETA', 'PVALUE_MLOG', 'RISK.ALLELE', 'SAMPLE.SIZE'] 
            outname_hg38=os.path.basename(file).replace("_hg38.txt","_hg38_with_alleles.txt")
            os.makedirs(os.path.join("data","lifted","transactors",tf), exist_ok=True)
            hg38_out.to_csv(os.path.join("data","lifted","transactors",tf,outname_hg38), index=False, sep="\t")

def main():
    # Check the number of arguments
    if len(sys.argv) != 5:
        print("Usage: python add_coords.py <snp-positions-file> <reference-file> <original-file> <data-type>")
        sys.exit(1)

    # Parse command-line arguments
    snp_positions_file_arg = sys.argv[1]
    reference_file_arg = sys.argv[2]
    original_file_arg = sys.argv[3]
    datatype_arg = sys.argv[4]

    # Call the function to process alleles
    add_genome_coordinates(snp_positions_file_arg, reference_file_arg, original_file_arg, datatype_arg)

if __name__ == "__main__":
    main()