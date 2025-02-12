import os
import pandas as pd

def get_unique_ids_from_csv(directory):
    unique_ids = set()
    
    for filename in os.listdir(directory):
        if filename.endswith(".csv"):  # Ensure it's a CSV file
            print("reading file", filename)
            filepath = os.path.join(directory, filename)
            df = pd.read_csv(filepath)
            
            if "ASB_quality" in df.columns and "ID" in df.columns:
                filtered_ids = df[df["ASB_quality"] == "High"]["ID"].dropna().unique()
                unique_ids.update(filtered_ids)
            else:
                print(f"Skipping {filename}: Missing required columns")
    
    return list(unique_ids)

def main():
    # Example usage
    directory_path = "/exports/igmm/eddie/ponting-lab/breeshey/data/baal-nf-all-tfs/compiled/"
    unique_ids_list = get_unique_ids_from_csv(directory_path)

    out = pd.DataFrame({ 'ID' : unique_ids_list})
    os.makedirs(os.path.join("data","intermediary"), exist_ok = True)
    out.to_csv("data/intermediary/baal-nf-high-quality-snps.txt", header = None, index = False)

if __name__ == "__main__":
    main()