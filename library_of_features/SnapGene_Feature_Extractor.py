import os
import argparse
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

def extract_features_from_dna_files(input_folders, output_file="feature_library.csv"):
    feature_data = []

    ALLOWED_TYPES_MAP = {
        "features_promoter_from_snapgene": {"promoter", "5' LTR", "LTR", "5'LTR"},
        "features_polya_from_snapgene": {"polyA_signal", "3' LTR", "LTR", "3'LTR", "terminator", "polyA_site"},
        "features_terminator_from_snapgene": {"polyA_signal", "3' LTR", "LTR", "3'LTR", "terminator", "polyA_site"}
    }

    STANDARD_TYPE_NAME = {
        "features_promoter_from_snapgene": "promoter",
        "features_polya_from_snapgene": "polyA_signal",
        "features_terminator_from_snapgene": "polyA_signal"
    }

    for folder_path in input_folders:
        normalized_path = os.path.abspath(folder_path).rstrip(os.sep)
        if not os.path.exists(normalized_path):
            print(f"Warning: Folder {normalized_path} does not exist. Skipping...")
            continue

        folder_name = os.path.basename(normalized_path)
        
        if folder_name not in ALLOWED_TYPES_MAP:
            print(f"Skipping folder '{folder_name}' as it does not match defined feature folder names.")
            continue
            
        allowed_types = ALLOWED_TYPES_MAP[folder_name]
        target_type_label = STANDARD_TYPE_NAME[folder_name]

        dna_files = [f for f in os.listdir(normalized_path) if f.endswith('.dna')]
        
        if not dna_files:
            print(f"No .dna files found in {normalized_path}.")
            continue

        print(f"Processing {len(dna_files)} files in: {folder_name} -> Mapping to Feature_Type: {target_type_label}")

        for filename in tqdm(dna_files, desc=f"Scanning {folder_name}"):
            file_path = os.path.join(normalized_path, filename)
            
            try:
                with open(file_path, "rb") as handle:
                    record = SeqIO.read(handle, "snapgene")
                    
                    for feature in record.features:
                        raw_feat_type = feature.type
                        
                        if raw_feat_type not in allowed_types:
                            continue
            
                        feat_name = feature.qualifiers.get('label', 
                                    feature.qualifiers.get('note', 
                                    feature.qualifiers.get('gene', ['Unknown'])))[0]
                    
                        feat_seq = str(feature.extract(record.seq))
                        
                        feature_data.append({
                            "Source_File": filename,
                            "Feature_Name": feat_name,
                            "Feature_Type": target_type_label,
                            "Sequence_Length": len(feat_seq),
                            "Sequence": feat_seq
                        })
            except Exception as e:
                print(f"Error processing {filename}: {e}")

    df = pd.DataFrame(feature_data)
    
    if not df.empty:
        df.to_csv(output_file, index=False, encoding='utf-8')
        print(f"\nExtraction complete! Total features collected: {len(df)}")
        print(f"Library saved to: {os.path.abspath(output_file)}")
    else:
        print("\nNo matching feature data was extracted. Please ensure folder names match exactly.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract features from specific SnapGene folders with strict categorization.")
    parser.add_argument("--input", "-i", nargs="+", required=True, help="Path(s) to folders (e.g. features_promoter_from_snapgene).")
    parser.add_argument("--output", "-o", default="feature_library.csv", help="Output CSV name.")
    
    args = parser.parse_args()
    extract_features_from_dna_files(args.input, args.output)
