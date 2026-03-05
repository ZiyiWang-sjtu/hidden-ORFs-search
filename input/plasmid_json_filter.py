import json
import argparse
from Bio import SeqIO
import io

def get_rev_comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq.upper()))

def load_fasta_robust(fasta_path):
    target_genes = []
    try:
        with open(fasta_path, 'rb') as f:
            content = f.read().decode('utf-8', errors='ignore')
        
        fasta_io = io.StringIO(content)
        for record in SeqIO.parse(fasta_io, "fasta"):
            clean_seq = "".join(str(record.seq).split()).upper()
            if clean_seq:
                target_genes.append({
                    "id": record.id,
                    "seq": clean_seq
                })
        return target_genes
    except Exception as e:
        print(f"Fail to read the FASTA file: {e}")
        return []

def main():
    parser = argparse.ArgumentParser(description="Filter plasmids containing target genes")
    parser.add_argument("-i", "--input", required=True, help="Input JSON file")
    parser.add_argument("-f", "--fasta", required=True, help="Target gene FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output JSON file")
    args = parser.parse_args()

    targets = load_fasta_robust(args.fasta)
    if not targets:
        print("No valid gene sequences found. Please check the FASTA file.")
        return
    print(f"Loaded {len(targets)} target sequences.")

    with open(args.input, 'r', encoding='utf-8') as f:
        data = json.load(f)
        plasmids = data.get("plasmids", [])

    results = []
    for p in plasmids:
        found = False
        full_seq_list = p.get("sequences", {}).get("public_addgene_full_sequences", [])
        
        for s_obj in full_seq_list:
            original_seq = s_obj.get("sequence", "").upper()
            circular_seq = original_seq + original_seq
            
            for t in targets:
                fwd = t["seq"]
                rev = get_rev_comp(fwd)
                
                if fwd in circular_seq or rev in circular_seq:
                    p["matched_by"] = t["id"]  
                    found = True
                    break
            if found: break
        
        if found:
            results.append(p)

    with open(args.output, 'w', encoding='utf-8') as f:
        json.dump({"count": len(results), "plasmids": results}, f, indent=2, ensure_ascii=False)

    print(f"Processing complete! Matched {len(results)} plasmids.")

if __name__ == "__main__":
    main()
