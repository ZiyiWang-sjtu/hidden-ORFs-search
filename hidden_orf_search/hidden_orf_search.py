import csv
import json
import re
import argparse
import sys
import io
import os
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio.Seq import Seq
from Bio import SeqIO
from typing import List, Tuple, Dict

# Clean DNA sequence (remove non-ATCGN characters and convert to uppercase)
def clean_sequence(seq_str: str) -> str:
    if not seq_str:
        return ""
    return re.sub(r'[^ATCGN]', '', str(seq_str).upper())

# Run BLASTn and return all hits above identity threshold
def find_all_blast(query_seq: str, subject_seq: str, identity_threshold: float = 96.0) -> List[Dict]:
    if not query_seq or not subject_seq:
        return []

    hits = []
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as q_file, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as s_file:
        
        q_file.write(f">query\n{query_seq}\n")
        s_file.write(f">subject\n{subject_seq}\n")
        q_path, s_path = q_file.name, s_file.name

    try:
        cmd = [
            "blastn",
            "-query", q_path,
            "-subject", s_path,
            "-perc_identity", str(identity_threshold),
            "-outfmt", "6 sstart send pident length", 
            "-task", "blastn-short" if len(query_seq) < 30 else "blastn"
        ]
        
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = process.communicate()

        for line in stdout.strip().split('\n'):
            if not line: continue
            parts = line.split('\t')
            s_start = int(parts[0])
            s_end = int(parts[1])
            pident = float(parts[2])
            length = int(parts[3])
            
            strand = "+" if s_start < s_end else "-"
            start_idx = min(s_start, s_end) - 1
            end_idx = max(s_start, s_end)
            
            hits.append({
                "start": start_idx, 
                "end": end_idx, 
                "strand": strand,
                "pident": pident, 
                "length": length
            })
    finally:
        if os.path.exists(q_path): os.remove(q_path)
        if os.path.exists(s_path): os.remove(s_path)
    return hits

# Safely read a file using multiple encoding attempts
def read_file_safe(file_path: str) -> str:
    encodings = ['utf-8-sig', 'gbk', 'utf-8', 'latin1']
    for enc in encodings:
        try:
            with open(file_path, 'r', encoding=enc) as f:
                return f.read()
        except: continue
    sys.exit(1)

# Load promoter and polyA feature sequences from CSV
def load_feature_library(csv_path: str) -> Dict[str, List[Dict]]:
    feature_db = {"promoter": [], "polyA": []}
    content = read_file_safe(csv_path)
    lines = content.splitlines()
    reader = csv.DictReader(lines)
    raw_fields = [f.strip() for f in reader.fieldnames]
    cols_map = {f.lower(): f for f in raw_fields}
    
    name_col = next((cols_map[v] for v in ['name', 'label', 'feature_name'] if v in cols_map), None)
    type_col = next((cols_map[v] for v in ['feature_type', 'type'] if v in cols_map), None)
    seq_col = next((cols_map[v] for v in ['sequence', 'seq'] if v in cols_map), None)

    for row in reader:
        f_type = (row.get(type_col) or "").strip().lower()
        f_name = (row.get(name_col) or "Unnamed").strip()
        seq = clean_sequence((row.get(seq_col) or "").strip())
        if not seq: continue
        
        entry = {"name": f_name, "seq": seq}
        if 'promoter' in f_type: 
            feature_db["promoter"].append(entry)
        elif any(kw in f_type for kw in ['polya', 'term', 'stop', 'ltr']): 
            feature_db["polyA"].append(entry)
    return feature_db

# Map feature sequences onto circular plasmid
def map_features_to_plasmid(plasmid_seq: str, feature_entries: List[Dict], threshold: float) -> List[Dict]:
    L = len(plasmid_seq)
    dseq = plasmid_seq + plasmid_seq 
    mapped_results = []
    
    for entry in feature_entries:
        f_seq = entry["seq"]
        q_len = len(f_seq)
        hits = find_all_blast(f_seq, dseq, threshold)
        for h in hits:
            if h["start"] < L and h["length"] >= (q_len * 0.9):
                mapped_results.append({
                    "name": entry["name"],
                    "start": h["start"],
                    "end": h["end"]
                })
    return mapped_results

# Extract full plasmid sequence from JSON
def extract_plasmid_seq(plasmid: Dict) -> str:
    seqs_dict = plasmid.get("sequences", {})
    for key in ["public_addgene_full_sequences", "public_user_full_sequences", "full_sequences"]:
        seq_list = seqs_dict.get(key)
        if seq_list and isinstance(seq_list, list) and len(seq_list) > 0:
            item = seq_list[0]
            return item.get("sequence", "") if isinstance(item, dict) else str(item)
    return ""

# Scan ORFs in both strands of circular sequence
def scan_orfs(dseq: str, L: int, min_aa: int) -> List[Dict]:
    start_codons, stop_codons = {"ATG"}, {"TAA", "TAG", "TGA"}
    orfs, total_len = [], len(dseq)
    for frame in (0, 1, 2):
        i = frame
        while i < total_len - 2:
            if dseq[i:i+3] in start_codons:
                j = i + 3
                while j < total_len - 2:
                    if dseq[j:j+3] in stop_codons:
                        length_aa = (j + 3 - i) // 3
                        if length_aa >= min_aa:
                            orfs.append({"strand": "+", "frame": frame + 1, "start": i, "end": j + 3, "length_aa": length_aa})
                        i = j + 3
                        break
                    j += 3
                else: i += 3
            else: i += 3
    s_rev = str(Seq(dseq).reverse_complement())
    for frame in (0, 1, 2):
        i = frame
        while i < total_len - 2:
            if s_rev[i:i+3] in start_codons:
                j = i + 3
                while j < total_len - 2:
                    if s_rev[j:j+3] in stop_codons:
                        actual_end, actual_start = total_len - i, total_len - (j + 3)
                        length_aa = (j + 3 - i) // 3
                        if length_aa >= min_aa:
                            orfs.append({"strand": "-", "frame": frame + 1, "start": actual_start, "end": actual_end, "length_aa": length_aa})
                        i = j + 3
                        break
                    j += 3
                else: i += 3
            else: i += 3
    return [o for o in orfs if o["start"] < L]

# Classify ORF overlap relative to target gene
def classify_overlap(orf: Tuple[int, int], target: Tuple[int, int]) -> str:
    a, b, c, d = orf[0], orf[1], target[0], target[1]
    if b <= c or a >= d: return "outside_target_gene"
    if a >= c and b <= d: return "inside_target_gene"
    if a <= c and b >= d: return "spanning_target_gene"
    return "partial_target_gene"

# Calculate upstream circular distance
def get_circular_upstream_dist(feature_end, orf_start, L):
    if orf_start >= feature_end:
        return orf_start - feature_end
    else:
        return orf_start + (L - feature_end)

# Calculate downstream circular distance
def get_circular_downstream_dist(orf_start, feature_start, L):
    if feature_start >= orf_start:
        return feature_start - orf_start
    else:
        return (L - orf_start) + feature_start

# Find closest circular feature (upstream or downstream)
def find_closest_circular_feature(feature_list, pos_mod, L, mode='upstream'):
    if not feature_list: return None, None
    best_feat, min_dist = None, float('inf')
    for f in feature_list:
        if mode == 'upstream':
            d = get_circular_upstream_dist(f["end"], pos_mod, L)
        else:
            d = get_circular_downstream_dist(pos_mod, f["start"], L)
        if d < min_dist:
            min_dist = d
            best_feat = f
    return best_feat, min_dist

# Analyze ORFs relative to target gene and regulatory elements
def plasmid_orf_analyzer(p_seq, t_seq, proms_mapped, polys_mapped, min_aa, max_dist_prom, max_dist_poly, threshold):
    L, dseq = len(p_seq), p_seq + p_seq 
    all_hits = find_all_blast(t_seq, dseq, threshold)
    target_hits = []
    
    for h in all_hits:
        if h["start"] < L and h["length"] >= (len(t_seq) * 0.9):
            target_hits.append({
                "range": (h["start"], h["end"]),
                "start_abs": h["start"],
                "end_abs": h["end"],
                "strand": h["strand"]
            })

    orfs, results = scan_orfs(dseq, L, min_aa), []
    for orf in orfs:
        start_dseq, end_dseq = orf["start"], orf["end"]
        overlap_status = "outside_target_gene"
        tg_frame = "-"
        
        for g in target_hits:
            status = classify_overlap((start_dseq, end_dseq), g["range"])
            
            if status != "outside_target_gene":
                overlap_status = status
                
                if orf["strand"] == "+": 
                    base_phase = (start_dseq - g["start_abs"]) % 3
                else:
                    base_phase = (g["end_abs"] - end_dseq) % 3

                display_phase = base_phase + 1
                prefix = "+" if orf["strand"] == g["strand"] else "-"
                tg_frame = f"{prefix}{display_phase}"
                break

        pos_mod = start_dseq % L
        closest_p, dist_p = find_closest_circular_feature(proms_mapped, pos_mod, L, 'upstream')
        closest_a, dist_a = find_closest_circular_feature(polys_mapped, end_dseq % L, L, 'downstream')
        
        has_promoter_nearby = (dist_p is not None and dist_p <= max_dist_prom)
        has_polyA_nearby = (dist_a is not None and dist_a <= max_dist_poly)
        
        expressible = (overlap_status != "outside_target_gene") and has_promoter_nearby and has_polyA_nearby
        
        is_risk = False
        if expressible:
            if tg_frame != "+1":
                is_risk = True
        
        results.append({
            "TGframe": tg_frame,
            "orf_frame": orf["frame"],   
            "strand": orf["strand"], 
            "start": pos_mod + 1, 
            "stop": (end_dseq - 1) % L + 1, 
            "length_aa": orf["length_aa"],
            "overlap": overlap_status,
            "has_promoter": has_promoter_nearby,
            "has_polyA": has_polyA_nearby,
            "dist_promoter": dist_p if dist_p is not None else "-",
            "dist_polyA": dist_a if dist_a is not None else "-",
            "expressible": expressible,
            "is_risk": is_risk
        })
    return results

# Export results to CSV
def export_to_csv(data: List[Dict], filename: str):
    if not data: return
    cols = [
        "plasmid_id", "plasmid_name", "TGframe", "orf_frame", "strand", "start", "stop", "length_aa", 
        "overlap", "has_promoter", "has_polyA", "dist_promoter", "dist_polyA", "expressible", "is_risk"
    ]
    with open(filename, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=cols)
        writer.writeheader()
        writer.writerows(data)

# Process a single plasmid and return ORF annotation
def process_single_plasmid(idx: int, plasmid: Dict, target_seq: str, feature_db: Dict, 
                           min_aa: int, max_prom: int, max_poly: int, identity: float) -> List[Dict]:
    p_seq_raw = extract_plasmid_seq(plasmid)
    if not p_seq_raw: return []
    p_seq = clean_sequence(p_seq_raw)
    
    proms_mapped = map_features_to_plasmid(p_seq, feature_db["promoter"], identity)
    polys_mapped = map_features_to_plasmid(p_seq, feature_db["polyA"], identity)
    
    results = plasmid_orf_analyzer(p_seq, target_seq, proms_mapped, polys_mapped, 
                                    min_aa, max_prom, max_poly, identity)
    for r in results:
        r["plasmid_id"] = plasmid.get("id", f"{idx}")
        r["plasmid_name"] = plasmid.get("name", "Unknown")
    return results

# Main execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", type=str, required=True)
    parser.add_argument("--fa", type=str, required=True)
    parser.add_argument("--features", type=str, required=True)
    parser.add_argument("--identity", type=float, default=96.0)
    parser.add_argument("--min_aa", type=int, default=100)
    parser.add_argument("--max_prom", type=int, default=10000, help="Maximum distance from the promoter")
    parser.add_argument("--max_poly", type=int, default=10000, help="Maximum distance from the polyA")
    parser.add_argument("--out", type=str, default="results.csv")
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    args = parser.parse_args()

    target_seq = clean_sequence(str(next(SeqIO.parse(io.StringIO(read_file_safe(args.fa)), "fasta")).seq))
    feature_db = load_feature_library(args.features)
    
    with open(args.json, 'r') as f:
        plasmids = json.load(f).get("plasmids", [])
    
    print(f"Starting processing of {len(plasmids)} plasmids")

    all_results = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(process_single_plasmid, i, p, target_seq, feature_db, 
                            args.min_aa, args.max_prom, args.max_poly, args.identity): p 
            for i, p in enumerate(plasmids, 1)
        }
        for i, future in enumerate(as_completed(futures), 1):
            try:
                res = future.result()
                if res: all_results.extend(res)
            except Exception as e: print(f"\n[Error] {e}")
            sys.stdout.write(f"\rProcessing progress: [{i}/{len(plasmids)}]"); sys.stdout.flush()

    print("\nSave data...")
    export_to_csv(all_results, args.out)
    print(f"Analysis complete! File saved to:: {args.out}")
