import os
import sys
import csv
from collections import defaultdict


# Add src to path to import config if needed
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import config

def parse_gtf_attributes(attribute_string):
    """Parse GTF attribute string into a dictionary without regex."""
    attributes = {}
    # Split by semicolon to get key-value pairs
    # Example: transcript_id "NM_001"; gene_id "GeneA";
    parts = attribute_string.strip().split(';')
    for part in parts:
        part = part.strip()
        if not part:
            continue
        # Split by first space
        if ' ' in part:
            key, value = part.split(' ', 1)
            # Remove quotes from value if present
            value = value.strip('"')
            attributes[key] = value
    return attributes

def load_gtf(gtf_path):
    """
    Load GTF file and organize by transcript_id -> feature_type -> list of intervals.
    Returns: dict[transcript_id][feature_type] = [(start, end), ...]
    """
    print(f"Loading GTF file: {gtf_path} ...")
    transcript_data = defaultdict(lambda: defaultdict(list))
    
    try:
        with open(gtf_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                
                feature_type = parts[2]
                start = int(parts[3])
                end = int(parts[4])
                attributes_str = parts[8]
                
                attributes = parse_gtf_attributes(attributes_str)
                transcript_id = attributes.get('transcript_id')
                
                if transcript_id:
                    transcript_data[transcript_id][feature_type].append((start, end))
                    
    except FileNotFoundError:
        print(f"Error: GTF file not found at {gtf_path}")
        sys.exit(1)
        
    print(f"Loaded {len(transcript_data)} transcripts from GTF.")
    return transcript_data

def validate_locations(tsv_path, gtf_data):
    """
    Validate that SNVs in the TSV fall within the specified transcript feature regions.
    """
    print(f"Validating variants in: {tsv_path} ...")
    
    if not os.path.exists(tsv_path):
        print(f"Error: TSV file not found at {tsv_path}")
        sys.exit(1)
        
    mismatches = []
    missing_transcripts = set()
    total_checked = 0
    valid_count = 0
    
    with open(tsv_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        
        # Check required columns
        required_cols = ['chrom', 'pos', 'transcript_id', 'feature_type']
        if not reader.fieldnames:
             print("Error: TSV file is empty or has no header.")
             sys.exit(1)
             
        for col in required_cols:
            if col not in reader.fieldnames:
                print(f"Error: Missing required column '{col}' in TSV.")
                sys.exit(1)
        
        for idx, row in enumerate(reader):
            total_checked += 1
            # row index in file is idx + 2 (header is 1)
            line_num = idx + 2
            
            try:
                chrom = row['chrom']
                pos = int(row['pos'])
                transcript_id = row['transcript_id']
                feature_type = row['feature_type']
            except ValueError:
                print(f"Warning: Invalid data at line {line_num}: {row}")
                continue
            
            # Skip if transcript_id is missing or dot
            if not transcript_id or transcript_id == '.':
                continue
                
            if transcript_id not in gtf_data:
                missing_transcripts.add(transcript_id)
                continue
                
            # Get intervals for this feature type
            # Normalize feature type if needed (e.g., 5'UTR -> 5UTR)
            intervals = gtf_data[transcript_id].get(feature_type, [])
            
            # If no exact match, try mapping
            if not intervals:
                if feature_type == "5'UTR":
                    intervals = gtf_data[transcript_id].get('5UTR', [])
                elif feature_type == "3'UTR":
                    intervals = gtf_data[transcript_id].get('3UTR', [])
            
            if not intervals:
                mismatches.append({
                    'row': line_num,
                    'transcript_id': transcript_id,
                    'pos': pos,
                    'feature_type': feature_type,
                    'reason': f"Feature '{feature_type}' not found for transcript in GTF"
                })
                continue
                
            # Check if pos is in any interval
            is_valid = False
            for start, end in intervals:
                if start <= pos <= end:
                    is_valid = True
                    break
            
            if is_valid:
                valid_count += 1
            else:
                mismatches.append({
                    'row': line_num,
                    'transcript_id': transcript_id,
                    'pos': pos,
                    'feature_type': feature_type,
                    'reason': f"Position {pos} not in {feature_type} intervals: {intervals}"
                })

    # Report results
    print("\n" + "="*40)
    print("Validation Results")
    print("="*40)
    print(f"Total variants checked: {total_checked}")
    print(f"Valid variants: {valid_count}")
    print(f"Mismatches found: {len(mismatches)}")
    print(f"Transcripts missing in GTF: {len(missing_transcripts)}")
    
    if missing_transcripts:
        print(f"(First 10 missing transcripts: {list(missing_transcripts)[:10]})")
        
    if mismatches:
        print("\nMismatch Details (First 20):")
        for m in mismatches[:20]:
            print(f"Line {m['row']}: {m['transcript_id']} ({m['feature_type']}) at {m['pos']} - {m['reason']}")
            
    if len(mismatches) == 0 and len(missing_transcripts) == 0:
        print("\nSUCCESS: All checked variants match their transcript features!")
    else:
        print("\nWARNING: Some discrepancies were found.")

if __name__ == "__main__":
    # Default paths based on project structure
    project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    default_tsv = os.path.join(project_root, "output", "matched_variants.tsv")
    default_gtf = os.path.join(project_root, "data", "hg38.ncbiRefSeq.gtf")
    
    tsv_file = sys.argv[1] if len(sys.argv) > 1 else default_tsv
    gtf_file = sys.argv[2] if len(sys.argv) > 2 else default_gtf
    
    print(f"Target TSV: {tsv_file}")
    print(f"Target GTF: {gtf_file}")
    
    if not os.path.exists(gtf_file):
        print(f"GTF file not found: {gtf_file}")
        # We continue to let load_gtf fail gracefully or user might have provided args
    
    gtf_data = load_gtf(gtf_file)
    validate_locations(tsv_file, gtf_data)
