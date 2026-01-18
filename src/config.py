# src/config.py
import os

# 專案根目錄 (假設 src 在 project_root/src)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# 輸入檔案路徑 (請根據實際情況修改)
BED_FILE = os.path.join(BASE_DIR, "data_0107", "MANE_transcript.bed")
VCF_FILE = os.path.join(BASE_DIR, "data_new", "clinvar_snv_new.vcf.gz")

# 輸出檔案路徑
OUTPUT_DIR = os.path.join(BASE_DIR, "output")
# 確保輸出目錄存在
os.makedirs(OUTPUT_DIR, exist_ok=True)

OUTPUT_FILE = os.path.join(OUTPUT_DIR, "matched_variants.tsv")

# 定義輸出的欄位順序 (Header)
CSV_HEADER = [
    "chrom", "pos", "ref", "alt",
    "transcript_id", "feature_type", "gene_name",
    "clnsig_category", "mc_category",
    "clnsig_raw", "mc_raw"
]

# all_sequence.tsv的Header
SEQ_CSV_HEADER = [
    "snv_key", 
    "chrom", "pos", "ref", "alt", 
    "feature_type", "mc_category", "clnsig_category", 
    "ref_seq", "alt_seq"
]

# prepare_cdhit.py 需要讀取的欄位
CDHIT_PREPARE_READ_COLS = [
    "snv_key", 
    "feature_type", "mc_category", "clnsig_category", 
    "alt_seq"
]

# FASTA檔路徑
REF_DIR = os.path.join(BASE_DIR, "data", "FASTA_data_from_HW1")

# SNV前後128 bases的檔案
SEQUENCES_FILE = os.path.join(OUTPUT_DIR, "all_sequences.tsv")

# 存CD-Hit暫存檔
INTERIM_DIR = os.path.join(BASE_DIR, "output", "interim")

# CD-HIT 相關檔案
CDHIT_FASTA_FILE = os.path.join(INTERIM_DIR, "for_cdhit.fasta")     # 轉檔後的 FASTA
CDHIT_OUTPUT_PREFIX = os.path.join(INTERIM_DIR, "cdhit_out")        # CD-HIT 輸出前綴