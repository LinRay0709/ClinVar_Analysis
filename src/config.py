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