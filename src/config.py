# src/config.py
import os


# ==========================================
# 路徑
# ==========================================

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

# FASTA檔路徑
REF_DIR = os.path.join(BASE_DIR, "data", "FASTA_data_from_HW1")

# SNV前後128 bases的檔案
SEQUENCES_FILE = os.path.join(OUTPUT_DIR, "all_sequences.tsv")

# 存CD-Hit暫存檔
INTERIM_DIR = os.path.join(BASE_DIR, "output", "interim")

# CD-HIT 相關檔案
CDHIT_FASTA_FILE = os.path.join(INTERIM_DIR, "for_cdhit.fasta")     # 轉檔後的 FASTA
CDHIT_OUTPUT_PREFIX = os.path.join(INTERIM_DIR, "cdhit_out")        # CD-HIT 輸出前綴

PROCESSED_DIR = os.path.join(BASE_DIR, "output", "processed")

# CD-HIT 產生的 Cluster 檔案 (.clstr)
# 注意: CD-HIT 會自動在 output prefix 後面加上 .clstr
CDHIT_CLSTR_FILE = CDHIT_OUTPUT_PREFIX + ".clstr"

# 切分好的資料集路徑
TRAIN_FILE = os.path.join(PROCESSED_DIR, "train.csv")
VAL_FILE   = os.path.join(PROCESSED_DIR, "val.csv")
TEST_FILE  = os.path.join(PROCESSED_DIR, "test.csv")

# ==========================================
# Header
# ==========================================

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

# ==========================================
# RDDL
# ==========================================

TASK_NAME = "ClinVar"          # RDDL 任務名稱
RDDL_DIR_NAME = "RDDL"         # RDDL 套件所在的資料夾名稱
DATA_SOURCE_DIR = "data"       # 原始 CSV 存放位置

# 設定 RDDL 任務目標路徑: ./RDDL/ClinVar/
RDDL_TASK_DIR = os.path.join(BASE_DIR, RDDL_DIR_NAME, TASK_NAME)
POS_DIR = os.path.join(RDDL_TASK_DIR, 'USER_pos')
NEG_DIR = os.path.join(RDDL_TASK_DIR, 'USER_neg')
SPLIT_INFO_DIR = os.path.join(RDDL_TASK_DIR, 'RDDL_splitting_info')

# 原始資料來源路徑
CSV_SOURCE_PATH = PROCESSED_DIR

# DNA 轉換字典
DNA_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3}