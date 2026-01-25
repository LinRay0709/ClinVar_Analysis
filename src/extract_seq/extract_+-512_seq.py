# src/extract_seq/extraction.py
import sys
import os
import pandas as pd
import time

# --- 路徑設定 ---
current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.dirname(current_dir)
project_root = os.path.dirname(src_dir)

if project_root not in sys.path:
    sys.path.append(project_root)

from src import config
# 引用隔壁鄰居的函式 (讀取 FASTA 的邏輯)
from src.extract_seq.ref_seq_consistency_check import load_fasta_sequence

# ==========================================
# 參數設定 (正式模式)
# ==========================================
IS_DEBUG = False            # 正式執行時改成 False
DEBUG_LIMIT = 10           # 只抓前 10 筆
CONTEXT_RADIUS = 512         # 前後各抓 5 bp (正式版請改 128)

# 輸出路徑邏輯
if IS_DEBUG:
    # Debug 模式下，避免覆蓋正式檔案，我們寫在另一個地方
    OUTPUT_PATH = os.path.join(config.OUTPUT_DIR, "debug_sequences.tsv")
else:
    # 正式模式下，讀取 config 中的設定
    OUTPUT_PATH = config.SEQUENCES_FILE
# ==========================================

def get_padded_sequence(seq_str, pos_0based, radius):
    """
    安全地抓取序列，遇到邊界自動補 'N'
    Args:
        seq_str: 整條染色體的序列字串
        pos_0based: 變異位點的索引 (0-based)
        radius: 前後要抓多長
    Returns:
        padded_seq: 補完 N 的序列
    """
    seq_len = len(seq_str)
    
    # 計算想要抓取的範圍
    start = pos_0based - radius
    end = pos_0based + radius + 1 -1 # Python 切片 end 是不包含的，所以要 +1，但抓共1024(512+1+511)
    
    # 初始化補白
    left_pad = ""
    right_pad = ""
    
    # 處理左邊界 (如果 start < 0)
    real_start = start
    if start < 0:
        left_pad = "@" * abs(start)
        real_start = 0
        
    # 處理右邊界 (如果 end > seq_len)
    real_end = end
    if end > seq_len:
        right_pad = "@" * (end - seq_len)
        real_end = seq_len
        
    # 抓取實際存在的序列
    actual_seq = seq_str[real_start:real_end]
    
    # 組合
    full_seq = left_pad + actual_seq + right_pad
    return full_seq

def main():
    print("="*60)
    print(f"序列提取程式 ({'DEBUG 模式' if IS_DEBUG else '正式模式'})")
    print(f"  - 抓取範圍: Ref 前後,總長度: {CONTEXT_RADIUS*2} bp")
    print(f"  - 預期總長: {CONTEXT_RADIUS * 2} bp")
    if IS_DEBUG:
        print(f"  - 限制筆數: 前 {DEBUG_LIMIT} 筆")
    print("="*60)

    # 1. 讀取配對資料
    if not os.path.exists(config.OUTPUT_FILE):
        print(f"錯誤: 找不到 {config.OUTPUT_FILE}")
        return

    df = pd.read_csv(config.OUTPUT_FILE, sep='\t', dtype={'chrom': str})
    
    # 2. 排序 (重要)
    def sort_key(x):
        try:
            return int(x)
        except:
            return 999
    df['_sort_key'] = df['chrom'].apply(sort_key)
    df = df.sort_values(by=['_sort_key', 'pos'])
    
    # --- DEBUG 模式的切片 ---
    if IS_DEBUG:
        print(f"[Debug] 僅保留前 {DEBUG_LIMIT} 筆資料進行測試...")
        df = df.head(DEBUG_LIMIT)
    # -----------------------

    # 準備寫入檔案
    # 使用 'w' 模式，並設定 header
    print(f"輸出檔案: {OUTPUT_PATH}")
    
    with open(OUTPUT_PATH, 'w') as f_out:
        # 寫入 Header
        header = config.SEQ_CSV_HEADER
        f_out.write("\t".join(header) + "\n")
        
        processed_count = 0
        
        # 3. 依染色體分組處理
        for chrom, group in df.groupby('chrom', sort=False):
            print(f"\n處理染色體: {chrom} (本批次共 {len(group)} 筆)")
            
            # 載入序列
            seq_str = load_fasta_sequence(chrom)
            if seq_str is None:
                continue
                
            for idx, row in group.iterrows():
                try:
                    pos_1based = int(row['pos'])
                    ref_base = row['ref']
                    alt_base = row['alt']
                    feat_type = str(row['feature_type'])
                    mc_cat = str(row['mc_category'])
                    clnsig_cat = str(row['clnsig_category'])
                    
                    # 建立唯一 ID
                    snv_key = f"{chrom}:{pos_1based}:{ref_base}:{alt_base}"
                    
                    # 轉換座標 (1-based -> 0-based)
                    pos_0based = pos_1based - 1
                    
                    # A. 抓取 Reference 序列 (包含 padding)
                    full_ref_seq = get_padded_sequence(seq_str, pos_0based, CONTEXT_RADIUS)
                    
                    # B. 建立 Alt 序列
                    # 因為確認過都是 SNV，我們直接把正中間那個鹼基換掉即可
                    # 驗證一下中間那個是不是真的等於 ref_base
                    mid_index = CONTEXT_RADIUS # 因為前面有 radius 個鹼基
                    
                    # 防呆檢查 (Debug 用)
                    extracted_center = full_ref_seq[mid_index]
                    if extracted_center != ref_base:
                        # 只有當不是 'N' 的時候才報錯 (因為如果剛好全是 N 就不準了)
                        if extracted_center != 'N':
                            print(f"  [警告] 序列中心不匹配! Key: {snv_key}, FASTA: {extracted_center}, VCF: {ref_base}")

                    # 組合 Alt 序列: 左邊 + ALT + 右邊
                    full_alt_seq = full_ref_seq[:mid_index] + alt_base + full_ref_seq[mid_index+1:]
                    
                    # 寫入檔案
                    line_items = [
                        snv_key,
                        chrom,
                        str(pos_1based),
                        ref_base,
                        alt_base,
                        feat_type,
                        mc_cat,
                        clnsig_cat,
                        full_ref_seq,
                        full_alt_seq
                    ]
                    f_out.write("\t".join(line_items) + "\n")
                    
                    processed_count += 1
                    
                except Exception as e:
                    print(f"  [錯誤] 處理 {snv_key} 時發生異常: {e}")
            
            # 釋放記憶體
            del seq_str

    print("\n" + "="*60)
    print(f"測試完成！已寫入 {processed_count} 筆資料。")
    print(f"請檢查輸出檔案: {OUTPUT_PATH}")
    print("="*60)

if __name__ == "__main__":
    main()