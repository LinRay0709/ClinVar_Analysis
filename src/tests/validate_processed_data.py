# src/tests/validate_processed_data.py
# 驗證 split_datasets.py 輸出的 CSV 檔案正確性
# 比對 ref_seq 和 alt_seq 是否與原始 genome FASTA 一致

import sys
import os
import pandas as pd

# --- 路徑設定 ---
current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.dirname(current_dir)
project_root = os.path.dirname(src_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

from src import config
from src.extract_seq.ref_seq_consistency_check import load_fasta_sequence

# ==========================================
# 使用者設定區
# ==========================================

# 選擇要檢查的檔案 (設為 True 啟用)
CHECK_TRAIN = False
CHECK_VAL = False
CHECK_TEST = True

# 要檢查哪一行 (0-based index，設為 None 則檢查全部)
TARGET_ROW = None

# 序列長度與中心點
SEQ_LENGTH = 1024
MID_IDX = 512

# ==========================================

def validate_row(row, chrom_seq):
    """
    驗證單一資料行
    回傳: (is_valid, error_message)
    """
    snv_key = row.name  # index 是 snv_key
    ref_seq = row['ref_seq']
    alt_seq = row['alt_seq']
    
    # 從 snv_key 解析資訊
    # 格式: chrom:pos:ref:alt
    parts = snv_key.split(":")
    pos_1based = int(parts[1])
    ref_base = parts[2]
    alt_base = parts[3]
    
    pos_0based = pos_1based - 1
    
    # 1. 從 genome 中取出預期的 ref_seq
    start = pos_0based - MID_IDX
    end = pos_0based + (SEQ_LENGTH - MID_IDX)
    
    # 處理邊界 (超出範圍補 @)
    left_pad = ""
    right_pad = ""
    real_start = start
    real_end = end
    
    if start < 0:
        left_pad = "@" * abs(start)
        real_start = 0
    if end > len(chrom_seq):
        right_pad = "@" * (end - len(chrom_seq))
        real_end = len(chrom_seq)
    
    expected_ref_seq = left_pad + chrom_seq[real_start:real_end] + right_pad
    
    # 2. 檢查 ref_seq 是否與 genome 一致
    if ref_seq != expected_ref_seq:
        # 找出差異位置
        diff_positions = [i for i in range(min(len(ref_seq), len(expected_ref_seq))) 
                          if ref_seq[i] != expected_ref_seq[i]]
        return False, f"ref_seq 與 genome 不一致 (差異位置: {diff_positions[:5]}...)"
    
    # 3. 檢查 alt_seq 是否等於 ref_seq 只把中心換成 alt_base
    expected_alt_seq = ref_seq[:MID_IDX] + alt_base + ref_seq[MID_IDX+1:]
    
    if alt_seq != expected_alt_seq:
        return False, f"alt_seq 不符合預期 (中心應為 '{alt_base}'，實際為 '{alt_seq[MID_IDX]}')"
    
    return True, "OK"

def validate_file(file_path, target_row=None):
    """驗證單一 CSV 檔案"""
    print(f"\n{'='*50}")
    print(f"驗證檔案: {os.path.basename(file_path)}")
    print(f"{'='*50}")
    
    if not os.path.exists(file_path):
        print(f"  [錯誤] 檔案不存在!")
        return
    
    df = pd.read_csv(file_path, index_col='snv_key')
    print(f"  -> 載入 {len(df)} 筆資料")
    
    # 決定要檢查的範圍
    if target_row is not None:
        if target_row >= len(df):
            print(f"  [錯誤] 指定行 {target_row} 超出範圍 (共 {len(df)} 行)")
            return
        df_to_check = df.iloc[[target_row]]
        print(f"  -> 只檢查第 {target_row} 行")
    else:
        df_to_check = df
        print(f"  -> 檢查全部 {len(df)} 行")
    
    # 按染色體分組處理 (減少 FASTA 載入次數)
    # 先從 snv_key 提取 chrom
    df_to_check = df_to_check.copy()
    df_to_check['chrom'] = df_to_check.index.map(lambda x: x.split(":")[0])
    
    errors = []
    
    for chrom, group in df_to_check.groupby('chrom'):
        print(f"\n  處理染色體 {chrom} ({len(group)} 筆)...")
        
        # 載入該染色體的序列
        chrom_seq = load_fasta_sequence(chrom)
        if chrom_seq is None:
            errors.append((None, f"chr{chrom}", f"無法載入染色體序列"))
            continue
        
        for idx, row in group.iterrows():
            is_valid, msg = validate_row(row, chrom_seq)
            if not is_valid:
                errors.append((idx, row.name, msg))
        
        # 釋放記憶體
        del chrom_seq
    
    # 輸出結果
    print(f"\n{'-'*40}")
    if errors:
        print(f"  [失敗] 發現 {len(errors)} 個錯誤:")
        for row_idx, snv_key, msg in errors[:10]:
            print(f"    {snv_key}: {msg}")
        if len(errors) > 10:
            print(f"    ... 還有 {len(errors) - 10} 個錯誤未顯示")
    else:
        print(f"  [成功] 全部通過驗證!")

def main():
    print("="*50)
    print("Processed Data 驗證程式 (使用 Genome FASTA)")
    print("="*50)
    print(f"參考序列目錄: {config.REF_DIR}")
    
    # 依設定檢查檔案
    if CHECK_TRAIN:
        validate_file(config.TRAIN_FILE, TARGET_ROW)
    
    if CHECK_VAL:
        validate_file(config.VAL_FILE, TARGET_ROW)
    
    if CHECK_TEST:
        validate_file(config.TEST_FILE, TARGET_ROW)
    
    print("\n" + "="*50)
    print("驗證完成")
    print("="*50)

if __name__ == "__main__":
    main()
