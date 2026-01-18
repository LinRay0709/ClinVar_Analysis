# src/split_datasets.py
import sys
import os
import random
import pandas as pd
import re

# --- 路徑設定 ---
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

from src import config

# ==========================================
# 參數設定
# ==========================================
SEED = 42  # 固定種子，確保每次切分結果一樣
SPLIT_RATIOS = (0.8, 0.1, 0.1)  # Train / Ensemble / Test

def parse_clstr_file(clstr_path):
    """
    解析 CD-HIT .clstr 檔案
    邏輯:
    1. 忽略非代表序列 (沒有 '*' 的行) -> 降低組內相似度
    2. 提取代表序列的 snv_key 和 label
    
    Returns:
        pos_keys: list of snv_keys (Pathogenic)
        neg_keys: list of snv_keys (Benign)
    """
    print(f"正在解析 Cluster 檔案: {clstr_path} ...")
    
    pos_keys = []
    neg_keys = []
    
    with open(clstr_path, 'r') as f:
        for line in f:
            # 範例行: 0	257nt, >chr1:12345:A:G|1... *
            if "*" in line:  # 只處理代表序列 (Representative)
                # 使用 Regex 提取 > 與 ... 之間的內容
                # 預期格式: >snv_key|label...
                # re: ()內為ㄧgroup
                match = re.search(r'>(.+)\|([01])', line)
                if match:
                    snv_key = match.group(1)
                    label = int(match.group(2))
                    
                    if label == 1:
                        pos_keys.append(snv_key)
                    else:
                        neg_keys.append(snv_key)
    
    return pos_keys, neg_keys

def split_data(keys, ratios):
    """
    將列表依照比例切分為三份
    """
    random.shuffle(keys)
    n = len(keys)
    
    n_train = int(n * ratios[0])
    n_val = int(n * ratios[1])
    # 剩下的都給 test，確保總數不變
    
    train_keys = keys[:n_train]
    val_keys = keys[n_train : n_train + n_val]
    test_keys = keys[n_train + n_val:]
    
    return train_keys, val_keys, test_keys

def load_and_reconstruct_from_fasta(fasta_path, valid_keys_set):
    """
    從 FASTA 讀取 alt_seq，並根據 snv_key 推導出 ref_seq
    
    Args:
        fasta_path: for_cdhit.fasta 的路徑
        valid_keys_set: 從 CD-HIT 篩選出來的黃金 Key 清單 (Set)
    
    Returns:
        pd.DataFrame: 包含 snv_key, label, alt_seq, ref_seq 的表格
    """
    print(f"正在從 FASTA 讀取並重建序列資料...")
    data_list = []
    
    with open(fasta_path, 'r') as f:
        current_key = None
        current_label = None
        
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                # 解析 Header: >chr1:12345:A:G|1
                content = line[1:] # 移除 '>'
                if "|" in content:
                    key_part, label_part = content.split("|")
                    
                    # 只有當這個 Key 在我們的黃金清單內，才處理它
                    if key_part in valid_keys_set:
                        current_key = key_part
                        current_label = int(label_part)
                    else:
                        current_key = None # 忽略這個 Cluster 的非代表序列
                else:
                    current_key = None
            else:
                # 這是序列行 (alt_seq)
                if current_key:
                    alt_seq = line
                    
                    # --- 核心邏輯：重建 Ref Seq ---
                    # 1. 解析 Key 取得 Reference Base
                    # Key 格式: chrom:pos:ref:alt
                    parts = current_key.split(":")
                    ref_base = parts[2]
                    alt_base = parts[3]
                    
                    # 2. 驗證中心點 (Sanity Check)
                    mid_idx = 128  # 257bp 的中心點是 index 128
                    if alt_seq[mid_idx] != alt_base:
                        # 只有當不是 N 時才報警 (有些資料 alt 可能是 N)
                        if alt_seq[mid_idx] != 'N':
                            print(f"[警告] Key {current_key} 的 Alt ({alt_base}) 與序列中心 ({alt_seq[mid_idx]}) 不符，跳過。")
                            current_key = None
                            continue

                    # 3. 生成 Ref Seq (把中間換回 ref_base)
                    ref_seq = alt_seq[:mid_idx] + ref_base + alt_seq[mid_idx+1:]
                    
                    # 4. 加入資料列表
                    data_list.append({
                        "snv_key": current_key,
                        "label": current_label,
                        "ref_seq": ref_seq, # 重建的 Ref
                        "alt_seq": alt_seq  # 原始的 Alt
                    })
                    
                    # 重置以免重複讀取
                    current_key = None

    df = pd.DataFrame(data_list)
    df.set_index("snv_key", inplace=True)
    return df

def main():
    print("="*60)
    print("資料集拆分程式 (Cluster-based Stratified Split)")
    print("="*60)
    
    # 0. 環境準備
    os.makedirs(config.PROCESSED_DIR, exist_ok=True)
    random.seed(SEED)

    # 1. 解析 CD-HIT 結果 (只取代表序列)
    if not os.path.exists(config.CDHIT_CLSTR_FILE):
        print(f"錯誤: 找不到 {config.CDHIT_CLSTR_FILE}")
        return

    pos_keys, neg_keys = parse_clstr_file(config.CDHIT_CLSTR_FILE)
    
    all_valid_keys = set(pos_keys + neg_keys)

    print(f"\n[去冗餘後統計]")
    print(f"  - Positive (Pathogenic) 代表序列: {len(pos_keys)}")
    print(f"  - Negative (Benign)     代表序列: {len(neg_keys)}")
    print(f"  - 總計: {len(pos_keys) + len(neg_keys)}")

    # 2. 分層拆分 (Stratified Split)
    print("\n[執行 80/10/10 拆分]")
    
    # 拆分 Positive
    p_train, p_val, p_test = split_data(pos_keys, SPLIT_RATIOS)
    # 拆分 Negative
    n_train, n_val, n_test = split_data(neg_keys, SPLIT_RATIOS)
    
    # 合併
    train_keys = set(p_train + n_train)
    val_keys = set(p_val + n_val)
    test_keys = set(p_test + n_test)
    
    print(f"  - Train Set : {len(train_keys)} (Pos:{len(p_train)}, Neg:{len(n_train)})")
    print(f"  - Val Set   : {len(val_keys)} (Pos:{len(p_val)}, Neg:{len(n_val)})")
    print(f"  - Test Set  : {len(test_keys)} (Pos:{len(p_test)}, Neg:{len(n_test)})")

    # 3. 讀取原始資料並寫入檔案
    df = load_and_reconstruct_from_fasta(config.CDHIT_FASTA_FILE, all_valid_keys)
    
    print(f"  -> 重建完成，資料庫中可用筆數: {len(df)}")
    
    # 檢查是否所有 Key 都有找到
    if len(df) != len(all_valid_keys):
        print(f"  [警告] Key 數量不匹配！預期 {len(all_valid_keys)}，實際重建 {len(df)}。")

    # 4. 存檔 (邏輯不變，但現在 df 來源變了)
    def save_subset(keys_set, output_path, label_name):
        valid_keys = list(keys_set.intersection(df.index))
        subset_df = df.loc[valid_keys].copy()
        subset_df.to_csv(output_path, index=True)
        print(f"  -> 已儲存 {label_name}: {output_path} ({len(subset_df)} 筆)")

    save_subset(train_keys, config.TRAIN_FILE, "Training Set")
    save_subset(val_keys, config.VAL_FILE, "Validation Set")
    save_subset(test_keys, config.TEST_FILE, "Test Set")

    print("\n" + "="*60)
    print("✅ 資料集製作完成 (FASTA-Reconstruction Mode)")
    print("資料來源: for_cdhit.fasta (alt_seq)")
    print("處理邏輯: 自動重建 ref_seq，確保僅中心點不同")
    print("注意: 此模式生成的檔案不包含 feature_type 等 Metadata")
    print("="*60)

if __name__ == "__main__":
    main()