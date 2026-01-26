# src/split_datasets.py
import sys
import os
import random
import pandas as pd
import re
from collections import defaultdict

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
SEQ_LENGTH = 1024  # 序列長度
MID_IDX = 512      # 中心點 index (0-based)

def parse_clstr_file_full(clstr_path):
    """
    解析 CD-HIT .clstr 檔案，回傳完整的 Cluster 資料
    
    Returns:
        clusters: dict[cluster_id] -> list of {"snv_key": str, "label": int}
        total_entries: int - 原始資料總筆數（未處理前）
    """
    print(f"正在解析 Cluster 檔案: {clstr_path} ...")
    
    clusters = defaultdict(list)
    current_cluster_id = None
    total_entries = 0  # 計算原始資料總筆數
    
    with open(clstr_path, 'r') as f:
        for line in f:
            line = line.strip()
            
            # 新 Cluster 開始
            if line.startswith(">Cluster"):
                # 格式: >Cluster 0
                current_cluster_id = int(line.split()[1])
            else:
                # 成員序列行
                # 範例: 0	1024nt, >chr1:12345:A:G|1... *
                # 或:   1	1024nt, >chr1:12346:C:T|0... at 95.5%
                match = re.search(r'>(.+)\|([01])', line)
                if match and current_cluster_id is not None:
                    snv_key = match.group(1)
                    label = int(match.group(2))
                    clusters[current_cluster_id].append({
                        "snv_key": snv_key,
                        "label": label
                    })
                    total_entries += 1  # 計算每筆資料
    
    print(f"  -> 解析完成，共 {len(clusters)} 個 Clusters")
    print(f"  -> 原始資料總筆數: {total_entries}")
    return clusters, total_entries

def get_mutation_type(snv_key):
    """
    從 snv_key 提取突變類型
    snv_key 格式: chrom:pos:ref:alt
    回傳: "ref->alt" (例如 "A->G")
    """
    parts = snv_key.split(":")
    ref = parts[2]
    alt = parts[3]
    return f"{ref}->{alt}"

def select_from_clusters(clusters):
    """
    根據新邏輯從 Clusters 中選取資料
    
    邏輯:
    1. 如果該 cluster 中 label 都是 1或0:
       -> 分突變類型，每種類型只選一筆
    2. 如果該 cluster 的 label 有 1 有 0:
       -> 分突變類型
       -> 如果同突變類型中 label 有 1 有 0，這些資料都不要
       -> 如果同突變類型中都是 1 或都是 0，每種類型只選一筆
    
    Returns:
        pos_keys: list of snv_keys (Pathogenic, label=1)
        neg_keys: list of snv_keys (Benign, label=0)
    """
    pos_keys = []
    neg_keys = []
    
    stats = {
        "pure_clusters": 0,
        "mixed_clusters": 0,
        "discarded_mutation_types": 0
    }
    
    for cluster_id, members in clusters.items():
        # 取得該 cluster 的所有 labels
        labels = set(m["label"] for m in members)
        
        # 按突變類型分組
        mutation_groups = defaultdict(list)
        for m in members:
            mut_type = get_mutation_type(m["snv_key"])
            mutation_groups[mut_type].append(m)
        
        if labels == {1}:
            # Case 1: 純 Pathogenic cluster
            stats["pure_clusters"] += 1
            for mut_type, group in mutation_groups.items():
                # 每種類型選一筆
                pos_keys.append(group[0]["snv_key"])
        
        elif labels == {0}:
            # Case 1b: 純 Benign cluster (邏輯相同)
            stats["pure_clusters"] += 1
            for mut_type, group in mutation_groups.items():
                neg_keys.append(group[0]["snv_key"])
        
        else:
            # Case 2: 混合 cluster (有 0 也有 1)
            stats["mixed_clusters"] += 1
            for mut_type, group in mutation_groups.items():
                group_labels = set(m["label"] for m in group)
                
                if len(group_labels) > 1:
                    # 同突變類型有衝突 -> 全部丟棄
                    stats["discarded_mutation_types"] += 1
                    continue
                else:
                    # 同突變類型 label 一致 -> 選一筆
                    selected = group[0]
                    if selected["label"] == 1:
                        pos_keys.append(selected["snv_key"])
                    else:
                        neg_keys.append(selected["snv_key"])
    
    print(f"\n[Cluster 選取統計]")
    print(f"  - 純 Label Clusters: {stats['pure_clusters']}")
    print(f"  - 混合 Label Clusters: {stats['mixed_clusters']}")
    print(f"  - 因衝突而丟棄的突變類型數: {stats['discarded_mutation_types']}")
    
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
    同時偵測並輸出重複的 snv_key 到 redundant_snvs.tsv
    
    Args:
        fasta_path: for_cdhit.fasta 的路徑
        valid_keys_set: 從 CD-HIT 篩選出來的黃金 Key 清單 (Set)
    
    Returns:
        pd.DataFrame: 包含 snv_key, label, alt_seq, ref_seq 的表格（已去重）
    """
    print(f"正在從 FASTA 讀取並重建序列資料...")
    
    # 先收集所有資料（包含重複的）
    all_data = defaultdict(list)  # snv_key -> list of records
    
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
                        current_key = None # 忽略不在清單內的序列
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
                    # 序列長度 1024，中心點是 index 512
                    if alt_seq[MID_IDX] != alt_base:
                        # 只有當不是 N 或 @ 時才報警
                        if alt_seq[MID_IDX] not in ('N', '@'):
                            print(f"[警告] Key {current_key} 的 Alt ({alt_base}) 與序列中心 ({alt_seq[MID_IDX]}) 不符，跳過。")
                            current_key = None
                            continue

                    # 3. 生成 Ref Seq (把中間換回 ref_base)
                    ref_seq = alt_seq[:MID_IDX] + ref_base + alt_seq[MID_IDX+1:]
                    
                    # 4. 加入資料字典
                    all_data[current_key].append({
                        "snv_key": current_key,
                        "label": current_label,
                        "ref_seq": ref_seq,
                        "alt_seq": alt_seq
                    })
                    
                    # 重置以免重複讀取
                    current_key = None

    # 偵測重複的 snv_key 並輸出
    redundant_records = []
    final_data_list = []
    
    for snv_key, records in all_data.items():
        if len(records) > 1:
            # 有重複，全部加入 redundant
            redundant_records.extend(records)
            # 只取第一筆加入最終資料
            final_data_list.append(records[0])
        else:
            # 沒有重複，直接加入
            final_data_list.append(records[0])
    
    # 輸出重複資料到 tsv
    if redundant_records:
        redundant_path = os.path.join(config.OUTPUT_DIR, "redundant_snvs.tsv")
        redundant_df = pd.DataFrame(redundant_records)
        redundant_df.to_csv(redundant_path, sep='\t', index=False)
        print(f"  [警告] 發現 {len(redundant_records)} 筆重複 snv_key 資料，已輸出至: {redundant_path}")
        print(f"  -> 涉及 {len([k for k, v in all_data.items() if len(v) > 1])} 個不同的 snv_key")
    else:
        print(f"  [資訊] 沒有發現重複的 snv_key。")

    df = pd.DataFrame(final_data_list)
    df.set_index("snv_key", inplace=True)
    return df

def main():
    print("="*60)
    print("資料集拆分程式 (Mutation-Type Based Selection)")
    print("="*60)
    print(f"序列長度: {SEQ_LENGTH} bp, 中心點 index: {MID_IDX}")
    
    # 0. 環境準備
    os.makedirs(config.PROCESSED_DIR, exist_ok=True)
    random.seed(SEED)

    # 1. 解析 CD-HIT 結果 (取得完整 Cluster 資料)
    if not os.path.exists(config.CDHIT_CLSTR_FILE):
        print(f"錯誤: 找不到 {config.CDHIT_CLSTR_FILE}")
        return

    clusters, total_raw_entries = parse_clstr_file_full(config.CDHIT_CLSTR_FILE)
    
    # 2. 根據新邏輯選取資料
    pos_keys, neg_keys = select_from_clusters(clusters)
    
    print(f"\n[原始 vs 篩選後]")
    print(f"  - CD-HIT .clstr 原始資料總筆數: {total_raw_entries}")
    
    all_valid_keys = set(pos_keys + neg_keys)

    print(f"\n[篩選後統計]")
    print(f"  - Positive (Pathogenic): {len(pos_keys)}")
    print(f"  - Negative (Benign)    : {len(neg_keys)}")
    print(f"  - 總計: {len(pos_keys) + len(neg_keys)}")

    # 3. 分層拆分 (Stratified Split)
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

    # 4. 讀取原始資料並寫入檔案
    df = load_and_reconstruct_from_fasta(config.CDHIT_FASTA_FILE, all_valid_keys)
    
    print(f"  -> 重建完成，資料庫中可用筆數: {len(df)}")
    
    # 檢查是否所有 Key 都有找到
    if len(df) != len(all_valid_keys):
        print(f"  [警告] Key 數量不匹配！預期 {len(all_valid_keys)}，實際重建 {len(df)}。")

    # 5. 存檔
    def save_subset(keys_set, output_path, label_name):
        valid_keys = list(keys_set.intersection(df.index))
        subset_df = df.loc[valid_keys].copy()
        subset_df.to_csv(output_path, index=True)
        print(f"  -> 已儲存 {label_name}: {output_path} ({len(subset_df)} 筆)")

    save_subset(train_keys, config.TRAIN_FILE, "Training Set")
    save_subset(val_keys, config.VAL_FILE, "Validation Set")
    save_subset(test_keys, config.TEST_FILE, "Test Set")

    print("\n" + "="*60)
    print("✅ 資料集製作完成 (Mutation-Type Based Selection Mode)")
    print(f"序列長度: {SEQ_LENGTH} bp")
    print("處理邏輯:")
    print("  1. 純 label cluster -> 每種突變類型選一筆")
    print("  2. 混合 label cluster -> 衝突突變類型丟棄，其餘每種選一筆")
    print("="*60)

if __name__ == "__main__":
    main()