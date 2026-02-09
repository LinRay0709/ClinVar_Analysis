# bridge_utils.py
import os
import pickle
import numpy as np
import pandas as pd
import config as cfg  # 引用設定檔

def one_hot_encode(sequence):
    """
    將 DNA 字串轉換為 One-Hot 編碼 (Numpy Array)
    格式: (Length, 4) -> 專為 Keras/TensorFlow 設計
    
    特殊處理:
    - A, C, G, T: 標準 One-Hot (1 在對應位置)
    - N: 均勻分布 [0.25, 0.25, 0.25, 0.25]
    - @: 全零 [0, 0, 0, 0] (表示超出基因體範圍)
    """
    seq_len = len(sequence)
    one_hot = np.zeros((seq_len, 4), dtype=np.float32)
    
    for i, base in enumerate(sequence):
        base = base.upper()
        if base in cfg.DNA_MAP:
            # 標準鹼基: A, C, G, T
            idx = cfg.DNA_MAP[base]
            one_hot[i, idx] = 1.0
        elif base == 'N':
            # N: 不確定鹼基，用均勻分布表示
            one_hot[i, :] = 0.25
        elif base == '@':
            # @: 超出基因體範圍，保持全零
            pass  # 已經是 0
        else:
            print(f"[警告] 遇到未知的鹼基: {base}")
    
    return one_hot

def save_pickle_and_get_path(row, output_dir):
    """
    將單筆資料轉成矩陣並存檔，回傳 RDDL 需要的相對路徑
    同時編碼 Ref 和 Alt，並以 (Length, 8) 的格式儲存
    """
    # 處理檔名
    if 'snv_key' in row:
        safe_key = str(row['snv_key']).replace(':', '_')
    else:
        safe_key = f"sample_{row.name}"
        
    filename = f"{safe_key}.pkl"
    filepath = os.path.join(output_dir, filename)
    
    # 編碼與存檔
    ref_seq = row['ref_seq']
    alt_seq = row['alt_seq']

    ref_tensor = one_hot_encode(ref_seq)  # Shape: (1024, 4)
    alt_tensor = one_hot_encode(alt_seq)  # Shape: (1024, 4)

    # axis=1 代表在「頻道」方向合併 -> (1024, 4+4) = (1024, 8)
    combined_tensor = np.concatenate([ref_tensor, alt_tensor], axis=1)
    
    # 'wb': 二進位模式
    with open(filepath, 'wb') as f:
        pickle.dump(combined_tensor, f)
        
    # 回傳相對路徑
    return filepath

def process_dataset(csv_filename, dataset_type):
    """
    讀取 CSV -> 轉檔 -> 回傳清單 DataFrame
    """
    csv_path = os.path.join(cfg.CSV_SOURCE_PATH, csv_filename)
    
    if not os.path.exists(csv_path):
        print(f"[警告] 找不到檔案: {csv_path}")
        return None
        
    print(f"正在處理 {dataset_type} ({csv_filename})...")
    
    df = pd.read_csv(csv_path)
    
    # 檢查檔案是否為空
    if len(df) == 0:
        print(f"[警告] 檔案為空: {csv_path}")
        return None
    
    # 檢查必要欄位
    required_cols = ['alt_seq', 'ref_seq', 'label']
    for col in required_cols:
        if col not in df.columns:
             raise ValueError(f"{csv_filename} 缺少必要欄位: {col}")
        
    file_list = []
    
    for idx, row in df.iterrows():
        label = int(row['label'])
        
        # 決定存到 USER_pos 還是 USER_neg
        if label == 1:
            target_dir = cfg.POS_DIR
            rel_dir = "USER_pos"
        else:
            target_dir = cfg.NEG_DIR
            rel_dir = "USER_neg"
            
        # 執行存檔
        abs_path = save_pickle_and_get_path(row, target_dir)

        # 配合RDDL要求加上weight
        weight = 1.0

        file_list.append([abs_path, label, weight])
        
    return pd.DataFrame(file_list, columns=['sample', 'label', 'weight'])

def stratified_split(df, n_splits=5, random_state=42):
    """
    將 DataFrame 依標籤比例分層切分為 n 份
    確保每個 fold 的正負樣本比例與原始資料一致
    """
    from sklearn.model_selection import StratifiedKFold
    
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    
    X = df.index
    y = df['label'].values
    
    folds = []
    for _, fold_idx in skf.split(X, y):
        folds.append(df.iloc[fold_idx].copy())
    
    return folds