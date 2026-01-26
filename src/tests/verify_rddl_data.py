# verify_rddl_data.py
import os
import sys
import pandas as pd
import numpy as np
import pickle

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

import config as cfg

def check_init_file():
    """檢查 Python Package 標記是否存在"""
    init_path = os.path.join(cfg.RDDL_TASK_DIR, '__init__.py')
    print(f"\n[檢查環境設定]")
    if os.path.exists(init_path):
        print(f"✅ __init__.py 存在: {init_path}")
    else:
        print(f"❌ __init__.py 缺失！run_train.py 可能會無法載入模型。")

def check_manifest_file(filename):
    """
    檢查清單檔案格式並回傳資料筆數
    """
    print(f"\n[檢查清單檔案] {filename} ...")
    filepath = os.path.join(cfg.SPLIT_INFO_DIR, filename)
    
    if not os.path.exists(filepath):
        print(f"❌ 錯誤: 找不到檔案 {filepath}")
        return 0, None
    
    try:
        # [修改] 讀取 CSV (預設 header=0，自動讀取標題)
        df = pd.read_csv(filepath)
        
        # 1. 檢查標題欄位
        expected_cols = ['path', 'label', 'weight']
        if list(df.columns) == expected_cols:
            print(f"✅ 標題格式正確: {expected_cols}")
        else:
            print(f"❌ 標題格式錯誤！預期: {expected_cols}, 實際: {list(df.columns)}")
            return 0, None

        # 2. 檢查路徑格式 (是否為絕對路徑)
        first_path = df.iloc[0]['path']
        if os.path.isabs(first_path):
            print(f"✅ 路徑格式正確 (絕對路徑)。")
        else:
            print(f"⚠️ 警告: 路徑似乎是相對路徑 ({first_path})，可能會導致 RDDL 找不到檔案。")

        # 3. 統計數量
        count = len(df)
        pos_count = len(df[df['label'] == 1])
        neg_count = len(df[df['label'] == 0])
        print(f"  -> 總筆數: {count} (Pos: {pos_count}, Neg: {neg_count})")
            
        return count, first_path # 回傳筆數和第一筆路徑供檢查
        
    except Exception as e:
        print(f"❌ 讀取失敗: {e}")
        return 0, None

def check_pickle_content(abs_path):
    """檢查 Pickle 內容 (使用絕對路徑)"""
    print(f"\n[檢查 Pickle 資料內容] ...")
    print(f"目標檔案: {abs_path}")
    
    if not os.path.exists(abs_path):
        print(f"❌ 錯誤: 找不到 Pickle 檔案！")
        return

    try:
        with open(abs_path, 'rb') as f:
            data = pickle.load(f)
            
        print(f"✅ 載入成功！")
        print(f"  - 資料形狀 (Shape): {data.shape}")
        
        # 檢查 Ref/Alt 是否都存在 (1024, 8)
        if data.shape == (1024, 8):
            print(f"✅ [通過] 形狀符合 Ref+Alt 設計 (1024, 8)")
        elif data.shape == (1024, 4):
            print(f"⚠️ [注意] 形狀為 (1024, 4)，只有單一序列。")
        else:
            print(f"⚠️ [警告] 形狀異常: {data.shape}")

        # 檢查數值 (新規則: 0, 0.25, 1 都是合法的)
        # 0: @字元 或 0位置的 one-hot
        # 0.25: N字元 (均勻分布)
        # 1: 標準鹼基 ACGT 的 one-hot
        unique_vals = np.unique(data)
        valid_vals = {0.0, 0.25, 1.0}
        
        # 使用 np.isclose 處理浮點數精度
        all_valid = all(any(np.isclose(v, valid) for valid in valid_vals) for v in unique_vals)
        
        if all_valid:
            print(f"✅ [通過] 數值符合編碼規則。")
            print(f"  - 出現的數值: {sorted(unique_vals)}")
            
            # 統計特殊字元
            count_025 = np.sum(np.isclose(data, 0.25))
            count_zeros_row = np.sum(np.all(data == 0, axis=1))  # 整行都是0 (對應 @)
            
            if count_025 > 0:
                print(f"  - 包含 N 編碼 (0.25): 共 {count_025} 個數值")
            if count_zeros_row > 0:
                # 因為 ref+alt 合併，所以檢查前4欄和後4欄
                ref_zeros = np.sum(np.all(data[:, :4] == 0, axis=1))
                alt_zeros = np.sum(np.all(data[:, 4:] == 0, axis=1))
                print(f"  - 包含 @ 編碼 (全零行): Ref 有 {ref_zeros} 個, Alt 有 {alt_zeros} 個")
        else:
            unexpected = [v for v in unique_vals if not any(np.isclose(v, valid) for valid in valid_vals)]
            print(f"❌ [失敗] 發現非法數值: {unexpected}")
            
    except Exception as e:
        print(f"❌ Pickle 解析失敗: {e}")

def main():
    print("="*60)
    print("RDDL 資料集健康檢查 (Upgrade Version)")
    print("="*60)
    
    # 1. 檢查 __init__.py
    check_init_file()

    # 2. 檢查各資料集並統計數量
    total_train = 0
    sample_pkl_path = None

    # 檢查 5 個 Training Folds
    for i in range(1, 6):
        count, path = check_manifest_file(f'training_{i}.csv')
        total_train += count
        if path and sample_pkl_path is None:
            sample_pkl_path = path

    # 檢查 Ensemble (Validation)
    ens_count, _ = check_manifest_file('ensemble.csv')
    
    # 檢查 Test
    test_count, _ = check_manifest_file('test.csv')

    print("\n" + "="*60)
    print("📊 資料數量總結 Report")
    print("-" * 30)
    print(f"Training Set (Total 5 Folds): {total_train} 筆")
    print(f"Ensemble Set (Validation)   : {ens_count} 筆")
    print(f"Test Set                    : {test_count} 筆")
    print("-" * 30)
    print(f"總資料量                    : {total_train + ens_count + test_count} 筆")
    print("="*60)

    # 3. 抽樣檢查一個 Pickle 檔
    if sample_pkl_path:
        check_pickle_content(sample_pkl_path)
    
    print("\n檢查結束。")

if __name__ == "__main__":
    main()