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

def check_manifest_file(filename):
    print(f"\n[檢查清單檔案] {filename} ...")
    filepath = os.path.join(cfg.SPLIT_INFO_DIR, filename)
    
    if not os.path.exists(filepath):
        print(f"❌ 錯誤: 找不到檔案 {filepath}")
        return None
    
    # 嘗試讀取 (RDDL 要求無 Header)
    # 我們讀取前兩行來看看
    try:
        df = pd.read_csv(filepath, header=None, nrows=5)
        print(f"✅ 讀取成功。前 5 筆資料預覽:")
        print(df)
        
        # 檢查欄位數 (應該要有 2 欄: 路徑, Label)
        if df.shape[1] != 2:
            print(f"⚠️ 警告: 欄位數量不正確！預期 2 欄，實際 {df.shape[1]} 欄。")
        else:
            print(f"✅ 格式正確 (2 欄)。")
            
        return df.iloc[0, 0] # 回傳第一筆資料的路徑，供後續檢查 pickle 用
        
    except Exception as e:
        print(f"❌ 讀取失敗: {e}")
        return None

def check_pickle_content(rel_path):
    print(f"\n[檢查 Pickle 資料內容] ...")
    # RDDL 的路徑是相對路徑 (例如 USER_pos/xxx.pkl)
    # 我們需要加上 RDDL_TASK_DIR
    full_path = os.path.join(cfg.RDDL_TASK_DIR, rel_path)
    
    print(f"目標檔案: {full_path}")
    
    if not os.path.exists(full_path):
        print(f"❌ 錯誤: 找不到 Pickle 檔案！請檢查 bridge_to_rddl.py 的路徑設定。")
        return

    try:
        with open(full_path, 'rb') as f:
            data = pickle.load(f)
            
        print(f"✅ 載入成功！")
        print(f"  - 資料型態 (Type): {type(data)}")
        print(f"  - 資料形狀 (Shape): {data.shape}")
        print(f"  - 數據類型 (Dtype): {data.dtype}")
        
        # === 關鍵檢查區 ===
        
        # 1. 檢查 Ref/Alt 是否都存在
        # 如果是 (257, 8)，代表 Ref(4) + Alt(4) 成功合併
        expected_shape = (257, 8) 
        if data.shape == expected_shape:
            print(f"✅ [通過] 形狀符合 Ref+Alt 設計: {data.shape} (257 bp, 8 channels)")
        elif data.shape == (257, 4):
            print(f"⚠️ [注意] 形狀為 (257, 4)，代表只包含了單一序列 (只有 Alt?)。")
        else:
            print(f"⚠️ [警告] 形狀異常: {data.shape}，請確認是否符合您的模型輸入層設計。")

        # 2. 檢查數值是否為 One-Hot (只有 0 和 1)
        unique_vals = np.unique(data)
        if np.all(np.isin(unique_vals, [0, 1])):
            print(f"✅ [通過] 數值檢查正常 (One-Hot Encoded)。包含數值: {unique_vals}")
        else:
            print(f"❌ [失敗] 發現非 0/1 的數值: {unique_vals}")
            
    except Exception as e:
        print(f"❌ Pickle 解析失敗: {e}")

def main():
    print("="*60)
    print("RDDL 資料集健康檢查 (Sanity Check)")
    print("="*60)
    
    # 1. 檢查 RDDL_splitting_info 裡的 training_1.csv
    # 這是 RDDL 開始訓練時讀的第一個檔案
    first_pkl_path = check_manifest_file('training_1.csv')
    
    if first_pkl_path:
        # 2. 順藤摸瓜，去檢查該路徑指向的 Pickle 檔
        check_pickle_content(first_pkl_path)
    
    print("\n" + "="*60)
    print("檢查結束。請根據上述結果確認是否符合預期。")
    print("="*60)

if __name__ == "__main__":
    main()