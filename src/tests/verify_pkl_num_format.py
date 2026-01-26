# src/tests/verify_pkl_num_format.py
# 遍歷所有 pickle 檔案，檢查數量、形狀、數值是否正確

import os
import sys
import pickle
import numpy as np
import glob

# --- 路徑設定 ---
current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.dirname(current_dir)
project_root = os.path.dirname(src_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

from src import config as cfg

# ==========================================
# 預期參數
# ==========================================
EXPECTED_SHAPE = (1024, 8)  # 預期形狀
VALID_VALUES = {0.0, 0.25, 1.0}  # 合法數值

def check_pickle_file(filepath):
    """
    檢查單一 pickle 檔案
    回傳: (is_valid, error_message)
    """
    try:
        with open(filepath, 'rb') as f:
            data = pickle.load(f)
        
        # 1. 檢查形狀
        if data.shape != EXPECTED_SHAPE:
            return False, f"形狀錯誤: {data.shape} (預期 {EXPECTED_SHAPE})"
        
        # 2. 檢查數值
        unique_vals = np.unique(data)
        invalid_vals = [v for v in unique_vals if not any(np.isclose(v, valid) for valid in VALID_VALUES)]
        
        if invalid_vals:
            return False, f"包含非法數值: {invalid_vals}"
        
        return True, "OK"
        
    except Exception as e:
        return False, f"讀取失敗: {e}"

def main():
    print("="*60)
    print("Pickle 檔案完整性驗證")
    print("="*60)
    print(f"預期形狀: {EXPECTED_SHAPE}")
    print(f"合法數值: {VALID_VALUES}")
    print("="*60)
    
    # 收集所有 pickle 檔案
    pos_dir = cfg.POS_DIR
    neg_dir = cfg.NEG_DIR
    
    pos_files = glob.glob(os.path.join(pos_dir, "*.pkl"))
    neg_files = glob.glob(os.path.join(neg_dir, "*.pkl"))
    
    all_files = pos_files + neg_files
    total_count = len(all_files)
    
    print(f"\n[檔案數量統計]")
    print(f"  - USER_pos: {len(pos_files)} 個檔案")
    print(f"  - USER_neg: {len(neg_files)} 個檔案")
    print(f"  - 總計: {total_count} 個檔案")
    
    if total_count == 0:
        print("\n[警告] 沒有找到任何 pickle 檔案！")
        print(f"  搜尋目錄: {pos_dir}")
        print(f"            {neg_dir}")
        return
    
    # 遍歷檢查
    print(f"\n[開始驗證] 正在檢查 {total_count} 個檔案...")
    
    errors = []
    valid_count = 0
    
    for i, filepath in enumerate(all_files):
        is_valid, msg = check_pickle_file(filepath)
        
        if is_valid:
            valid_count += 1
        else:
            errors.append((os.path.basename(filepath), msg))
        
        # 進度顯示 (每 500 個檔案顯示一次)
        if (i + 1) % 500 == 0:
            print(f"  已檢查: {i + 1}/{total_count} ({(i+1)/total_count*100:.1f}%)")
    
    # 輸出結果
    print("\n" + "="*60)
    print("驗證結果")
    print("="*60)
    print(f"總檔案數: {total_count}")
    print(f"通過驗證: {valid_count}")
    print(f"驗證失敗: {len(errors)}")
    
    if errors:
        print(f"\n[錯誤清單] (前 20 筆):")
        for filename, msg in errors[:20]:
            print(f"  ❌ {filename}: {msg}")
        if len(errors) > 20:
            print(f"  ... 還有 {len(errors) - 20} 個錯誤未顯示")
        print("\n結果: ❌ 驗證失敗")
    else:
        print("\n結果: ✅ 全部通過驗證！")

if __name__ == "__main__":
    main()
