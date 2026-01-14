# tests/verify_consistency.py
import sys
import os
import pandas as pd

# --- 路徑設定 (讓測試程式能找到 src 裡的 config) ---
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

import config
# -----------------------------------------------

def main():
    print("="*50)
    print("邏輯一致性驗證工具")
    print("="*50)
    print(f"正在讀取檔案: {config.OUTPUT_FILE} ...")
    
    if not os.path.exists(config.OUTPUT_FILE):
        print(f"錯誤: 找不到檔案 {config.OUTPUT_FILE}")
        print("請先執行 main.py 產生資料。")
        return

    # 讀取 TSV 檔案
    # header=0 表示第一列是標題
    try:
        df = pd.read_csv(config.OUTPUT_FILE, sep='\t')
    except Exception as e:
        print(f"讀取失敗: {e}")
        return

    total_rows = len(df)
    print(f"總資料筆數: {total_rows}")
    print("-" * 30)

    # 1. 計算第六行 (feature_type) 為 'CDS' 的數量
    # 對應欄位名稱: feature_type
    cds_count = (df['feature_type'] == 'CDS').sum()
    
    # 2. 計算第八行 (clnsig_category) 為 'benign' 的數量
    # 對應欄位名稱: clnsig_category
    benign_count = (df['clnsig_category'] == 'benign').sum()
    
    # 3. 計算第九行 (mc_category) 為 'missense_variant' 的數量
    # 對應欄位名稱: mc_category
    missense_count = (df['mc_category'] == 'missense_variant').sum()

    # 輸出結果
    print(f"1. CDS 區域變異數量 (Feature Type):")
    print(f"   {cds_count:>8}  ({cds_count/total_rows*100:.2f}%)")
    
    print(f"\n2. Benign 分類數量 (Clinical Significance):")
    print(f"   {benign_count:>8}  ({benign_count/total_rows*100:.2f}%)")
    
    print(f"\n3. Missense Variant 數量 (Molecular Consequence):")
    print(f"   {missense_count:>8}  ({missense_count/total_rows*100:.2f}%)")
    print("="*50)

    # 加碼檢查：列出所有欄位的名稱，確認「第六行、第八行」是否對應正確
    print("\n[欄位索引檢查]")
    for idx, col in enumerate(df.columns):
        # idx 是 0-based，所以顯示時 +1 變成人類習慣的 1-based
        print(f"第 {idx+1} 行: {col}")

if __name__ == "__main__":
    main()