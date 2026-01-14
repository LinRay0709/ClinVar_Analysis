# tests/check_multiallelic.py
import sys
import os
import pandas as pd

# --- 1. 設定路徑以匯入 src.config ---
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

import config
# ----------------------------------

def check_multiallelic():
    print("="*60)
    print("第 4 行 (ALT) 多重等位基因檢查工具")
    print("="*60)
    
    # 檢查檔案是否存在
    if not os.path.exists(config.OUTPUT_FILE):
        print(f"錯誤: 找不到輸出檔案 {config.OUTPUT_FILE}")
        return

    # 讀取 TSV
    print(f"正在讀取: {config.OUTPUT_FILE} ...")
    df = pd.read_csv(config.OUTPUT_FILE, sep='\t', dtype=str) 
    # dtype=str 確保所有欄位都當作字串讀取，避免 pos 被轉成數字

    total_count = len(df)
    
    # --- 核心檢查邏輯 ---
    # 檢查 'alt' 欄位 (第 4 行) 是否包含逗號 ','
    # na=False 是為了防止欄位如果是空值(NaN)時報錯，雖然理論上不該有空值
    multi_allelic_mask = df['alt'].str.contains(',', na=False)
    
    multi_count = multi_allelic_mask.sum()
    single_count = total_count - multi_count
    
    # --- 輸出統計結果 ---
    print("\n[統計結果]")
    print(f"總變異筆數 : {total_count}")
    print("-" * 30)
    print(f"單一 ALT   : {single_count:>8} 筆 ({single_count/total_count*100:.2f}%)")
    print(f"多重 ALT   : {multi_count:>8} 筆 ({multi_count/total_count*100:.2f}%)")
    
    # --- 如果有多重變異，列出範例 ---
    if multi_count > 0:
        print("\n[發現多重變異！以下是前 5 筆範例]")
        print("-" * 60)
        # 顯示 chrom, pos, ref, alt 四個欄位
        examples = df[multi_allelic_mask][['chrom', 'pos', 'ref', 'alt']].head(5)
        print(examples.to_string(index=False))
        print("-" * 60)
        print("說明: 這些變異在 VCF 中原本就是多等位基因 (Multi-allelic)。")
    else:
        print("\n[結果]")
        print("檢查通過：目前資料中沒有發現包含逗號的多重變異。")

if __name__ == "__main__":
    check_multiallelic()