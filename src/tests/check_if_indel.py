# tests/check_variant_types.py
# 確認是否都是SNV(ref and alt == 1) -> 都是
import sys
import os
import pandas as pd

# --- 路徑設定 ---
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

import config

def main():
    print("="*60)
    print("變異類型分佈檢查 (SNP vs Indel)")
    print("="*60)
    
    if not os.path.exists(config.OUTPUT_FILE):
        print(f"找不到檔案: {config.OUTPUT_FILE}")
        return

    print(f"正在讀取: {config.OUTPUT_FILE} ...")
    # dtype=str 確保序列不會被當成奇怪的格式讀入
    df = pd.read_csv(config.OUTPUT_FILE, sep='\t', dtype=str)
    
    # 填充空值以免報錯 (雖然理論上不該有)
    df['ref'] = df['ref'].fillna('')
    df['alt'] = df['alt'].fillna('')
    
    # 計算長度
    df['ref_len'] = df['ref'].apply(len)
    df['alt_len'] = df['alt'].apply(len)
    
    # --- 分類邏輯 ---
    # 1. SNP: ref 和 alt 都只有 1 個鹼基
    is_snp = (df['ref_len'] == 1) & (df['alt_len'] == 1)
    
    # 2. Indel: 長度不一樣 (包含 Insertion 和 Deletion)
    is_indel = (df['ref_len'] != df['alt_len'])
    
    # 3. MNP (Multi-Nucleotide Polymorphism): 長度一樣但都大於 1 (例如 AT -> GC)
    is_mnp = (df['ref_len'] == df['alt_len']) & (df['ref_len'] > 1)
    
    snps = df[is_snp]
    indels = df[is_indel]
    mnps = df[is_mnp]
    
    # --- 輸出報告 ---
    total = len(df)
    print(f"\n[統計結果] 總筆數: {total}")
    print("-" * 40)
    print(f"SNP (單鹼基替換) : {len(snps):>8} ({len(snps)/total*100:.2f}%)")
    print(f"Indel (插入/缺失)  : {len(indels):>8} ({len(indels)/total*100:.2f}%)")
    print(f"MNP (多鹼基替換)   : {len(mnps):>8} ({len(mnps)/total*100:.2f}%)")
    print("-" * 40)
    
    # --- 展示範例 ---
    if len(indels) > 0:
        print("\n[發現 Indel！前 5 筆範例]")
        print(indels[['chrom', 'pos', 'ref', 'alt']].head(5).to_string(index=False))
        print("\n注意：這些變異將導致 Alt 序列總長度不等於 257 bp。")
        
    if len(mnps) > 0:
        print("\n[發現 MNP！前 5 筆範例]")
        print(mnps[['chrom', 'pos', 'ref', 'alt']].head(5).to_string(index=False))

    if len(indels) == 0 and len(mnps) == 0:
        print("\n[結論] 資料集中僅包含純 SNP，輸出長度將會完全固定。")

if __name__ == "__main__":
    main()