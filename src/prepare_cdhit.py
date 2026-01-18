# src/prepare_cdhit.py
import sys
import os
import pandas as pd

import config

# ==========================================
# 參數設定
# ==========================================

# 篩選條件 (這些屬於邏輯參數，可以留在這裡，或是也移到 config)
TARGET_FEATURE = "CDS"
TARGET_MC = "synonymous_variant"
GROUP_PATHOGENIC = ["Pathogenic", "pathogenic"]
GROUP_BENIGN = ["Benign", "benign"]
# ==========================================

def main():
    print("="*60)
    print("CD-HIT 準備程式: TSV -> FASTA")
    print("="*60)
    
    # 0. 確保輸出目錄存在
    # [修改] 使用 config.INTERIM_DIR
    os.makedirs(config.INTERIM_DIR, exist_ok=True)

    # 1. 讀取資料
    input_file = config.SEQUENCES_FILE
    if not os.path.exists(input_file):
        print(f"錯誤: 找不到輸入檔案 {input_file}")
        return

    print(f"正在讀取: {input_file} ...")
    
    # [修改] 使用 config.CDHIT_PREPARE_READ_COLS
    df = pd.read_csv(
        input_file, 
        sep='\t', 
        usecols=config.CDHIT_PREPARE_READ_COLS, 
        dtype=str
    )
    
    total_raw = len(df)
    print(f"原始資料筆數: {total_raw}")

    # 2. 進行篩選 (Filter)
    print("\n[篩選條件]")
    print(f"  - Feature Type: {TARGET_FEATURE}")
    print(f"  - MC Category : {TARGET_MC}")
    print(f"  - Clinical Sig: Pathogenic/Benign")

    # 2.1 篩選 CDS 和 Synonymous
    # 轉小寫比較保險，或是用 str.contains
    mask_cds = df['feature_type'] == TARGET_FEATURE
    # 這裡使用 contains 來容許像 "synonymous_variant&splice_region" 這種複合情況
    mask_syn = df['mc_category'] == TARGET_MC    
    df_filtered = df[mask_cds & mask_syn].copy()
    print(f"  -> 符合 CDS & Synonymous 的筆數: {len(df_filtered)}")

    # 2.2 標記良惡性 (Labeling)
    def assign_label(clnsig):
        if pd.isna(clnsig): return None
        # 移除多餘空白並標準化
        sig = str(clnsig).strip()
        
        # 檢查是否屬於 Pathogenic 群
        if sig in GROUP_PATHOGENIC:
            return 1
        # 檢查是否屬於 Benign 群
        if sig in GROUP_BENIGN:
            return 0
        return None

    df_filtered['label'] = df_filtered['clnsig_category'].apply(assign_label)
    
    # 移除無法分類(CLNSIG不為良/惡性)的資料 (Label 為 None)
    before_dropna = len(df_filtered)
    df_final = df_filtered.dropna(subset=['label'])
    after_dropna = len(df_final)
    count_dropna = before_dropna - after_dropna

    # 轉換 label 為整數 (0 或 1)
    df_final['label'] = df_final['label'].astype(int)

    # 3. 去重複 (Deduplication)
    # 根據 snv_key 去重。因為基因組序列只看座標，同一座標的序列是一樣的
    # 把重複的存入redundant_snvs.tsv
    duplicate_mask = df_final.duplicated(subset=['snv_key'], keep='first')
    redundant_df = df_final[duplicate_mask]
    if len(redundant_df) > 0:
        redundant_path = os.path.join(config.INTERIM_DIR, "redundant_snvs.tsv")
        print(f"  -> 將冗餘資料備份至: {redundant_path}")
        redundant_df.to_csv(redundant_path, sep='\t', index=False)

    clean_df = df_final[~duplicate_mask]
    before_dedup = len(df_final)
    df_final = clean_df
    after_dedup = len(df_final)
    
    # 統計
    count_pos = (df_final['label'] == 1).sum()
    count_neg = (df_final['label'] == 0).sum()
    
    print("\n[最終統計]")
    print(f"  - dropna前/後筆數: {before_dropna}/{after_dropna}(差距:{count_dropna})")
    print(f"  - 篩選後筆數 (含重複): {before_dedup}")
    print(f"  - 去重後筆數 (Unique SNVs): {after_dedup} (移除 {before_dedup - after_dedup} 筆重複)")
    print(f"  - Pathogenic (Label=1): {count_pos}")
    print(f"  - Benign     (Label=0): {count_neg}")

    if after_dedup == 0:
        print("\n[警告] 篩選後沒有資料！請檢查 mc_category 或 clnsig_category 的關鍵字是否正確。")
        # 印出一些範例幫助除錯
        print("資料中的 mc_category 範例:", df['mc_category'].unique()[:5])
        print("資料中的 clnsig_category 範例:", df['clnsig_category'].unique()[:5])
        return

    # 4. 寫入 FASTA
    # [修改] 使用 config.CDHIT_FASTA_FILE
    print(f"\n正在寫入 FASTA: {config.CDHIT_FASTA_FILE} ...")
    with open(config.CDHIT_FASTA_FILE, 'w') as f_out:
        for idx, row in df_final.iterrows():
            header = f">{row['snv_key']}|{row['label']}"
            sequence = row['alt_seq']
            f_out.write(f"{header}\n")
            f_out.write(f"{sequence}\n")

    print(f"完成！請使用以下指令執行 CD-HIT:")
    print("-" * 60)
    # [修改] 使用 config.CDHIT_FASTA_FILE 和 config.CDHIT_OUTPUT_PREFIX
    print(f"cd-hit-est -i {config.CDHIT_FASTA_FILE} -o {config.CDHIT_OUTPUT_PREFIX} -c 0.9 -n 8 -M 16000 -T 8")
    print("-" * 60)

if __name__ == "__main__":
    main()