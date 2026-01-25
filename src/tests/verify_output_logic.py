# tests/verify_output_logic.py
#檢查all_sequence.tsv的ref_seq和alt_seq是否都是257bases -> 中間是否等於ref/alt -> ref_seq和alt_seq去掉中間後是否相同

import sys
import os
import pandas as pd

# 路徑設定
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

import config

def main():
    print("="*60)
    print("輸出結果邏輯驗證 (Automated Logic Check)")
    print("="*60)
    
    output_file = config.SEQUENCES_FILE
    if not os.path.exists(output_file):
        print(f"找不到檔案: {output_file}")
        return

    print(f"正在讀取: {output_file}")
    # 讀取時將所有欄位視為字串，確保不會遺漏
    df = pd.read_csv(output_file, sep='\t', dtype=str)
    
    total = len(df)
    errors = 0
    
    print(f"開始檢查 {total} 筆資料...")
    
    for idx, row in df.iterrows():
        try:
            # 1. 檢查長度
            if len(row['ref_seq']) != 1024 or len(row['alt_seq']) != 1024:
                print(f"[錯誤] 行 {idx+2}: 序列長度不正確 (Ref: {len(row['ref_seq'])}, Alt: {len(row['alt_seq'])})")
                errors += 1
                if errors > 5: break
                continue

            # 2. 檢查中心點 (0-based index 512 is the 512th base)
            center_idx = 512
            ref_center = row['ref_seq'][center_idx]
            alt_center = row['alt_seq'][center_idx]
            
            if ref_center != row['ref']:
                print(f"[錯誤] 行 {idx+2}: Ref 序列中心 '{ref_center}' 與 Ref 欄位 '{row['ref']}' 不符")
                errors += 1
            
            if alt_center != row['alt']:
                # 這裡可能會抓到那 33 個 N，確認是否為預期行為
                print(f"[錯誤] 行 {idx+2}: Alt 序列中心 '{alt_center}' 與 Alt 欄位 '{row['alt']}' 不符")
                errors += 1

            # 3. 檢查周邊序列一致性 (Context Identity)
            # Ref 序列挖掉中間 vs Alt 序列挖掉中間，應該要一樣
            ref_context = row['ref_seq'][:center_idx] + row['ref_seq'][center_idx+1:]
            alt_context = row['alt_seq'][:center_idx] + row['alt_seq'][center_idx+1:]
            
            if ref_context != alt_context:
                print(f"[錯誤] 行 {idx+2}: Ref 和 Alt 的周邊序列不一致！(這不該發生在 SNP)")
                errors += 1
                
        except Exception as e:
            print(f"[異常] 行 {idx+2} 發生未預期錯誤: {e}")
            errors += 1
            if errors > 5: break

    print("-" * 40)
    if errors == 0:
        print("✅ 驗證通過！所有序列的長度、中心點與周邊一致性皆正確。")
    else:
        print(f"❌ 驗證失敗，發現 {errors} 個錯誤。")

    # 順便看一下那 33 個 N 是誰
    n_alts = df[df['alt'] == 'N']
    if len(n_alts) > 0:
        print(f"\n[資訊] 發現 {len(n_alts)} 筆 Alt 為 'N' 的資料。範例:")
        print(n_alts[['chrom', 'pos', 'ref', 'alt']].head(3).to_string(index=False))

     # 新增功能：檢查 ref_seq 或 alt_seq 含有 @ 符號的資料
    at_in_ref = df[df['ref_seq'].str.contains('@', na=False)]
    at_in_alt = df[df['alt_seq'].str.contains('@', na=False)]
    at_rows = pd.concat([at_in_ref, at_in_alt]).drop_duplicates()
    
    if len(at_rows) > 0:
        print(f"\n[警告] 發現 {len(at_rows)} 筆 ref_seq 或 alt_seq 含有 '@' 符號的資料:")
        print(at_rows[['snv_key', 'chrom', 'pos', 'ref', 'alt']].to_string(index=False))
    else:
        print("\n[資訊] 沒有發現 ref_seq 或 alt_seq 含有 '@' 符號的資料。")

    # 新增功能：檢查 ref_seq 或 alt_seq 含有小寫字母的資料
    # 使用正則表達式 [a-z] 檢查是否有小寫字母
    lower_in_ref = df[df['ref_seq'].str.contains('[a-z]', na=False, regex=True)]
    lower_in_alt = df[df['alt_seq'].str.contains('[a-z]', na=False, regex=True)]
    lower_rows = pd.concat([lower_in_ref, lower_in_alt]).drop_duplicates()
    
    if len(lower_rows) > 0:
        print(f"\n[警告] 發現 {len(lower_rows)} 筆 ref_seq 或 alt_seq 含有小寫字母的資料:")
        print(lower_rows[['snv_key', 'chrom', 'pos', 'ref', 'alt']].head(10).to_string(index=False))
        if len(lower_rows) > 10:
            print(f"  ... 還有 {len(lower_rows) - 10} 筆未顯示")
    else:
        print("\n[資訊] 沒有發現 ref_seq 或 alt_seq 含有小寫字母的資料。")

if __name__ == "__main__":
    main()