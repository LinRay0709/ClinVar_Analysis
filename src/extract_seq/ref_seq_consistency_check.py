# src/extract_seq/ref_seq_consistency_check.py
import sys
import os
import pandas as pd

# --- 路徑設定 (關鍵修改) ---
# 1. 目前位置: .../project/src/extract_seq/
current_dir = os.path.dirname(os.path.abspath(__file__))
# 2. 上一層: .../project/src/
src_dir = os.path.dirname(current_dir)
# 3. 再上一層 (專案根目錄): .../project/
project_root = os.path.dirname(src_dir)

# 將專案根目錄加入搜尋路徑，這樣才能 import src.config
if project_root not in sys.path:
    sys.path.append(project_root)

from src import config
# -------------------------

def load_fasta_sequence(chrom_name):
    """
    讀取單一染色體的 FASTA 檔案到記憶體字串。
    依賴 config.REF_DIR
    """
    # 組合完整路徑
    # 假設 TSV 是 "1", 檔案是 "chr1.fa"
    target_file = os.path.join(config.REF_DIR, f"chr{chrom_name}.fa")
    
    if not os.path.exists(target_file):
        print(f"  [警告] 找不到 FASTA 檔案: {target_file}")
        # 嘗試印出 config.REF_DIR 幫助除錯
        print(f"  (搜尋目錄為: {config.REF_DIR})")
        return None

    print(f"  正在載入序列: {os.path.basename(target_file)} ...", end="", flush=True)
    
    try:
        with open(target_file, 'r') as f:
            # 讀取 header (>chr1)
            header = f.readline()
            #如果第一行不是註解跳回開頭
            if not header.startswith('>'):
                f.seek(0)
            
            # 讀取序列 (去換行, 轉大寫)
            lines = []
            for line in f:
                cleaned_line = line.strip().upper()
                if cleaned_line:
                    lines.append(cleaned_line)
        sequence = "".join(lines)
            
        print(f" 完成 (長度: {len(sequence):,})")
        return sequence
        
    except Exception as e:
        print(f"  [錯誤] 讀取失敗: {e}")
        return None

def main():
    print("="*60)
    print("Reference Genome 一致性檢查 (Verify Ref Base)")
    print("="*60)
    print(f"Reference 目錄: {config.REF_DIR}")

    # 1. 讀取 TSV
    if not os.path.exists(config.OUTPUT_FILE):
        print(f"找不到輸入檔案: {config.OUTPUT_FILE}")
        return

    print("正在讀取變異資料表...")
    #強制chrom欄要用str讀,避免把1,2...讀成int
    df = pd.read_csv(config.OUTPUT_FILE, sep='\t', dtype={'chrom': str})
    
    # 2. 資料排序 (為了 IO 效能，必須排序)
    def sort_key(x):
        try:
            return int(x)
        except:
            return 999 
    #新增一行_sort_key,值用sort_key去算        
    df['_sort_key'] = df['chrom'].apply(sort_key)
    # 排序邏輯: 先排染色體(按照_sort_key的1~22,999...排列) -> 再排位置(同染色體內按位置)
    df = df.sort_values(by=['_sort_key', 'pos'])
    
    total_mismatches = 0
    total_checked = 0
    
    # 3. 依染色體分組處理
    #df.groupby:叫Pandas按照chrom這欄分組,每次給chrom和group(屬於此chrom的所有資料,是個小的df)
    # sort=False 是因為我們上面已經手動排好順序了
    for chrom, group in df.groupby('chrom', sort=False):
        print(f"\n處理染色體: {chrom} (共 {len(group)} 筆變異)")
        
        # 載入該chrom的序列(str)
        seq_str = load_fasta_sequence(chrom)
        
        if seq_str is None:
            continue
            
        mismatch_count = 0
        
        for idx, row in group.iterrows():
            pos_1based = int(row['pos'])
            ref_base = row['ref']
            
            # 轉換座標 (1-based -> 0-based)
            pos_0based = pos_1based - 1
            
            # 邊界檢查
            if pos_0based >= len(seq_str):
                print(f"  [異常] 座標超出範圍! Pos: {pos_1based}")
                mismatch_count += 1
                continue
                
            # 比對
            genome_base = seq_str[pos_0based]
            ref_to_check = ref_base[0] 
            
            if genome_base != ref_to_check:
                if mismatch_count < 5:
                    print(f"  [Mismatch] {chrom}:{pos_1based} | VCF: {ref_base} vs FASTA: {genome_base}")
                mismatch_count += 1
                
        if mismatch_count == 0:
            print("  -> OK (全數匹配)")
        else:
            print(f"  -> 發現 {mismatch_count} 筆不匹配！")
            
        total_checked += len(group)
        total_mismatches += mismatch_count
        
        # 釋放記憶體
        del seq_str

    print("\n" + "="*60)
    print("檢查總結")
    print(f"總檢查筆數: {total_checked}")
    print(f"總不匹配數: {total_mismatches}")
    
    if total_mismatches == 0:
        print("結果: [PASS] 完美一致！")
    else:
        print("結果: [FAIL] 請檢查 Reference 版本或座標系統。")

if __name__ == "__main__":
    main()