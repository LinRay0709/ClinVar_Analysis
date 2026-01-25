# src/tests/check_pos_radius.py
# 檢查 matched_variants.tsv 中的 pos 是否在有效範圍內
# 有效範圍: (radius + 1) <= pos <= (染色體長度 - radius)
# 避免提取序列時超出邊界

import sys
import os

# --- 路徑設定 ---
current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.dirname(current_dir)
project_root = os.path.dirname(src_dir)

if project_root not in sys.path:
    sys.path.append(project_root)

from src import config
from src.extract_seq.ref_seq_consistency_check import load_fasta_sequence

# ==========================================
# 參數設定
# ==========================================
RADIUS = 512  # 前後各抓 512 bp

def main():
    print("=" * 60)
    print("SNV 位置邊界檢查 (Position Boundary Check)")
    print("=" * 60)
    print(f"  - Radius: {RADIUS}")
    print(f"  - 有效範圍: pos >= {RADIUS + 1} 且 pos <= (染色體長度 - {RADIUS} + 1)")
    print("=" * 60)

    # 1. 讀取配對資料
    if not os.path.exists(config.OUTPUT_FILE):
        print(f"錯誤: 找不到 {config.OUTPUT_FILE}")
        return

    # 使用 csv 模組避免依賴 pandas
    import csv
    
    # 讀取資料並按染色體分組
    chrom_data = {}  # chrom -> [(pos, line_num), ...]
    
    with open(config.OUTPUT_FILE, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        
        for line_num, row in enumerate(reader, start=2):  # start=2 因為 header 是第 1 行
            chrom = row['chrom']
            pos = int(row['pos'])
            
            if chrom not in chrom_data:
                chrom_data[chrom] = []
            chrom_data[chrom].append((pos, line_num))
    
    total_count = sum(len(v) for v in chrom_data.values())
    print(f"\n讀取完成，共 {total_count} 筆資料，涵蓋 {len(chrom_data)} 個染色體。\n")

    # 2. 依染色體檢查
    total_out_of_range = 0
    out_of_range_details = []  # [(chrom, pos, line_num, reason), ...]
    
    # 按染色體名自然排序 (1, 2, ..., 22, X, Y, ...)
    def sort_key(x):
        try:
            return (0, int(x))
        except:
            return (1, x)
    
    sorted_chroms = sorted(chrom_data.keys(), key=sort_key)
    
    for chrom in sorted_chroms:
        positions = chrom_data[chrom]
        print(f"處理染色體: {chrom} (共 {len(positions)} 筆變異)")
        
        # 載入該染色體的序列以獲取長度
        seq_str = load_fasta_sequence(chrom)
        
        if seq_str is None:
            print(f"  [跳過] 無法載入序列，跳過此染色體。")
            continue
        
        chrom_len = len(seq_str)
        
        # 計算有效範圍 (1-based)
        # pos 必須滿足: pos - radius >= 1 且 pos + radius < chrom_len
        # 即: pos >= radius + 1 且 pos < chrom_len - radius + 1
        min_valid_pos = RADIUS + 1
        max_valid_pos = chrom_len - RADIUS + 1
        
        chrom_out_of_range = 0
        
        for pos, line_num in positions:
            reason = None
            
            if pos < min_valid_pos:
                reason = f"pos ({pos}) < 最小有效位置 ({min_valid_pos})"
            elif pos > max_valid_pos:
                reason = f"pos ({pos}) > 最大有效位置 ({max_valid_pos})"
            
            if reason:
                chrom_out_of_range += 1
                out_of_range_details.append((chrom, pos, line_num, reason))
        
        if chrom_out_of_range == 0:
            print(f"  -> OK (全部在有效範圍內)")
        else:
            print(f"  -> 發現 {chrom_out_of_range} 筆超出有效範圍！")
        
        total_out_of_range += chrom_out_of_range
        
        # 釋放記憶體
        del seq_str

    # 3. 輸出總結
    print("\n" + "=" * 60)
    print("檢查總結")
    print("=" * 60)
    print(f"總檢查筆數: {total_count}")
    print(f"超出有效範圍: {total_out_of_range}")
    
    if total_out_of_range == 0:
        print("\n結果: [PASS] 所有 SNV 位置皆在有效範圍內！")
    else:
        print(f"\n結果: [WARNING] 有 {total_out_of_range} 筆 SNV 位置超出有效範圍。")
        print("\n超出範圍的詳細資料 (前 20 筆):")
        print("-" * 60)
        for i, (chrom, pos, line_num, reason) in enumerate(out_of_range_details[:20]):
            print(f"  Line {line_num}: chr{chrom}:{pos} - {reason}")
        
        if len(out_of_range_details) > 20:
            print(f"  ... 還有 {len(out_of_range_details) - 20} 筆未顯示")

if __name__ == "__main__":
    main()
