# src/processing.py
import pysam
import csv
import sys
import time
from typing import List

# 匯入我們定義好的模組
from . import config
from . import parsers
from .models import MatchedVariant

def process_bed_regions():
    """
    主處理函式：
    1. 讀取 BED 檔
    2. 對每個區間查詢 VCF
    3. 解析並寫入 TSV
    """
    print(f"開始處理...")
    print(f"BED 來源: {config.BED_FILE}")
    print(f"VCF 來源: {config.VCF_FILE}")
    print(f"輸出目標: {config.OUTPUT_FILE}")
    
    # 1. 準備 VCF 讀取器
    try:
        vcf_reader = pysam.VariantFile(config.VCF_FILE)
    except Exception as e:
        print(f"錯誤: 無法讀取 VCF 檔案。請檢查路徑與索引 (.tbi)。\n詳細錯誤: {e}")
        return

    # 統計計數器
    stats = {
        "processed_regions": 0,
        "matched_variants": 0,
        "start_time": time.time()
    }

    # 2. 開啟輸出檔案 (使用 TSV 格式)
    # newline='' 是為了防止在 Windows/某些環境下出現多餘空行
    with open(config.OUTPUT_FILE, 'w', newline='', encoding='utf-8') as out_f:
        writer = csv.writer(out_f, delimiter='\t')
        
        # 寫入標題列 (Header)
        writer.writerow(config.CSV_HEADER)
        
        # 3. 讀取 BED 檔案
        with open(config.BED_FILE, 'r') as bed_f:
            for line in bed_f:
                # 簡單的進度顯示
                stats["processed_regions"] += 1
                if stats["processed_regions"] % 5000 == 0:
                    elapsed = time.time() - stats["start_time"]
                    print(f"已處理 {stats['processed_regions']} 個 BED 區間... "
                          f"(目前配對到 {stats['matched_variants']} 個變異, 耗時 {elapsed:.1f}s)")

                parts = line.strip().split('\t')
                if len(parts) < 4:
                    continue
                
                # 解析 BED 欄位
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                transcript_raw = parts[3]
                
                # 使用 parser 解析 transcript info
                tr_info = parsers.parse_transcript_info(transcript_raw)
                
                # 4. 在 VCF 中搜尋重疊變異
                try:
                    # fetch 回傳的是一個 iterator
                    for record in vcf_reader.fetch(chrom, start, end):
                        # 建立唯一的 SNV 識別 (ref/alt)
                        # 注意: 這裡簡化處理，若有多個 ALT 會被視為同一筆 record
                        # 若需要嚴格拆分多個 ALT，需額外迴圈處理 record.alts
                        
                        # 提取並解析 VCF 資訊
                        mc_field = record.info.get("MC")
                        clnsig_field = record.info.get("CLNSIG")
                        
                        mc_list = parsers.parse_mc_value(mc_field)
                        mc_cat = parsers.categorize_mc(mc_list)
                        clnsig_cat = parsers.categorize_clnsig(clnsig_field)
                        
                        # 建立資料物件
                        variant = MatchedVariant(
                            chrom=record.chrom,
                            pos=record.pos,
                            ref=record.ref,
                            alt=",".join(record.alts) if record.alts else ".",
                            transcript_id=tr_info['transcript_id'],
                            feature_type=tr_info['feature_type'],
                            gene_name=tr_info['gene_name'],
                            clnsig_category=clnsig_cat,
                            mc_category=mc_cat,
                            clnsig_raw=str(clnsig_field) if clnsig_field else None,
                            mc_raw=str(mc_field) if mc_field else None
                        )
                        
                        # 寫入檔案
                        writer.writerow(variant.to_csv_row())
                        stats["matched_variants"] += 1
                        
                except ValueError as ve:
                    # 通常發生在染色體名稱不匹配 (例如 BED寫 chr1 但 VCF寫 1)
                    # 這裡選擇忽略錯誤並繼續，避免程式中斷
                    continue

    vcf_reader.close()
    
    # 結束報告
    total_time = time.time() - stats["start_time"]
    print("\n" + "="*40)
    print("處理完成!")
    print(f"總耗時: {total_time:.2f} 秒")
    print(f"總掃描 BED 區間: {stats['processed_regions']}")
    print(f"總提取變異數量: {stats['matched_variants']}")
    print(f"結果已儲存至: {config.OUTPUT_FILE}")
    print("="*40)