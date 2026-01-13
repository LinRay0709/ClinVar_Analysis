import sys
from collections import defaultdict

def analyze_bed_consistency(bed_file):
    print(f"正在讀取並分析: {bed_file} ...")
    
    # 資料結構: Key=GeneID, Value=Set(TranscriptIDs)
    # 使用 Set 可以自動去除重複 (因為一個 transcript 會有很多行 exon/CDS)
    gene_to_transcripts = defaultdict(set)
    unique_transcripts = set()
    
    line_count = 0
    
    try:
        with open(bed_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split('\t')
                
                # 確保 BED 至少有 4 欄
                if len(parts) < 4:
                    continue
                
                # 第 4 欄格式預期為: "NM_001.3|GeneID:123|CDS"
                name_col = parts[3]
                
                try:
                    # 根據您上一段程式碼的邏輯，用 '|' 分隔
                    # split 限制切 2 刀，確保取出前兩項即可，後面不重要
                    extracted_parts = name_col.split('|')
                    
                    if len(extracted_parts) >= 2:
                        transcript_id = extracted_parts[0] # NM_...
                        gene_id = extracted_parts[1]       # GeneID...
                        
                        # 紀錄
                        gene_to_transcripts[gene_id].add(transcript_id)
                        unique_transcripts.add(transcript_id)
                        line_count += 1
                    else:
                        # 格式不符的情況
                        pass
                        
                except Exception as e:
                    print(f"警告: 解析第 {line_count} 行發生錯誤: {name_col}, 錯誤: {e}")
                    continue

        # --- 輸出統計結果 ---
        total_genes = len(gene_to_transcripts)
        total_transcripts = len(unique_transcripts)
        
        print("-" * 40)
        print(f"📊 分析報告 (Analysis Report)")
        print("-" * 40)
        print(f"處理行數 (BED lines): {line_count}")
        print(f"總 Gene 數量:       {total_genes}")
        print(f"總 Transcript 數量: {total_transcripts}")
        print("-" * 40)

        # --- 關鍵檢查: 是否有一基因對應多 Transcript ---
        # 找出那些 set 長度 > 1 的 gene
        multi_transcript_genes = {
            gid: tids 
            for gid, tids in gene_to_transcripts.items() 
            if len(tids) > 1
        }
        
        if len(multi_transcript_genes) == 0:
            print("✅ 完美！驗證通過。")
            print("所有 Gene 都只對應唯一的 Transcript (1:1 Mapping)。")
            print("這代表您的 MANE Select 篩選與去重邏輯運作正常。")
        else:
            print(f"⚠️ 注意！發現 {len(multi_transcript_genes)} 個 Gene 擁有多個 Transcripts。")
            print("這可能是因為版本號不同 (如 .3 vs .4) 或篩選未完全。")
            print("\n前 10 個異常範例:")
            for i, (gid, tids) in enumerate(multi_transcript_genes.items()):
                #if i >= 10: break
                print(f"  - Gene: {gid}")
                print(f"    Transcripts: {', '.join(tids)}")

    except FileNotFoundError:
        print(f"❌ 錯誤: 找不到檔案 {bed_file}")
    except Exception as e:
        print(f"❌ 發生未預期的錯誤: {e}")

if __name__ == "__main__":
    # 請確認這是您剛剛生成的 BED 檔案路徑
    bed_file_path = '/home/czlin/ClinVar_project/data_0107/MANE_transcript.bed'
    
    analyze_bed_consistency(bed_file_path)