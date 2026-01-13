import pysam

def main():
    bed_file = "/home/czlin/ClinVar_project/data/targets.bed"
    vcf_file = "/home/czlin/ClinVar_project/data/clinvar_snv_pysam.vcf.gz" # 注意：必須讀取 .gz
    output_file = "/home/czlin/ClinVar_project/data/clinvar_refseq_pysam_result.txt"
    
    print(f"開始比對: {bed_file} vs {vcf_file}")
    
    # 開啟 VCF 檔案 (必須有 .tbi 索引檔在同目錄下)
    try:
        vcf_reader = pysam.VariantFile(vcf_file)
    except Exception as e:
        print(f"錯誤: 無法讀取 VCF，請確認是否已壓縮並建立索引 (.tbi)。\n{e}")
        return

    match_count = 0
    unmatch_count = 0
    
    with open(bed_file, 'r') as bed_in, open(output_file, 'w') as fout:
        # 寫入檔頭 (Header)
        fout.write("Chrom\tPos\tRef\tAlt\tClinVar_ID\tTranscript_Info\n")
        
        for line in bed_in:
            parts = line.strip().split('\t')
            # 解析 BED: chrom, start, end, name
            chrom = parts[0] #格式: 1~22,X,Y,M
            start = int(parts[1])
            end = int(parts[2])
            transcript_info = parts[3] #f"{transcript_id}|{gene_id})"
            
            # --- 核心搜尋邏輯 (fetch) ---
            try:
                # fetch(染色體, 起始, 結束)
                # 這比自己寫迴圈快上萬倍
                # chrom格式(chr1 vs 1)要跟要讀的vcf.gz檔內chrom格式ㄧ致
                for record in vcf_reader.fetch(chrom, start, end):
                    # 提取變異資訊
                    pos = record.pos
                    ref = record.ref
                    # 處理多重 Alt (通常 SNV 過濾後只有一個，但以防萬一)
                    if record.alts:
                        alt = ",".join(record.alts)
                    else:
                        print(f'clinvar中的{record}沒有alt')
                        alt = "."
                    if record.id:
                        var_id = record.id
                    else:
                        print(f'clinvar中的{record}沒有var_id')
                        var_id = "."
                    
                    # 寫入結果
                    fout.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{var_id}\t{transcript_info}\n")
                    match_count += 1
                    
            except ValueError:
                # 有時候 BED 的染色體 (如 chrM) 在 VCF 裡沒有紀錄，fetch 會報錯，跳過即可
                print(f'BED中transcript沒有找到對應clinvar record:{line}\n')
                unmatch_count += 1
                continue
                
    print(f"完成！共找到 {match_count} 個重疊的變異-轉錄本配對。")
    print(f"和 {unmatch_count} 個沒有record的transcript。")

if __name__ == "__main__":
    main()