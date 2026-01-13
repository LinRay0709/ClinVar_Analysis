import pysam

def main():
    # 輸入與輸出檔名
    input_vcf = "/home/czlin/ClinVar_project/data/clinvar.vcf.gz"  # 若是壓縮檔請改為 "clinvar.vcf.gz"
    output_vcf = "/home/czlin/ClinVar_project/data/clinvar_snv_pysam.vcf"
    
    print(f"開始篩選 SNV: {input_vcf} -> {output_vcf}")
    
    # 開啟 VCF 檔案
    # mode 'r'讀取純文字, 'rb'讀取壓縮檔。pysam 會自動偵測，通常用 'r' 即可
    vcf_in = pysam.VariantFile(input_vcf)
    
    # 建立輸出檔案，直接複製輸入檔的檔頭 (header)
    vcf_out = pysam.VariantFile(output_vcf, 'w', header=vcf_in.header)
    
    count_total = 0
    count_kept = 0
    
    for record in vcf_in:
        count_total += 1
        
        # --- SNV 判斷邏輯 ---
        # 1. Reference 鹼基長度必須為 1
        if len(record.ref) != 1:
            continue
            
        # 2. 檢查所有的 Alternative 鹼基
        # ClinVar 有時一行會有多個 Alt (例如 A -> C, G)
        # 我們只保留 "所有 Alt 長度都為 1" 的紀錄
        is_snv = True
        
        if not record.alts:
            continue

        for alt in record.alts:
            if len(alt) != 1:
                is_snv = False
                break
        
        if is_snv:
            vcf_out.write(record)
            count_kept += 1
            
        if count_total % 100000 == 0:
            print(f"已處理 {count_total} 筆資料...", end='\r')
            
    vcf_in.close()
    vcf_out.close()
    
    print(f"\n完成！")
    print(f"原始筆數: {count_total}")
    print(f"SNV筆數 : {count_kept}")
    print(f"篩選率 : {count_kept/count_total:.2%}")

if __name__ == "__main__":
    main()