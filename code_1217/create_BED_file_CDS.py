import sys

def parse_attributes(attr_str):
    """解析 GTF 第9欄的屬性字串，回傳字典"""
    attributes = {}
    # 屬性通常格式為: gene_id "ABC"; transcript_id "NM_123";
    parts = attr_str.strip().split(';')
    for part in parts:
        part = part.strip()
        if not part:
            continue
        # 分割鍵值 (例如: transcript_id "NM_123")
        try:
            #以space分隔,只切一刀
            #split如果沒東西切會回傳原字串,這行會ValueError
            key, value = part.split(' ',1)
            #移除包住value的 " "
            value = value.strip('"')
            attributes[key] = value
        #如果上面格式不符合,pass
        except ValueError:
            print(f'警告:缺少屬性欄位:{parts}')
        
    return attributes

def main():
    gtf_file = '/home/czlin/ClinVar_project/data/hg38.ncbiRefSeq.gtf'
    bed_file = '/home/czlin/ClinVar_project/data_new/targets_CDS.bed'

    #僅處理chr1~22,chrX/Y/M
    valid_chroms = set()
    for i in range(1, 23):
        valid_chroms.add(f"chr{i}")
    valid_chroms.add("chrX")
    valid_chroms.add("chrY")
    valid_chroms.add("chrM")
    
    print(f"開始處理: {gtf_file} -> {bed_file}")
    
    count = 0
    count_of_incomplete = 0
    with open(gtf_file, 'r') as fin, open(bed_file, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                continue
                
            parts = line.strip().split('\t')

            if len(parts) < 9:
                count_of_incomplete+=1
                continue
            
            #過濾掉奇怪染色體
            chrom = parts[0]
            if chrom not in valid_chroms:
                continue

            # 我們只關心 '5UTR', 'CDS', '3UTR' 區域
            # RefSeq GTF 在第3欄標示 feature type
            feature_type = parts[2]
            
            # 定義要抓取的 feature types
            target_features = {'5UTR', 'CDS', '3UTR'}
            
            if feature_type in target_features:
                
                # GTF (1-based) -> BED (0-based)
                # BED start = GTF start - 1
                start = int(parts[3]) - 1
                end = int(parts[4])
                
                # 解析屬性取得 transcript_id
                #attribute=以種類為key的dict
                attributes = parse_attributes(parts[8])

                try:
                    transcript_id = attributes['transcript_id']
                    gene_id = attributes['gene_id']
                except Exception as e:
                    print(f'警告:{parts}的attributes抓不到應該存在的屬性:{e}')
                    continue
                
                # 組合 BED 輸出的 name 欄位: "NM_001.3|GeneA|CDS"
                # 包含 feature_type 以便後續區分區域類型
                name = f"{transcript_id}|{gene_id}|{feature_type}"
                
                # 寫入 BED: chrom, start, end, name
                fout.write(f"{chrom.replace('chr','')}\t{start}\t{end}\t{name}\n")
                count += 1
                
    print(f"完成！共輸出 {count} 個 5UTR/CDS/3UTR 區間。")
    print(f'不完整的數量:{count_of_incomplete}')

if __name__ == "__main__":
    main()