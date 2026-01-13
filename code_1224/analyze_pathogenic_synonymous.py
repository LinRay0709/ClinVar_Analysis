#!/usr/bin/env python3
"""
分析只配對到單一 CDS 區段的 SNV 變異。

此腳本的目的是：
1. 找出只配對到一個 CDS 區段的 SNV (排除同時配對到多個 CDS 的變異)
2. 分析這些 SNV 的 CLNSIG 良惡性分布
3. 分析各 CLNSIG 分類下的 MC (Molecular Consequence) 分布
"""

import pysam
from collections import defaultdict, Counter


def parse_mc_value(mc_field) -> list:
    """
    解析 MC (Molecular Consequence) 欄位值。
    
    MC 格式範例: "SO:0001627|intron_variant" 或多個值的 tuple
    
    Args:
        mc_field: MC 欄位值 (可能是 str 或 tuple)
    
    Returns:
        list: 解析出的 consequence 類型列表
    """
    consequences = []
    
    if mc_field is None:
        return consequences
    
    # 轉換為列表處理
    if isinstance(mc_field, str):
        mc_values = [mc_field]
    else:
        mc_values = list(mc_field)
    
    for mc in mc_values:
        # MC 格式: "SO:0001627|intron_variant"
        # 取 | 後面的部分作為 consequence type
        if '|' in mc:
            consequence = mc.split('|', 1)[1]
            consequences.append(consequence)
        else:
            consequences.append(mc)
    
    return consequences


def categorize_mc(mc_consequences: list) -> str:
    """
    將 MC 值列表分類為 conflicting/missense_variant/synonymous_variant/other。
    
    分類邏輯:
    - conflicting: 同時含有 missense_variant 和 synonymous_variant
    - missense_variant: 含有 missense_variant (無 synonymous)
    - synonymous_variant: 含有 synonymous_variant (無 missense)
    - other: 其他所有值
    
    Args:
        mc_consequences: MC consequence 類型列表
    
    Returns:
        str: 分類結果
    """
    if not mc_consequences:
        return "other"
    
    has_missense = "missense_variant" in mc_consequences
    has_synonymous = "synonymous_variant" in mc_consequences
    
    # 檢查是否同時有 missense 和 synonymous
    if has_missense and has_synonymous:
        return "conflicting"
    
    # 只有 missense
    if has_missense:
        return "missense_variant"
    
    # 只有 synonymous
    if has_synonymous:
        return "synonymous_variant"
    
    return "other"


def categorize_clnsig(clnsig_field) -> str:
    """
    將 CLNSIG 值分類為 pathogenic/likely pathogenic/benign/likely benign/other。
    
    處理:
    - tuple 值 (多個值)
    - 以 '|' 分隔的多個值
    - 以 '/' 分隔的復合狀況 (如 "Pathogenic/Likely_pathogenic")
    
    優先順序:
    1. pathogenic (如果含有 "pathogenic")
    2. benign (如果含有 "benign")
    3. likely pathogenic (如果含有 "likely_pathogenic")
    4. likely benign (如果含有 "likely_benign")
    5. other
    
    Args:
        clnsig_field: CLNSIG 欄位值 (可能是 str 或 tuple)
    
    Returns:
        str: 分類結果 ("pathogenic", "likely pathogenic", "benign", "likely benign", 或 "other")
    """
    if clnsig_field is None:
        return "other"
    
    # 如果是 tuple，先用 '|' 連接成字串
    if isinstance(clnsig_field, tuple):
        clnsig_str = "|".join(clnsig_field)
    else:
        clnsig_str = clnsig_field
    
    # 先用 '|' 分割多個值
    pipe_separated = clnsig_str.split('|')
    
    # 再用 '/' 分割復合狀況，並展開成列表
    all_values = []
    for value in pipe_separated:
        slash_separated = value.split('/')
        all_values.extend(slash_separated)
    
    # 轉為小寫比對
    values_lower = [v.strip().lower() for v in all_values]
    
    # 檢查各類別
    has_pathogenic = "pathogenic" in values_lower
    has_benign = "benign" in values_lower
    
    # 優先檢查 pathogenic
    if has_pathogenic:
        return "pathogenic"
    
    # 再檢查 benign
    if has_benign:
        return "benign"
    
    # 檢查 likely_pathogenic
    if "likely_pathogenic" in values_lower:
        return "likely pathogenic"
    
    # 檢查 likely_benign
    if "likely_benign" in values_lower:
        return "likely benign"
    
    return "other"


def analyze_clnsig_distribution(clnsig_counter: Counter, title: str):
    """
    分析並顯示 CLNSIG 的分布情況 (5 類)。
    
    Args:
        clnsig_counter: Counter of CLNSIG categories
        title: 標題
    """
    print("\n" + "=" * 70)
    print(f"CLNSIG 分布分析 - {title}")
    print("=" * 70)
    print("\n分類定義:")
    print("  pathogenic        : Pathogenic (優先)")
    print("  benign            : Benign (優先)")
    print("  likely pathogenic : Likely_pathogenic")
    print("  likely benign     : Likely_benign")
    print("  other             : 所有其他值")
    
    total = sum(clnsig_counter.values())
    
    print(f"\n共 {total} 個 SNV")
    print("-" * 50)
    
    if total == 0:
        print("  (無資料)")
        return
    
    print(f"  {'Category':<20} {'Count':>10} {'Percentage':>12}")
    print("  " + "-" * 45)
    
    # 按固定順序顯示: pathogenic, likely pathogenic, benign, likely benign, other
    for category in ['pathogenic', 'likely pathogenic', 'benign', 'likely benign', 'other']:
        count = clnsig_counter.get(category, 0)
        pct = (count / total * 100)
        print(f"  {category:<20} {count:>10} {pct:>11.2f}%")
    
    print("  " + "-" * 45)
    print(f"  {'Total':<20} {total:>10}")


def analyze_mc_by_clnsig(mc_by_clnsig: dict, title: str):
    """
    分析並顯示各 CLNSIG 分類下的 MC 分布情況。
    
    Args:
        mc_by_clnsig: {clnsig_category: Counter of MC categories}
        title: 標題
    """
    print("\n" + "=" * 70)
    print(f"各 CLNSIG 分類下的 MC 分布 - {title}")
    print("=" * 70)
    print("\nMC 分類定義:")
    print("  conflicting        : 同時含有 missense_variant 和 synonymous_variant")
    print("  missense_variant   : 含有 missense_variant (無 synonymous)")
    print("  synonymous_variant : 含有 synonymous_variant (無 missense)")
    print("  other              : 其他所有值")
    
    for clnsig_category in ['pathogenic', 'likely pathogenic', 'benign', 'likely benign', 'other']:
        mc_counter = mc_by_clnsig.get(clnsig_category, Counter())
        total = sum(mc_counter.values())
        
        print(f"\n【{clnsig_category}】SNV 的 MC 分布 (共 {total} 個 SNV)")
        print("-" * 55)
        
        if total == 0:
            print("  (無資料)")
            continue
        
        print(f"  {'MC Category':<25} {'Count':>10} {'Percentage':>12}")
        print("  " + "-" * 50)
        
        # 按固定順序顯示
        for mc_category in ['conflicting', 'missense_variant', 'synonymous_variant', 'other']:
            count = mc_counter.get(mc_category, 0)
            pct = (count / total * 100)
            print(f"  {mc_category:<25} {count:>10} {pct:>11.2f}%")
        
        print("  " + "-" * 50)
        print(f"  {'Total':<25} {total:>10}")


def main():
    bed_file = "/home/czlin/ClinVar_project/data_new/targets_CDS.bed"
    vcf_file = "/home/czlin/ClinVar_project/data_new/clinvar_snv_new.vcf.gz"
    
    print(f"分析只配對到單一 CDS 區段的 SNV")
    print(f"BED 檔案: {bed_file}")
    print(f"VCF 檔案: {vcf_file}")
    print("=" * 70)
    
    # 開啟 VCF 檔案 (必須有 .tbi 索引檔在同目錄下)
    try:
        vcf_reader = pysam.VariantFile(vcf_file)
    except Exception as e:
        print(f"錯誤: 無法讀取 VCF，請確認是否已壓縮並建立索引 (.tbi)。\n{e}")
        return

    # 第一輪：收集所有 SNV 配對到的 CDS 區段
    # {snv_key: [transcript_info, ...]} - 記錄每個 SNV 配對到的 CDS 區段
    snv_cds_matches = defaultdict(list)
    
    # 追蹤每個 SNV 的 MC 原始值和 CLNSIG 原始值
    snv_mc_raw = {}
    snv_clnsig_raw = {}
    
    # 記錄 SNV 配對到的 BED 區段資訊 (用於驗證)
    snv_bed_info = defaultdict(list)
    
    bed_region_count = 0
    cds_region_count = 0
    
    print("\n正在掃描所有配對...")
    
    with open(bed_file, 'r') as bed_in:
        for line in bed_in:
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
                
            bed_region_count += 1
            
            # 解析 BED: chrom, start, end, name
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            transcript_info = parts[3]
            
            # 從 transcript_info 解析 feature_type
            info_parts = transcript_info.split('|')
            if len(info_parts) >= 3:
                feature_type = info_parts[2]  # CDS, 5UTR, 或 3UTR
            else:
                continue
            
            # 只處理 CDS 區段
            if feature_type != 'CDS':
                continue
            
            cds_region_count += 1
            
            # 使用 fetch 找重疊的 SNV
            try:
                for record in vcf_reader.fetch(chrom, start, end):
                    pos = record.pos
                    ref = record.ref
                    alt = ",".join(record.alts) if record.alts else "."
                    
                    # 建立 SNV 唯一識別鍵
                    snv_key = f"{chrom}:{pos}:{ref}:{alt}"
                    
                    # 記錄此 SNV 配對到的 CDS 區段
                    snv_cds_matches[snv_key].append(transcript_info)
                    
                    # 記錄 BED 區段資訊
                    bed_region_str = f"{chrom}:{start}-{end}"
                    snv_bed_info[snv_key].append((bed_region_str, transcript_info))
                    
                    # 只在第一次遇到這個 SNV 時記錄 MC 和 CLNSIG
                    if snv_key not in snv_mc_raw:
                        snv_mc_raw[snv_key] = record.info.get("MC", None)
                        snv_clnsig_raw[snv_key] = record.info.get("CLNSIG", None)
                    
            except ValueError:
                continue
    
    vcf_reader.close()
    
    print(f"已處理 {bed_region_count} 個 BED 區間")
    print(f"其中 CDS 區段: {cds_region_count} 個")
    print(f"總共找到 {len(snv_cds_matches)} 個配對到 CDS 的不重複 SNV")
    
    # 第二輪：篩選只配對到單一 CDS 區段的 SNV
    single_cds_snvs = []
    multi_cds_snvs = []
    
    for snv_key, cds_list in snv_cds_matches.items():
        if len(cds_list) == 1:
            single_cds_snvs.append(snv_key)
        else:
            multi_cds_snvs.append(snv_key)
    
    print(f"\n篩選結果:")
    print(f"  只配對到單一 CDS 的 SNV: {len(single_cds_snvs)} 個")
    print(f"  配對到多個 CDS 的 SNV: {len(multi_cds_snvs)} 個")
    
    # 輸出前三個配對供驗證
    print("\n" + "=" * 80)
    print("配對正確性驗證 (前 3 個只配對到單一 CDS 的 SNV)")
    print("=" * 80)
    
    for i, snv_key in enumerate(single_cds_snvs[:3], 1):
        bed_infos = snv_bed_info[snv_key]
        cds_list = snv_cds_matches[snv_key]
        
        print(f"\n配對 {i}:")
        print(f"  SNV 座標       : {snv_key}")
        print(f"  配對的 CDS 數量: {len(cds_list)}")
        for bed_region, transcript in bed_infos[:1]:
            print(f"  BED 區段座標   : {bed_region}")
            print(f"  Transcript     : {transcript}")
    
    # 統計 CLNSIG 分布
    clnsig_counter = Counter()
    
    # 統計 MC 分布 (按 CLNSIG 分類)
    mc_by_clnsig = {
        'pathogenic': Counter(),
        'likely pathogenic': Counter(),
        'benign': Counter(),
        'likely benign': Counter(),
        'other': Counter()
    }
    
    # 分析每個只配對到單一 CDS 的 SNV
    for snv_key in single_cds_snvs:
        # 分類 CLNSIG
        clnsig_category = categorize_clnsig(snv_clnsig_raw.get(snv_key))
        clnsig_counter[clnsig_category] += 1
        
        # 分類 MC
        mc_consequences = parse_mc_value(snv_mc_raw.get(snv_key))
        mc_category = categorize_mc(mc_consequences)
        mc_by_clnsig[clnsig_category][mc_category] += 1
    
    # 輸出統計結果
    print("\n" + "=" * 70)
    print("SNV 統計結果")
    print("=" * 70)
    print(f"\n{'類型':<25} {'Count':>10}")
    print("-" * 40)
    print(f"{'只配對到單一 CDS 的 SNV':<25} {len(single_cds_snvs):>10}")
    print(f"{'配對到多個 CDS 的 SNV':<25} {len(multi_cds_snvs):>10}")
    print("-" * 40)
    print(f"{'總計':<25} {len(snv_cds_matches):>10}")
    
    # 輸出 CLNSIG 分布分析
    analyze_clnsig_distribution(clnsig_counter, "只配對到單一 CDS 的 SNV")
    
    # 輸出 MC 分布分析 (按 CLNSIG 分類)
    analyze_mc_by_clnsig(mc_by_clnsig, "只配對到單一 CDS 的 SNV")


if __name__ == "__main__":
    main()
