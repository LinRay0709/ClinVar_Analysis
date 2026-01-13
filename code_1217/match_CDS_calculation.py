import pysam
from collections import defaultdict, Counter


# CLNSIG 分類定義 (case-insensitive)
PATHOGENIC_VALUES = {"pathogenic", "likely_pathogenic"}
BENIGN_VALUES = {"benign", "likely_benign"}


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
    將 CLNSIG 值分類為 Pathogenic/Benign/Other。
    
    處理:
    - tuple 值 (多個值)
    - 以 '|' 分隔的多個值
    - 以 '/' 分隔的復合狀況 (如 "Pathogenic/Likely_pathogenic")
    
    Args:
        clnsig_field: CLNSIG 欄位值 (可能是 str 或 tuple)
    
    Returns:
        str: 分類結果 ("Pathogenic", "Benign", 或 "Other")
    """
    if clnsig_field is None:
        return "Other"
    
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
    
    # 優先檢查 pathogenic
    for value in values_lower:
        if value in PATHOGENIC_VALUES:
            return "Pathogenic"
    
    # 再檢查 benign
    for value in values_lower:
        if value in BENIGN_VALUES:
            return "Benign"
    
    return "Other"


def analyze_clnsig_distribution(clnsig_by_feature: dict):
    """
    分析並顯示各 feature type 中 CLNSIG 的分布情況 (Benign/Pathogenic/Other)。
    
    Args:
        clnsig_by_feature: {feature_type: Counter of CLNSIG categories}
    """
    print("\n" + "=" * 70)
    print("CLNSIG 分布分析 (Benign / Pathogenic / Other)")
    print("=" * 70)
    print("\n分類定義:")
    print("  Pathogenic : Pathogenic, Likely_pathogenic")
    print("  Benign     : Benign, Likely_benign")
    print("  Other      : 所有其他值")
    
    for feature_type in ['CDS', '5UTR', '3UTR']:
        clnsig_counter = clnsig_by_feature.get(feature_type, Counter())
        total = sum(clnsig_counter.values())
        
        print(f"\n【{feature_type}】區域 (共 {total} 個 SNV)")
        print("-" * 50)
        
        if total == 0:
            print("  (無資料)")
            continue
        
        print(f"  {'Category':<20} {'Count':>10} {'Percentage':>12}")
        print("  " + "-" * 45)
        
        # 按固定順序顯示: Pathogenic, Benign, Other
        for category in ['Pathogenic', 'Benign', 'Other']:
            count = clnsig_counter.get(category, 0)
            pct = (count / total * 100)
            print(f"  {category:<20} {count:>10} {pct:>11.2f}%")
        
        print("  " + "-" * 45)
        print(f"  {'Total':<20} {total:>10}")


def analyze_mc_by_clnsig(mc_by_clnsig: dict):
    """
    分析並顯示 CDS 區域中，各 CLNSIG 分類下的 MC 分布情況。
    
    Args:
        mc_by_clnsig: {clnsig_category: Counter of MC categories}
    """
    print("\n" + "=" * 70)
    print("CDS 區域 - 各 CLNSIG 分類下的 MC 分布")
    print("=" * 70)
    print("\nMC 分類定義:")
    print("  conflicting        : 同時含有 missense_variant 和 synonymous_variant")
    print("  missense_variant   : 含有 missense_variant (無 synonymous)")
    print("  synonymous_variant : 含有 synonymous_variant (無 missense)")
    print("  other              : 其他所有值")
    
    for clnsig_category in ['Pathogenic', 'Benign', 'Other']:
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
    
    print(f"開始比對: {bed_file} vs {vcf_file}")
    
    # 開啟 VCF 檔案 (必須有 .tbi 索引檔在同目錄下)
    try:
        vcf_reader = pysam.VariantFile(vcf_file)
    except Exception as e:
        print(f"錯誤: 無法讀取 VCF，請確認是否已壓縮並建立索引 (.tbi)。\n{e}")
        return

    # 追蹤每個 SNV 配對到的 feature types (list 類型，記錄每次配對)
    # {snv_key: list of feature_types} - 同一 feature type 可能出現多次
    matched_snv_features = defaultdict(list)
    
    # 追蹤每個 SNV 的 MC 分類
    # {snv_key: "missense_variant" | "synonymous_variant" | "other"}
    snv_mc_category = {}
    
    # 追蹤每個 SNV 的 CLNSIG 分類
    # {snv_key: "Pathogenic" | "Benign" | "Other"}
    snv_clnsig_category = {}
    
    # 統計各 feature type 中的 CLNSIG 分布
    # {feature_type: Counter of CLNSIG categories}
    clnsig_by_feature = {
        'CDS': Counter(),
        '5UTR': Counter(),
        '3UTR': Counter()
    }
    
    # 統計 CDS 區域中，各 CLNSIG 分類下的 MC 分布
    # {clnsig_category: Counter of MC categories}
    mc_by_clnsig = {
        'Pathogenic': Counter(),
        'Benign': Counter(),
        'Other': Counter()
    }
    
    bed_region_count = 0
    
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
            
            unmatch_count = 0

            # 使用 fetch 找重疊的 SNV
            try:
                for record in vcf_reader.fetch(chrom, start, end):
                    pos = record.pos
                    ref = record.ref
                    alt = ",".join(record.alts) if record.alts else "."
                    
                    # 建立 SNV 唯一識別鍵
                    snv_key = f"{chrom}:{pos}:{ref}:{alt}"
                    
                    # 記錄此 SNV 配對到的 feature type (每次配對都記錄)
                    matched_snv_features[snv_key].append(feature_type)
                    
                    # 只在第一次遇到這個 SNV 時解析 MC 和 CLNSIG
                    if snv_key not in snv_clnsig_category:
                        # 解析並分類 MC
                        mc_field = record.info.get("MC", None)
                        mc_consequences = parse_mc_value(mc_field)
                        mc_category = categorize_mc(mc_consequences)
                        snv_mc_category[snv_key] = mc_category
                        
                        # 解析並分類 CLNSIG
                        clnsig_field = record.info.get("CLNSIG", None)
                        clnsig_category = categorize_clnsig(clnsig_field)
                        snv_clnsig_category[snv_key] = clnsig_category
                    
            except ValueError:
                unmatch_count += 1
                continue
    
    vcf_reader.close()
    
    print(f"已處理 {bed_region_count} 個 BED 區間")
    
    # 統計每個 SNV 的分類 (ㄧSNV多個分類的話重複計算)
    snv_by_feature = {
        'CDS': 0,
        '5UTR': 0,
        '3UTR': 0
    }

    total_matched = 0

    # 統計 CLNSIG 分布 (按 feature type) 和 CDS 區域的 MC 分布 (按 CLNSIG)
    for snv_key, features in matched_snv_features.items():
        mc_category = snv_mc_category.get(snv_key, "other")
        clnsig_category = snv_clnsig_category.get(snv_key, "Other")
        
        for feature_type in features:
            if feature_type in snv_by_feature:
                snv_by_feature[feature_type] += 1
                total_matched+=1
                
                # 統計該 feature type 的 CLNSIG 分類
                clnsig_by_feature[feature_type][clnsig_category] += 1
                
                # 只在 CDS 區域統計 MC 分布 (按 CLNSIG 分類)
                if feature_type == 'CDS':
                    mc_by_clnsig[clnsig_category][mc_category] += 1
    
    
    
    print(f'沒配對到的transcript數量:{unmatch_count}')
    # 輸出 SNV 分類統計結果
    print("\n" + "=" * 50)
    print("SNV 分類統計結果")
    print("=" * 50)
    print(f"\n{'Category':<15} {'Count':>10} {'Percentage':>12}")
    print("-" * 40)
    
    for category in ['CDS', '5UTR', '3UTR']:
        count = snv_by_feature[category]
        pct = (count / total_matched * 100) if total_matched > 0 else 0
        print(f"{category:<15} {count:>10} {pct:>11.2f}%")
    
    print("-" * 40)
    print(f"{'TOTAL':<15} {total_matched:>10} {'100.00%':>12}")
    
    # 輸出 CLNSIG 分布分析 (所有 feature types)
    analyze_clnsig_distribution(clnsig_by_feature)
    
    # 輸出 CDS 區域的 MC 分布 (按 CLNSIG 分類)
    analyze_mc_by_clnsig(mc_by_clnsig)


if __name__ == "__main__":
    main()
