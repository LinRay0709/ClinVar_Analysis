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


def analyze_clnsig_distribution(clnsig_by_feature: dict):
    """
    分析並顯示各 feature type 中 CLNSIG 的分布情況 (5 類)。
    
    Args:
        clnsig_by_feature: {feature_type: Counter of CLNSIG categories}
    """
    print("\n" + "=" * 70)
    print("CLNSIG 分布分析 (5 類)")
    print("=" * 70)
    print("\n分類定義:")
    print("  pathogenic        : Pathogenic (優先)")
    print("  benign            : Benign (優先)")
    print("  likely pathogenic : Likely_pathogenic")
    print("  likely benign     : Likely_benign")
    print("  other             : 所有其他值")
    
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
        
        # 按固定順序顯示: pathogenic, likely pathogenic, benign, likely benign, other
        for category in ['pathogenic', 'likely pathogenic', 'benign', 'likely benign', 'other']:
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
    bed_file = "/home/czlin/ClinVar_project/data_0107/MANE_transcript.bed"
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
    
    # 追蹤每個 SNV 的 MC 原始值
    # {snv_key: mc_field_raw}
    snv_mc_raw = {}
    
    # 追蹤每個 SNV 的 CLNSIG 分類
    # {snv_key: "pathogenic" | "likely pathogenic" | "benign" | "likely benign" | "other"}
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
        'pathogenic': Counter(),
        'likely pathogenic': Counter(),
        'benign': Counter(),
        'likely benign': Counter(),
        'other': Counter()
    }
    
    # 記錄前三個配對供驗證使用
    # 格式: [(snv_key, bed_region_str, feature_type), ...]
    first_three_matches = []
    
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
            
            # 使用 fetch 找重疊的 SNV
            try:
                for record in vcf_reader.fetch(chrom, start, end):
                    pos = record.pos
                    ref = record.ref
                    alt = ",".join(record.alts) if record.alts else "."
                    
                    # 建立 SNV 唯一識別鍵
                    snv_key = f"{chrom}:{pos}:{ref}:{alt}"
                    
                    # 記錄前三個配對供驗證
                    if len(first_three_matches) < 3:
                        bed_region_str = f"{chrom}:{start}-{end}"
                        first_three_matches.append((snv_key, bed_region_str, feature_type, transcript_info))
                    
                    # 記錄此 SNV 配對到的 feature type (每次配對都記錄)
                    matched_snv_features[snv_key].append(feature_type)
                    
                    # 只在第一次遇到這個 SNV 時解析 MC 和 CLNSIG
                    if snv_key not in snv_clnsig_category:
                        # 解析並分類 MC
                        mc_field = record.info.get("MC", None)
                        mc_consequences = parse_mc_value(mc_field)
                        mc_category = categorize_mc(mc_consequences)
                        snv_mc_category[snv_key] = mc_category
                        snv_mc_raw[snv_key] = mc_field  # 儲存原始 MC 值
                        
                        # 解析並分類 CLNSIG
                        clnsig_field = record.info.get("CLNSIG", None)
                        clnsig_category = categorize_clnsig(clnsig_field)
                        snv_clnsig_category[snv_key] = clnsig_category
                    
            except ValueError:
                continue
    
    vcf_reader.close()
    
    print(f"已處理 {bed_region_count} 個 BED 區間")

    # ==========================================
    # [新增功能] 檢查與統計衝突情況
    # ==========================================
    print("\n" + "=" * 60)
    print("正在執行 Feature Type 衝突檢查 (Data Integrity Check)...")
    print("=" * 60)
    
    total_snvs = len(matched_snv_features)
    conflict_count = 0
    conflict_examples = []
    
    # 遍歷所有配對到的 SNV
    for snv_key, features in matched_snv_features.items():
        # 使用 set 去除重複 (例如 ['CDS', 'CDS'] -> {'CDS'} 是合法的)
        unique_features = set(features)
        
        # 如果 set 長度 > 1，代表同時被歸類為不同的區域 (例如 {'CDS', '5UTR'})
        if len(unique_features) > 1:
            conflict_count += 1
            if len(conflict_examples) < 10: # 只記錄前 10 個範例
                conflict_examples.append((snv_key, unique_features))
    
    # 輸出檢查報告
    print(f"總配對 SNV 數量: {total_snvs}")
    
    if total_snvs > 0:
        conflict_rate = (conflict_count / total_snvs) * 100
        print(f"發現衝突 SNV 數量: {conflict_count} ({conflict_rate:.4f}%)")
    else:
        print("發現衝突 SNV 數量: 0")

    if conflict_count > 0:
        print("\n[警告] 以下 SNV 同時屬於多個互斥區域 (前 10 筆範例):")
        for snv, feats in conflict_examples:
            print(f"  - {snv}: {feats}")
        print("\n提示: 這些 SNV 目前仍在後續統計中，根據您原本的邏輯，它們會被重複計算。")
    else:
        print("\n[完美] 檢查通過！沒有任何 SNV 同時屬於 CDS/UTR 衝突區域。")
        
    print("=" * 60 + "\n")
    
    # 輸出前三個配對供驗證
    print("\n" + "=" * 80)
    print("配對正確性驗證 (前 3 個配對)")
    print("=" * 80)
    for i, (snv_key, bed_region, feature_type, transcript_info) in enumerate(first_three_matches, 1):
        print(f"\n配對 {i}:")
        print(f"  SNV 座標      : {snv_key}")
        print(f"  BED 區段座標 : {bed_region}")
        print(f"  Feature Type : {feature_type}")
        print(f"  Transcript   : {transcript_info}")
    
    # 統計每個 SNV 的分類 (ㄧSNV多個分類的話重複計算)
    snv_by_feature = {
        'CDS': 0,
        '5UTR': 0,
        '3UTR': 0
    }
    
    # 統計 CLNSIG 分布 (按 feature type) 和 CDS 區域的 MC 分布 (按 CLNSIG)
    # 同時收集 CDS + pathogenic + MC other 的變異
    cds_pathogenic_mc_other = []  # [(snv_key, mc_raw), ...]
    
    for snv_key, features in matched_snv_features.items():
        mc_category = snv_mc_category.get(snv_key, "other")
        clnsig_category = snv_clnsig_category.get(snv_key, "other")
        
        for feature_type in features:
            if feature_type in snv_by_feature:
                snv_by_feature[feature_type] += 1
                
                # 統計該 feature type 的 CLNSIG 分類
                clnsig_by_feature[feature_type][clnsig_category] += 1
                
                # 只在 CDS 區域統計 MC 分布 (按 CLNSIG 分類)
                if feature_type == 'CDS':
                    mc_by_clnsig[clnsig_category][mc_category] += 1
                    
                    # 收集 CDS + pathogenic + MC other 的變異
                    if clnsig_category == 'pathogenic' and mc_category == 'other':
                        mc_raw = snv_mc_raw.get(snv_key)
                        mc_raw_str = str(mc_raw) if mc_raw else '(missing)'
                        cds_pathogenic_mc_other.append((snv_key, mc_raw_str))
    
    total_matched = len(matched_snv_features)
    
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
    
    # 輸出 CDS + pathogenic + MC other 的變異到 txt 檔案
    output_file = "/home/czlin/ClinVar_project/data_new/cds_pathogenic_mc_other.txt"
    
    # 去重 (同一個 SNV 可能配對到多個 CDS 區段)
    unique_records = list(set(cds_pathogenic_mc_other))
    
    with open(output_file, 'w') as f:
        # 寫入標題行
        f.write("SNV_Key\tMC_Raw\n")
        
        # 寫入每個變異
        for snv_key, mc_raw in unique_records:
            f.write(f"{snv_key}\t{mc_raw}\n")
    
    print("\n" + "=" * 70)
    print(f"CDS 區域 pathogenic 且 MC 為 other 的變異已輸出至:")
    print(f"  {output_file}")
    print(f"  共 {len(unique_records)} 筆不重複記錄")
    print("=" * 70)


if __name__ == "__main__":
    main()
