#!/usr/bin/env python3
"""
分析只存在於 5'UTR/3'UTR 區域的 SNV 變異（排除同時存在於 CDS 的變異）。

此腳本的目的是：
1. 找出配對到 5'UTR 或 3'UTR 的 SNV
2. 排除同時也配對到 CDS 的 SNV
3. 分析這些 SNV 的 CLNSIG 良惡性分布
4. 分析這些 SNV 的 MC (Molecular Consequence) 分布
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
        print(f'not tuple: {mc_field}')
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
    print("CLNSIG 分布分析 (5 類) - 僅 UTR 區域 (排除 CDS)")
    print("=" * 70)
    print("\n分類定義:")
    print("  pathogenic        : Pathogenic (優先)")
    print("  benign            : Benign (優先)")
    print("  likely pathogenic : Likely_pathogenic")
    print("  likely benign     : Likely_benign")
    print("  other             : 所有其他值")
    
    for feature_type in ['5UTR', '3UTR']:
        clnsig_counter = clnsig_by_feature.get(feature_type, Counter())
        total = sum(clnsig_counter.values())
        
        print(f"\n【{feature_type}】區域 - 僅 UTR (共 {total} 個 SNV)")
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


def analyze_mc_by_clnsig(mc_by_clnsig: dict, feature_type: str):
    """
    分析並顯示 UTR 區域中，各 CLNSIG 分類下的 MC 分布情況。
    
    Args:
        mc_by_clnsig: {clnsig_category: Counter of MC categories}
        feature_type: 區域類型 (5UTR 或 3UTR)
    """
    print("\n" + "=" * 70)
    print(f"{feature_type} 區域 - 各 CLNSIG 分類下的 MC 分布 (排除 CDS)")
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
    bed_file = "/home/czlin/ClinVar_project/data_new/all_transcript.bed"
    vcf_file = "/home/czlin/ClinVar_project/data_new/clinvar_snv_new.vcf.gz"
    
    print(f"分析僅存在於 5'UTR/3'UTR 的 SNV (排除同時存在於 CDS 的變異)")
    print(f"BED 檔案: {bed_file}")
    print(f"VCF 檔案: {vcf_file}")
    print("=" * 70)
    
    # 開啟 VCF 檔案 (必須有 .tbi 索引檔在同目錄下)
    try:
        vcf_reader = pysam.VariantFile(vcf_file)
    except Exception as e:
        print(f"錯誤: 無法讀取 VCF，請確認是否已壓縮並建立索引 (.tbi)。\n{e}")
        return

    # 第一輪：收集所有 SNV 配對到的 feature types
    # {snv_key: set of feature_types}
    snv_all_features = defaultdict(set)
    
    # 追蹤每個 SNV 的 MC 和 CLNSIG (原始值，稍後分類)
    snv_mc_raw = {}
    snv_clnsig_raw = {}
    
    # 記錄 SNV 配對到的 BED 區段資訊 (用於驗證)
    # {snv_key: [(bed_region_str, feature_type, transcript_info), ...]}
    snv_bed_info = defaultdict(list)
    
    bed_region_count = 0
    
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
            
            # 使用 fetch 找重疊的 SNV
            try:
                for record in vcf_reader.fetch(chrom, start, end):
                    pos = record.pos
                    ref = record.ref
                    alt = ",".join(record.alts) if record.alts else "."
                    
                    # 建立 SNV 唯一識別鍵
                    snv_key = f"{chrom}:{pos}:{ref}:{alt}"
                    
                    # 記錄此 SNV 配對到的 feature type (用 set 去重)
                    snv_all_features[snv_key].add(feature_type)
                    
                    # 記錄 BED 區段資訊
                    bed_region_str = f"{chrom}:{start}-{end}"
                    snv_bed_info[snv_key].append((bed_region_str, feature_type, transcript_info))
                    
                    # 只在第一次遇到這個 SNV 時記錄 MC 和 CLNSIG
                    if snv_key not in snv_mc_raw:
                        snv_mc_raw[snv_key] = record.info.get("MC", None)
                        snv_clnsig_raw[snv_key] = record.info.get("CLNSIG", None)
                    
            except ValueError:
                continue
    
    vcf_reader.close()
    
    print(f"已處理 {bed_region_count} 個 BED 區間")
    print(f"總共找到 {len(snv_all_features)} 個不重複 SNV")
    
    # 第二輪：篩選只存在於 UTR（不含 CDS）的 SNV
    utr_only_snvs = {
        '5UTR': [],  # 只在 5UTR (可能也在 3UTR，但不在 CDS)
        '3UTR': []   # 只在 3UTR (可能也在 5UTR，但不在 CDS)
    }
    
    for snv_key, features in snv_all_features.items():
        # 排除同時存在於 CDS 的 SNV
        if 'CDS' in features:
            continue
        
        # 根據 feature types 分類
        if '5UTR' in features:
            utr_only_snvs['5UTR'].append(snv_key)
        if '3UTR' in features:
            utr_only_snvs['3UTR'].append(snv_key)
    
    # 計算僅 UTR 的 SNV 數量 (去重)
    all_utr_only = set(utr_only_snvs['5UTR']) | set(utr_only_snvs['3UTR'])
    
    print(f"\n篩選結果:")
    print(f"  僅 5'UTR 的 SNV: {len(utr_only_snvs['5UTR'])} 個")
    print(f"  僅 3'UTR 的 SNV: {len(utr_only_snvs['3UTR'])} 個")
    print(f"  僅 UTR (不含 CDS) 的不重複 SNV 總數: {len(all_utr_only)} 個")
    
    # 輸出前三個配對供驗證
    print("\n" + "=" * 80)
    print("配對正確性驗證 (前 3 個僅 UTR 的配對)")
    print("=" * 80)
    
    validation_count = 0
    for snv_key in list(all_utr_only)[:3]:
        validation_count += 1
        features = snv_all_features[snv_key]
        bed_infos = snv_bed_info[snv_key]
        
        print(f"\n配對 {validation_count}:")
        print(f"  SNV 座標       : {snv_key}")
        print(f"  Feature Types  : {', '.join(features)}")
        print(f"  配對的 BED 區段:")
        for bed_region, feat_type, transcript in bed_infos[:2]:  # 最多顯示 2 個
            print(f"    - {bed_region} ({feat_type}) - {transcript}")
        if len(bed_infos) > 2:
            print(f"    ... 還有 {len(bed_infos) - 2} 個配對")
    
    # 統計 CLNSIG 分布 (按 feature type)
    clnsig_by_feature = {
        '5UTR': Counter(),
        '3UTR': Counter()
    }
    
    # 統計 MC 分布 (按 CLNSIG 分類，分別對 5UTR 和 3UTR)
    mc_by_clnsig_5utr = {
        'pathogenic': Counter(),
        'likely pathogenic': Counter(),
        'benign': Counter(),
        'likely benign': Counter(),
        'other': Counter()
    }
    
    mc_by_clnsig_3utr = {
        'pathogenic': Counter(),
        'likely pathogenic': Counter(),
        'benign': Counter(),
        'likely benign': Counter(),
        'other': Counter()
    }
    
    # 收集含有 missense/synonymous MC 的變異
    # 格式: [(snv_key, feature_type, clnsig_category, mc_category, mc_raw), ...]
    missense_synonymous_snvs = []
    
    # 分析每個僅 UTR 的 SNV
    for feature_type in ['5UTR', '3UTR']:
        for snv_key in utr_only_snvs[feature_type]:
            # 分類 CLNSIG
            clnsig_category = categorize_clnsig(snv_clnsig_raw.get(snv_key))
            clnsig_by_feature[feature_type][clnsig_category] += 1
            
            # 分類 MC
            mc_consequences = parse_mc_value(snv_mc_raw.get(snv_key))
            mc_category = categorize_mc(mc_consequences)
            
            if feature_type == '5UTR':
                mc_by_clnsig_5utr[clnsig_category][mc_category] += 1
            else:
                mc_by_clnsig_3utr[clnsig_category][mc_category] += 1
            
            # 收集含有 missense/synonymous 的變異
            if mc_category in ['missense_variant', 'synonymous_variant', 'conflicting']:
                mc_raw_str = str(snv_mc_raw.get(snv_key, ''))
                missense_synonymous_snvs.append((snv_key, feature_type, clnsig_category, mc_category, mc_raw_str))
    
    # 輸出統計結果
    print("\n" + "=" * 70)
    print("SNV 統計結果 (僅 UTR，排除 CDS)")
    print("=" * 70)
    print(f"\n{'Feature Type':<15} {'Count':>10} {'Percentage':>12}")
    print("-" * 40)
    
    total_utr = len(all_utr_only)
    for feature_type in ['5UTR', '3UTR']:
        count = len(utr_only_snvs[feature_type])
        pct = (count / total_utr * 100) if total_utr > 0 else 0
        print(f"{feature_type:<15} {count:>10} {pct:>11.2f}%")
    
    print("-" * 40)
    print(f"{'TOTAL (unique)':<15} {total_utr:>10}")
    print("\n註: 總數為不重複 SNV 數，單一 SNV 可能同時存在於 5'UTR 和 3'UTR")
    
    # 輸出 CLNSIG 分布分析
    analyze_clnsig_distribution(clnsig_by_feature)
    
    # 輸出 MC 分布分析 (按 CLNSIG 分類)
    analyze_mc_by_clnsig(mc_by_clnsig_5utr, "5UTR")
    analyze_mc_by_clnsig(mc_by_clnsig_3utr, "3UTR")
    
    # 輸出含有 missense/synonymous MC 的變異到 txt 檔案
    output_file = "/home/czlin/ClinVar_project/data_new/utr_missense_synonymous_snvs.txt"
    with open(output_file, 'w') as f:
        # 寫入標題行
        f.write("SNV_Key\tFeature_Type\tCLNSIG_Category\tMC_Category\tMC_Raw\n")
        
        # 寫入每個變異
        for snv_key, feature_type, clnsig_category, mc_category, mc_raw in missense_synonymous_snvs:
            f.write(f"{snv_key}\t{feature_type}\t{clnsig_category}\t{mc_category}\t{mc_raw}\n")
    
    print("\n" + "=" * 70)
    print(f"含有 missense/synonymous MC 的變異已輸出至:")
    print(f"  {output_file}")
    print(f"  共 {len(missense_synonymous_snvs)} 筆記錄")
    print("=" * 70)


if __name__ == "__main__":
    main()
