#!/usr/bin/env python3
"""
Analyze MC (Molecular Consequence) field values from ClinVar SNV VCF file.

This script reads a VCF file using pysam and counts all MC values,
displaying their distribution.
"""

import pysam
from collections import Counter


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


def analyze_mc(input_vcf: str) -> tuple[int, Counter, Counter]:
    """
    Analyze MC field values from a VCF file.

    Args:
        input_vcf: Path to input VCF file

    Returns:
        Tuple of (total_records, mc_counter, records_with_mc_count)
    """
    total_count = 0
    mc_counter = Counter()  # Count each MC value
    records_with_mc = 0
    records_without_mc = 0

    # Open input VCF file
    vcf_in = pysam.VariantFile(input_vcf, "r")

    # Iterate through all records
    for record in vcf_in:
        total_count += 1

        # Check if MC field exists
        if "MC" in record.info:
            mc_field = record.info["MC"]
            mc_consequences = parse_mc_value(mc_field)
            
            if mc_consequences:
                records_with_mc += 1
                for mc in mc_consequences:
                    mc_counter[mc] += 1
            else:
                records_without_mc += 1
                mc_counter["(empty)"] += 1
        else:
            records_without_mc += 1
            mc_counter["(missing)"] += 1

    vcf_in.close()

    return total_count, mc_counter, records_with_mc, records_without_mc


def main():
    # Define input file path
    input_file = "/home/czlin/ClinVar_project/data_new/clinvar_snv_new.vcf.gz"

    print(f"Analyzing MC values in: {input_file}")
    print("=" * 70)

    # Analyze MC values
    total, mc_counts, with_mc, without_mc = analyze_mc(input_file)

    # Print summary
    print(f"\nTotal records analyzed: {total}")
    print(f"Records with MC field: {with_mc}")
    print(f"Records without MC field: {without_mc}")

    # Print MC value distribution
    print("\n" + "=" * 70)
    print("MC (Molecular Consequence) 值分布")
    print("=" * 70)
    print(f"\nTotal unique MC values: {len(mc_counts)}")
    print(f"\n{'MC Value':<45} {'Count':>10} {'Percentage':>12}")
    print("-" * 70)

    total_mc_values = sum(mc_counts.values())
    
    # Sort by count (descending) and display
    for mc_type, count in mc_counts.most_common():
        pct = (count / total_mc_values * 100) if total_mc_values > 0 else 0
        # Truncate long MC values for display
        display_name = mc_type[:42] + "..." if len(mc_type) > 45 else mc_type
        print(f"{display_name:<45} {count:>10} {pct:>11.2f}%")

    print("-" * 70)
    print(f"{'TOTAL':<45} {total_mc_values:>10}")


if __name__ == "__main__":
    main()
