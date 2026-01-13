#!/usr/bin/env python3
"""
Analyze and categorize CLNSIG field values from ClinVar SNV VCF file.

This script reads a VCF file using pysam and categorizes variants into
six groups based on their CLNSIG values:
- conflicting: contains both Benign and Pathogenic
- benign: includes Benign only
- likely benign: includes Likely_benign only
- pathogenic: includes Pathogenic only
- likely pathogenic: includes Likely_pathogenic only
- other: all other values
"""

import pysam
from collections import Counter


def categorize_clnsig(clnsig_value: str) -> str:
    """
    Categorize CLNSIG value into one of six groups.
    
    Handles:
    - Multiple values separated by '|'
    - Compound conditions separated by '/' (e.g., "Pathogenic/Likely_pathogenic")
    
    Priority order (if multiple labels exist):
    1. conflicting (if contains both "pathogenic" and "benign")
    2. pathogenic (if contains "pathogenic")
    3. benign (if contains "benign")
    4. likely pathogenic (if contains "likely_pathogenic")
    5. likely benign (if contains "likely_benign")
    6. other

    Args:
        clnsig_value: CLNSIG string value (may contain '|' and '/')

    Returns:
        Category string: "conflicting", "benign", "likely benign", "pathogenic", "likely pathogenic", or "other"
    """
    if not clnsig_value:
        return "other"
    
    # First split by '|' to get individual values
    pipe_separated = clnsig_value.split('|')
    
    # Then split each value by '/' and flatten the list
    all_values = []
    for value in pipe_separated:
        slash_separated = value.split('/')
        all_values.extend(slash_separated)
    
    # Convert all values to lowercase for comparison
    values_lower = [v.strip().lower() for v in all_values]

    # Check for conflicting first (highest priority): both benign and pathogenic present
    has_pathogenic = "pathogenic" in values_lower
    has_benign = "benign" in values_lower
    if has_pathogenic and has_benign:
        return "conflicting"

    # Check for pathogenic
    if has_pathogenic:
        return "pathogenic"
    
    # Check for benign
    if has_benign:
        return "benign"

    # Check for likely_pathogenic
    for value in values_lower:
        if value == "likely_pathogenic":
            return "likely pathogenic"

    # Check for likely_benign
    for value in values_lower:
        if value == "likely_benign":
            return "likely benign"

    # Everything else is "other"
    return "other"


def analyze_clnsig(input_vcf: str) -> tuple[int, Counter, Counter]:
    """
    Analyze and categorize CLNSIG field values from a VCF file.

    Args:
        input_vcf: Path to input VCF file

    Returns:
        Tuple of (total_records, category_counter, original_value_counter)
    """
    total_count = 0
    category_counter = Counter()
    original_counter = Counter()  # Keep track of original values for reference

    # Open input VCF file
    vcf_in = pysam.VariantFile(input_vcf, "r")

    # Iterate through all records
    for record in vcf_in:
        total_count += 1

        # Check if CLNSIG field exists
        if "CLNSIG" in record.info:
            clnsig = record.info["CLNSIG"]
            
            # Check if CLNSIG is a tuple (multiple values)
            if isinstance(clnsig, tuple):
                # For original_counter: record each value separately
                for value in clnsig:
                    original_counter[value] += 1
                # For categorization: join with '|' to process as single string
                clnsig_str = "|".join(clnsig)
            else:
                # Single string value
                original_counter[clnsig] += 1
                clnsig_str = clnsig

            # Categorize the variant (splits by '|' and '/')
            category = categorize_clnsig(clnsig_str)
            category_counter[category] += 1
        else:
            category_counter["other"] += 1
            original_counter["(missing)"] += 1

    vcf_in.close()

    return total_count, category_counter, original_counter


def main():
    # Define input file path
    input_file = "/home/czlin/ClinVar_project/data_new/clinvar_snv_new.vcf.gz"

    print(f"Analyzing and categorizing CLNSIG values in: {input_file}")
    print("=" * 60)

    # Analyze CLNSIG values
    total, category_counts, original_counts = analyze_clnsig(input_file)

    # Print categorized results
    print(f"\nTotal variants analyzed: {total}")
    print("\n" + "=" * 60)
    print("CATEGORIZED RESULTS (6 groups)")
    print("=" * 60)
    print(f"\n{'Category':<20} {'Count':>10} {'Percentage':>12}")
    print("-" * 45)

    # Display in order: pathogenic, likely pathogenic, benign, likely benign, conflicting, other
    for category in ["pathogenic", "likely pathogenic", "benign", "likely benign", "conflicting", "other"]:
        count = category_counts.get(category, 0)
        percentage = (count / total * 100) if total > 0 else 0
        print(f"{category:<20} {count:>10} {percentage:>11.2f}%")

    print("-" * 45)
    print(f"{'TOTAL':<20} {total:>10} {'100.00%':>12}")

    # Print category definitions
    print("\n" + "=" * 60)
    print("Category Definitions:")
    print("-" * 60)
    print("  pathogenic        : Pathogenic only")
    print("  likely pathogenic : Likely_pathogenic only")
    print("  benign            : Benign only")
    print("  likely benign     : Likely_benign only")
    print("  conflicting       : Contains both Benign and Pathogenic")
    print("  other             : All other values")

    # Print original value breakdown for reference
    print("\n" + "=" * 60)
    print("ORIGINAL CLNSIG VALUES (for reference)")
    print("=" * 60)
    print(f"\n{'Original Value':<40} {'Count':>10}")
    print("-" * 55)

    for value, count in original_counts.most_common():
        print(f"{value:<40} {count:>10}")

    print("-" * 55)


if __name__ == "__main__":
    main()
