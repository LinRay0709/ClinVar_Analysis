#!/usr/bin/env python3
"""
Analyze and categorize CLNSIG field values from ClinVar SNV VCF file.

This script reads a VCF file using pysam and categorizes variants into
three groups based on their CLNSIG values:
- Benign: includes Benign and Likely_benign
- Pathogenic: includes Pathogenic and Likely_pathogenic
- Other: all other values
"""

import pysam
from collections import Counter


# Define category mappings (case-insensitive matching)
PATHOGENIC_VALUES = {"pathogenic", "likely_pathogenic"}
BENIGN_VALUES = {"benign", "likely_benign"}


def categorize_clnsig(clnsig_value: str) -> str:
    """
    Categorize CLNSIG value into one of three groups.
    
    Handles:
    - Multiple values separated by '|'
    - Compound conditions separated by '/' (e.g., "Pathogenic/Likely_pathogenic")

    Args:
        clnsig_value: CLNSIG string value (may contain '|' and '/')

    Returns:
        Category string: "Pathogenic", "Benign", or "Other"
    """
    if not clnsig_value:
        return "Other"
    
    # First split by '|' to get individual values
    pipe_separated = clnsig_value.split('|')
    
    # Then split each value by '/' and flatten the list
    all_values = []
    for value in pipe_separated:
        slash_separated = value.split('/')
        all_values.extend(slash_separated)
    
    # Convert all values to lowercase for comparison
    values_lower = [v.strip().lower() for v in all_values]

    # Check for pathogenic first (prioritize pathogenic over benign if both exist)
    for value in values_lower:
        if value in PATHOGENIC_VALUES:
            return "Pathogenic"

    # Check for benign
    for value in values_lower:
        if value in BENIGN_VALUES:
            return "Benign"

    # Everything else is "Other"
    return "Other"


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
            category_counter["Other"] += 1
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
    print("CATEGORIZED RESULTS (3 groups)")
    print("=" * 60)
    print(f"\n{'Category':<20} {'Count':>10} {'Percentage':>12}")
    print("-" * 45)

    # Display in order: Pathogenic, Benign, Other
    for category in ["Pathogenic", "Benign", "Other"]:
        count = category_counts.get(category, 0)
        percentage = (count / total * 100) if total > 0 else 0
        print(f"{category:<20} {count:>10} {percentage:>11.2f}%")

    print("-" * 45)
    print(f"{'TOTAL':<20} {total:>10} {'100.00%':>12}")

    # Print category definitions
    print("\n" + "=" * 60)
    print("Category Definitions:")
    print("-" * 60)
    print("  Pathogenic : Pathogenic, Likely_pathogenic")
    print("  Benign     : Benign, Likely_benign")
    print("  Other      : All other values")

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
