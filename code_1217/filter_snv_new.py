#!/usr/bin/env python3
"""
Filter SNV (Single Nucleotide Variants) from ClinVar VCF file.

This script reads a clinvar.vcf.gz file using pysam, filters records where
INFO field CLNVC == "single_nucleotide_variant", and outputs matching records
to a new VCF file.
"""

import pysam


def filter_snv(input_vcf: str, output_vcf: str) -> tuple[int, int]:
    """
    Filter SNV records from a VCF file.

    Args:
        input_vcf: Path to input VCF file (can be .vcf or .vcf.gz)
        output_vcf: Path to output VCF file

    Returns:
        Tuple of (total_records, snv_records)
    """
    total_count = 0
    snv_count = 0

    # Open input VCF file
    vcf_in = pysam.VariantFile(input_vcf, "r")

    # Open output VCF file with the same header as input
    vcf_out = pysam.VariantFile(output_vcf, "w", header=vcf_in.header)

    # Iterate through all records
    for record in vcf_in:
        total_count += 1

        # Check if CLNVC field exists and equals "single_nucleotide_variant"
        if "CLNVC" in record.info:
            clnvc = record.info["CLNVC"]
            if clnvc == "single_nucleotide_variant":
                snv_count += 1
                vcf_out.write(record)

    # Close files
    vcf_in.close()
    vcf_out.close()

    return total_count, snv_count


def main():
    # Define input and output file paths
    input_file = "/home/czlin/ClinVar_project/data/clinvar.vcf.gz"
    output_file = "/home/czlin/ClinVar_project/data_new/clinvar_snv_new.vcf.gz"

    print(f"Processing: {input_file}")
    print(f"Output to: {output_file}")
    print("-" * 50)

    # Filter SNV records
    total, snv = filter_snv(input_file, output_file)

    # Print results
    print(f"Total records processed: {total}")
    print(f"SNV records found: {snv}")
    print(f"Filtering complete!")


if __name__ == "__main__":
    main()
