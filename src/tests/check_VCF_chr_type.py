# src/check_vcf_header.py
import sys
import os

dir_path = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.abspath(os.path.join(dir_path, ".."))
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

import pysam
import config


try:
    vcf = pysam.VariantFile(config.VCF_FILE)
    print("VCF 檔案中的前 10 個染色體名稱:")
    for i, contig in enumerate(vcf.header.contigs):
        if i >= 10: break
        print(f"  - {contig}")
except Exception as e:
    print(e)