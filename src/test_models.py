# src/test_models.py
from models import MatchedVariant

# 1. 模擬一筆測試資料 (假裝這是從 VCF 讀到的)
print("--- 測試 1: 建立完整的變異資料 ---")
variant = MatchedVariant(
    chrom="chr1",
    pos=123456,
    ref="A",
    alt="G",
    transcript_id="NM_001.2",
    feature_type="CDS",
    gene_name="TEST_GENE",
    clnsig_category="pathogenic",
    mc_category="missense_variant",
    # 這裡我們傳入原始資料
    clnsig_raw="Pathogenic/Likely_pathogenic",
    mc_raw="SO:0001583|missense_variant"
)

# 2. 測試 dataclass 自動產生的顯示功能 (__repr__)
print("物件內容 (自動格式化):")
print(variant) 
# 如果你看到像 MatchedVariant(chrom='chr1'...) 這樣的輸出，代表 dataclass 成功了

# 3. 測試轉 CSV 功能
print("\n轉換為 CSV 列表:")
csv_row = variant.to_csv_row()
print(csv_row)

# 4. 測試 Optional 功能 (故意不傳入 raw data)
print("\n--- 測試 2: 建立缺漏資料 (測試 Optional) ---")
variant_missing = MatchedVariant(
    chrom="chr2",
    pos=999,
    ref="T",
    alt="C",
    transcript_id="NM_002.1",
    feature_type="3UTR",
    gene_name="GENE_B",
    clnsig_category="other",
    mc_category="other"
    # 注意：這裡我們故意不傳入 clnsig_raw 和 mc_raw
)

print("缺漏資料物件:")
print(variant_missing)
print(f"檢查預設值: clnsig_raw = {variant_missing.clnsig_raw}") 
# 預期應該要印出 None