# src/models.py
from dataclasses import dataclass
from typing import Optional

@dataclass
class MatchedVariant:
    """
    定義一筆配對到的變異資料結構。
    """
    # SNV 基本資訊
    chrom: str
    pos: int
    ref: str
    alt: str
    
    # Transcript 資訊 (來自 BED)
    transcript_id: str
    feature_type: str  # CDS, 5UTR, 3UTR
    gene_name: str     # 如果 BED 裡有包含基因名稱
    
    # ClinVar 分析結果 (來自 VCF)
    clnsig_category: str
    mc_category: str
    
    # 原始資料 (用於除錯或後續詳細分析)
    clnsig_raw: Optional[str] = None
    mc_raw: Optional[str] = None

    def to_csv_row(self) -> list:
        """
        輔助函式：將物件轉為列表，方便寫入 CSV/TSV
        """
        return [
            self.chrom,
            str(self.pos),
            self.ref,
            self.alt,
            self.transcript_id,
            self.feature_type,
            self.gene_name,
            self.clnsig_category,
            self.mc_category,
            str(self.clnsig_raw),
            str(self.mc_raw)
        ]