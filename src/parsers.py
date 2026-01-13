# src/parsers.py

def parse_transcript_info(transcript_info: str) -> dict:
    """
    解析 BED 檔第四欄的資訊。
    假設格式範例: "NM_001005484.2|OR4F5|CDS" (需依實際 BED 內容調整)
    """
    parts = transcript_info.split('|')
    
    # 預設回傳值，避免 IndexOutOfRange 錯誤
    result = {
        "transcript_id": ".",
        "gene_name": ".",
        "feature_type": "."
    }
    
    if len(parts) >= 1:
        result["transcript_id"] = parts[0]
    if len(parts) >= 2:
        result["gene_name"] = parts[1]
    if len(parts) >= 3:
        result["feature_type"] = parts[2]  # 例如 CDS, 5UTR
        
    return result

def parse_mc_value(mc_field) -> list:
    """解析 VCF INFO 中的 MC 欄位"""
    consequences = []
    if mc_field is None:
        return consequences
    
    # 處理 tuple 或 string
    mc_values = [mc_field] if isinstance(mc_field, str) else list(mc_field)
    
    for mc in mc_values:
        # MC 格式通常為 "SO:000xxx|missense_variant"
        if '|' in mc:
            consequences.append(mc.split('|', 1)[1])
        else:
            consequences.append(mc)
    return consequences

def categorize_mc(mc_consequences: list) -> str:
    """將 MC 列表分類 (邏輯沿用舊程式)"""
    if not mc_consequences:
        return "other"
    
    has_missense = "missense_variant" in mc_consequences
    has_synonymous = "synonymous_variant" in mc_consequences
    
    if has_missense and has_synonymous:
        return "conflicting"
    if has_missense:
        return "missense_variant"
    if has_synonymous:
        return "synonymous_variant"
    return "other"

def categorize_clnsig(clnsig_field) -> str:
    """將 CLNSIG 分類 (邏輯沿用舊程式)"""
    if clnsig_field is None:
        return "other"
    
    # 處理 Tuple 並轉為字串
    clnsig_str = "|".join(clnsig_field) if isinstance(clnsig_field, tuple) else clnsig_field
    
    # 分割所有可能的值 ('|' 和 '/')
    pipe_separated = clnsig_str.split('|')
    all_values = []
    for value in pipe_separated:
        all_values.extend(value.split('/'))
        
    values_lower = [v.strip().lower() for v in all_values]
    
    # 優先順序判定
    if "pathogenic" in values_lower:
        return "pathogenic"
    if "benign" in values_lower:
        return "benign"
    if "likely_pathogenic" in values_lower:
        return "likely pathogenic"
    if "likely_benign" in values_lower:
        return "likely benign"
    
    return "other"