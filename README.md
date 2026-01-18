# Variant Sequence Extraction Pipeline (VCF to ML-Ready Dataset)

這是一個生物資訊資料處理專案，目標是從 VCF 變異檔與 ClinVar 資料庫中提取 Single Nucleotide Variants (SNV)，並擷取變異點前後的基因序列（Context Sequence），最終產生適用於深度學習模型訓練的結構化資料集。

env用clinvar_env

## 📁 專案架構 (Project Structure)

目前專案的檔案結構與功能說明如下：

project_root/
├── main.py                     # [Entry Point] 程式入口，負責執行 ETL 與資料配對流程
├── src/                        # 核心原始碼
│   ├── config.py               # 全域設定檔 (路徑、常數、欄位定義)
│   ├── models.py               # 資料模型定義 (Variant, MatchResult 等 Class)
│   ├── parsers.py              # 負責 VCF 與 ClinVar 資料的解析邏輯
│   ├── processing.py           # 核心 ETL 流程 (讀取、篩選、配對、寫入)
│   └── extract_seq/            # [Sub-package] 序列處理模組
│       ├── __init__.py
│       ├── ref_seq_consistency_check.py  # 檢查 Reference FASTA 與 VCF ref base 是否一致
│       └── extract_+-128_seq.py          # 正式程式：抓取前後 128bp 序列並整合標註資訊
├── tests/                      # 測試與驗證腳本
│   ├── check_paths.py                  # 檢查環境路徑設定
│   ├── check_VCF_chr_type.py           # 檢查 VCF 染色體命名格式
│   ├── check_if_indel.py               # 統計資料中是否包含 Indel (目前確認全為 SNV)
│   ├── check_if_multiple_alt.py        # 檢查是否存在多重等位基因 (Multi-allelic)
│   ├── verify_consistency_with_0121py.py # 與舊版程式邏輯進行一致性比對
│   └── verify_output_logic.py          # [最終驗證] 自動檢查輸出序列的長度與中心點正確性
└── data/                       # 資料存放區 (通常不納入 git)
    ├── reference/              # 存放 hg38 FASTA 檔案 (Symbolic Links)
    └── output/                 # 產出的 TSV 檔案

```

## 🚀 執行流程 (Workflow)

請依照以下順序執行程式以產出最終資料：

### 1. 環境準備

確保 `data/reference/` 資料夾下已建立 hg38 染色體序列的 Symbolic Links (`chr1.fa`, `chr2.fa`...)。

### 2. 資料清洗與配對 (ETL)

讀取原始 VCF 與 BED 檔案，進行篩選與 ClinVar 標註配對。

```bash
python main.py

```

* **產出**: `output/matched_variants.tsv`

### 3. 資料一致性檢查 (Optional)

在提取序列前，確保座標系統與 Reference Genome 版本一致。

```bash
python src/extract_seq/ref_seq_consistency_check.py

```

### 4. 序列提取 (Sequence Extraction)

根據配對結果，抓取變異點前後各 128bp 序列，並整合 Feature Type, MC, CLNSIG 等標籤。

```bash
python src/extract_seq/extract_+-128_seq.py

```

* **產出**: `output/all_sequences.tsv` (最終訓練資料)

### 5. 結果驗證

檢查最終輸出的序列邏輯（長度是否為 257bp、中心點是否匹配）。

```bash
python tests/verify_output_logic.py

```

---

## 📊 資料輸出格式

最終檔案 `all_sequences.tsv` 為 Tab 分隔文件，包含以下欄位：

| 欄位名稱 | 說明 |
| --- | --- |
| `snv_key` | 唯一識別碼 (Format: `chrom:pos:ref:alt`) |
| `chrom`, `pos` | 染色體與位置 (1-based) |
| `ref`, `alt` | 參考鹼基與變異鹼基 |
| `feature_type` | 基因特徵區域 (e.g., CDS, UTR) |
| `mc_category` | 分子後果 (e.g., missense_variant) |
| `clnsig_category` | 臨床意義 (e.g., pathogenic, benign) |
| `ref_seq` | 參考序列 (257bp, Center is Ref) |
| `alt_seq` | 變異序列 (257bp, Center is Alt) |

---

## 📝 開發進度 (Progress Log)

* [x] **Phase 1: 架構重構與 ETL**
* [x] 建立 `src/` 模組化架構
* [x] 完成 `processing.py` 處理 22 萬筆配對資料
* [x] 通過 `verify_consistency` 驗證


* [x] **Phase 2: 資料驗證**
* [x] 確認無 Multi-allelic 變異
* [x] 確認僅包含 SNV (無 Indel)
* [x] 建立 `check_paths` 與 `check_VCF` 工具


* [x] **Phase 3: 序列提取**
* [x] 建立 Reference Genome 連接 (`data/reference`)
* [x] 完成座標轉換與邊界處理 (Padding)
* [x] 整合標註資訊至序列檔案
* [x] 通過 UCSC Genome Browser 人工抽樣驗證




---

*Last Updated: 2026-01-17*