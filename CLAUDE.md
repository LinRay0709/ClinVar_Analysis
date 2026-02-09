# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a bioinformatics pipeline for extracting Single Nucleotide Variants (SNVs) from VCF files and ClinVar database, then generating context sequences (±128bp/512bp) for deep learning model training. The project processes ~220K variants for pathogenicity prediction.

**Environment**: Use the `clinvar_env` conda environment.

## Core Pipeline Workflow

The pipeline follows a strict sequential order:

### 1. ETL & Variant Matching
```bash
python main.py
```
- Reads BED regions and VCF variants
- Matches variants with ClinVar annotations
- Output: `output/matched_variants.tsv`

### 2. Reference Consistency Check (Optional)
```bash
python src/extract_seq/ref_seq_consistency_check.py
```
- Validates VCF reference bases against hg38 FASTA files

### 3. Sequence Extraction
```bash
python src/extract_seq/extract_+-512_seq.py
```
- Extracts 1024bp sequences (512bp context on each side)
- Creates both reference and alternate sequences
- Output: `output/all_sequences.tsv`
- Format: snv_key, chrom, pos, ref, alt, feature_type, mc_category, clnsig_category, ref_seq, alt_seq

### 4. CD-HIT Preparation & Deduplication
```bash
python src/prepare_cdhit.py
```
- Converts TSV to FASTA format for CD-HIT
- Output: `output/interim/for_cdhit.fasta`

Run CD-HIT externally (command in prepare_cdhit.py), then:

### 5. Dataset Splitting
```bash
python src/split_datasets.py
```
- Parses CD-HIT clusters (.clstr file)
- Implements mutation-type based selection logic:
  - Pure label clusters: select one per mutation type
  - Mixed label clusters: discard conflicting types, keep consistent ones
- Splits into 80/10/10 (train/val/test) stratified by label
- Reconstructs ref_seq from alt_seq by swapping center base
- Output: `output/processed/{train,val,test}.csv`

### 6. RDDL Format Conversion
```bash
cd src
python bridge_to_rddl.py
```
- Converts CSV to RDDL format (pickle files)
- Creates directory structure: `RDDL/ClinVar/{USER_pos,USER_neg,RDDL_splitting_info}`
- Splits training set into 5 folds for cross-validation
- Encodes DNA sequences: A=0, C=1, G=2, T=3

### 7. Model Training
```bash
cd RDDL
python run_train.py -name ClinVar -m <METHOD> -g <GPU_NUM>
```
- Trains 5-fold cross-validation
- Methods: defined by balancing_methods() in training_func
- Requires: USER_model.py and hyperparameters JSON files
- Outputs: models, learning curves, ROC/PRC plots in `RDDL_outputs/`

### 8. Ensemble Prediction
```bash
python run_ensemble.py -name ClinVar
```
- Uses trained models for ensemble prediction on test set

## Key Architecture Principles

### Module Structure
- `src/config.py`: Centralized paths and constants
- `src/models.py`: Data classes (MatchedVariant, etc.)
- `src/parsers.py`: VCF/ClinVar annotation parsing logic
- `src/processing.py`: Core ETL using pysam for VCF queries
- `src/extract_seq/`: Sequence extraction subpackage
- `src/tests/`: Validation and verification scripts

### Sequence Handling Logic
- **1024bp mode** (current): 512bp context, center at index 512 (0-based)
- alt_seq has variant base at center, ref_seq has reference base
- Padding with '@' for sequences near chromosome boundaries
- Always verify center position matches expected base

### Data Deduplication Strategy
The project uses CD-HIT for sequence clustering, then applies custom logic:
- Groups sequences by mutation type (e.g., "A->G")
- Within each cluster, handles label conflicts:
  - If same mutation type has both pathogenic (1) and benign (0) labels → discard all
  - If labels are consistent → select one representative per mutation type
- This prevents data leakage and reduces redundancy

### RDDL Integration
- RDDL (Resampling-based Discrimination of Deep Learning) is the training framework
- Data must be converted to pickle format with specific directory structure
- Positive (pathogenic) → `USER_pos/`
- Negative (benign) → `USER_neg/`
- Split info → `RDDL_splitting_info/`

## Testing & Validation

Run tests from `src/tests/`:
- `check_paths.py`: Verify file paths exist
- `check_VCF_chr_type.py`: Validate chromosome naming (chr1 vs 1)
- `check_if_indel.py`: Confirm only SNVs (no insertions/deletions)
- `check_if_multiple_alt.py`: Check for multi-allelic variants
- `verify_output_logic.py`: Validate sequence lengths and center positions
- `validate_matched_tsv_loc.py`: Cross-reference coordinates
- `verify_rddl_data.py`: Check RDDL pickle files

## Important Notes

### File Paths
- Reference genome: `data/FASTA_data_from_HW1/` (hg38 chromosomes as symbolic links)
- Input data: `data_0107/` (BED) and `data_new/` (VCF)
- Output: `output/` directory (created automatically)

### Coordinate Systems
- VCF uses 1-based positions
- BED uses 0-based half-open intervals [start, end)
- FASTA access uses 0-based indexing
- Be careful when converting between systems

### Label Encoding
- Pathogenic variants: label = 1
- Benign variants: label = 0
- Stored in FASTA headers as `>snv_key|label`

### Performance Considerations
- VCF must be bgzipped and indexed (.tbi)
- Processing ~220K variants takes minutes
- CD-HIT reduces dataset from 220K to ~8K sequences
- Use progress indicators (already implemented) for long-running operations
