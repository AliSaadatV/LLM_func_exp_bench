# LLM PS3/BS3 Functional Evidence Benchmark

Benchmarking large language models (LLMs) for extracting ACMG/AMP PS3/BS3 functional evidence from scientific literature. This repository evaluates `gpt-4o-mini` and `o4-mini` on abstract and full-text PDF classification tasks using ClinGen-curated variants as ground truth.

## Overview

This benchmark evaluates LLMs on four tasks:
1. **Abstract-level functional experiment detection** (binary classification)
2. **PS3/BS3 evidence level classification** from full-text PDFs
3. **Evidence strength classification** (4-class: supporting, moderate, strong, very_strong)
4. **Joint level+strength classification** (8-class)

## Quick Results

### PS3/BS3 Binary Classification (Full-Text)
- **o4-mini**: 96.3% accuracy, F1=0.979, 91.6% coverage
- **gpt-4o-mini**: 92.6% accuracy, F1=0.960, 99.4% coverage

### Abstract Classification
- **o4-mini**: 76.7% accuracy, F1=0.791
- **gpt-4o-mini**: 74.7% accuracy, F1=0.781

*Note: Evidence strength classification is more challenging (~34% accuracy for both models).*

## Setup

### Prerequisites
- Python 3.12
- Conda (recommended)

### Installation

1. Clone the repository:
```bash
git clone <repository-url>
cd llm-ps3-bs3-functional-evidence-benchmark
```

2. Create conda environment:
```bash
conda env create -f requirements.yml
conda activate ai
```

3. Install additional dependencies:
```bash
pip install pydantic tenacity tqdm seaborn scikit-learn numpy
```

4. Create `.env` file with your API keys:
```
OPENAI_API_KEY=your-key-here
NCBI_API_KEY=your-key-here
NCBI_EMAIL=your-email-here
```

## Usage

### Execution Order

Run notebooks in sequence:

1. **Preprocess ClinGen data** (`notebooks/01_preprocess_clingen.ipynb`)
   - Extracts PS3/BS3 evidence from ClinGen summaries
   - Fetches PubMed abstracts
   - Downloads PDFs to `res/pdfs/`
   - Output: `data/ps3_bs3_df_processed.csv`

2. **Abstract classification benchmark** (`notebooks/02_abstract_class_bench.ipynb`)
   - Binary classification: functional experiment detection
   - Output: `data/abstract_class_bench_functional_labels.csv`

3. **Full-text classification benchmark** (`notebooks/03_full_text_class_bench.ipynb`)
   - PS3/BS3 level, strength, and joint classification
   - Variant matching evaluation
   - LLM-as-judge coverage scoring
   - Output: `data/ps3_bs3_df_processed_withVariantIDs_llmClassifier_llmJudgeScored.csv`


## Data

### Input Files
- `data/clingen_curated_variants.txt`: ClinGen-curated variants with ACMG/AMP evidence codes
- `data/cgbench_pubmed_id_to_text.csv`: PubMed abstracts for negative sampling

### Output Files
- `data/ps3_bs3_df_processed.csv`: Processed ClinGen data with extracted PMIDs and abstracts
- `data/ps3_bs3_df_processed_withVariantIDs_llmClassifier_llmJudgeScored.csv`: Final results with LLM predictions and judge scores
- `res/figures/*.png`: Evaluation plots and visualizations

## Models Evaluated

- **gpt-4o-mini** (OpenAI)
- **o4-mini** (OpenAI)
- **gpt-4.1** (OpenAI, used as LLM judge)

## Citation

If you use this benchmark, please cite:

```bibtex
@software{llm_ps3_bs3_benchmark,
  title = {LLM PS3/BS3 Functional Evidence Benchmark},
  author = {[Your Name]},
  year = {2025},
  url = {https://github.com/[username]/llm-ps3-bs3-functional-evidence-benchmark}
}
```
