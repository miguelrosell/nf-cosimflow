# nf-cosimflow

**A Nextflow subworkflow to calculate cosine similarity and generate heatmaps from gene expression data.**

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3eb049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

## 📖 Overview

**nf-cosimflow** was originally designed as a standalone pipeline but has been refactored into a modular **subworkflow** to comply with nf-core standards. 

It takes a gene expression matrix (CSV or XLSX) as input and performs the following:
1.  **Preprocessing:** Filters genes based on mean expression (optional).
2.  **Calculation:** Computes pairwise similarity between samples using Cosine, Pearson, or Spearman metrics.
3.  **Visualization:** Generates a high-quality heatmap (`.png`) and saves the similarity matrix (`.csv`).

This tool is designed to be easily integrated into larger pipelines like `nf-core/rnaseq` or `nf-core/differentialabundance` for exploratory data analysis (EDA) and QC.

## 🚀 Features

* **Modular Design:** Built as a local module (`compute_cosine`) and subworkflow (`run_cosimflow`).
* **Reproducible:** Fully containerized environment (supports Docker and Conda).
* **Flexible Inputs:** Accepts standard expression matrices (genes as rows, samples as columns).
* **Customizable:** Choose your metric (`cosine`, `pearson`, `spearman`) and filtering thresholds.

## 🛠️ Usage

You can run the workflow directly from this repository to test it:

```bash
nextflow run main.nf \
    --input data/test_expression.csv \
    -profile docker

```

### Parameters

| Parameter | Description | Default |
| --- | --- | --- |
| `--input` | Path to the expression matrix (CSV/XLSX). | `null` |
| `--method` | Similarity metric: `cosine`, `pearson`, or `spearman`. | `cosine` |
| `--min_gene_mean` | Filter out genes with mean count < value. | `0.0` |
| `--outdir` | Directory to save results. | `results` |

## 📂 Output

All results are saved in the `results/` directory:

* `*_matrix.csv`: The numeric similarity matrix.
* `*_heatmap.png`: A visual heatmap of the sample clustering.

## 🧩 Integration

To use this in your own Nextflow pipeline:

1. Copy `modules/local/compute_cosine.nf` and `subworkflows/local/run_cosimflow.nf` to your project.
2. Import it in your workflow:

```groovy
include { RUN_COSIMFLOW } from './subworkflows/local/run_cosimflow.nf'

workflow {
    // Create input channel: [ meta, file ]
    ch_input = Channel.fromPath(params.input).map { file -> [ [id: file.baseName], file ] }
    
    RUN_COSIMFLOW(ch_input)
}

```

## 👥 Credits

Developed by **Miguel Rosell** (@miguelrosell).
Feel free to open an issue if you find any bugs or have suggestions!
