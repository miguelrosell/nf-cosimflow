# nf-cosimflow
Welcome to the workflow used to calculate cosine similarity in Nextflow

nf-cosimflow is a lightweight and reproducible Nextflow pipeline that computes cosine similarity between multiple biological samples based on their gene expression profiles. This pipeline helps bioinformaticians quickly assess sample similarity, reproducibility, and data consistency in transcriptomics or any other quantitative dataset.

🚀 Features

- Reads an expression matrix (e.g., RNA-seq normalized counts, proteomics, metabolomics…)
- Computes pairwise cosine similarity between samples
- Outputs both: a cosine similarity matrix (CSV) and a heatmap visualization (PNG)
- Fully reproducible and container-ready (Docker/Singularity)
- Ideal for data QC, batch comparison, and exploratory data analysis


🧩 Input format

Your input file must be a CSV file with this structure:

gene	sample1	sample2	sample3
TP53	5.2	4.9	5.1
BRCA1	3.4	3.7	3.5
EGFR	8.0	7.9	8.1
...	...	...	...


⚙️ How do i use it? 

1. Install Nextflow with: curl -s https://get.nextflow.io | bash mv nextflow ~/bin/
2. Clone this repository with: git clone https://github.com/<your-username>/nf-cosimflow.git cd nf-cosimflow
3. Run the pipeline with: nextflow run main.nf --input test_expression.csv
4. Results. After execution, all outputs are saved in the results/ directory: cosine_matrix.csv – pairwise cosine similarity and cosine_heatmap.png – heatmap of sample similarities

🧠 Parameters

You can adjust parameters in your nextflow.config file or directly in the command line: nextflow run main.nf --input my_data.csv --min_gene_mean 1.0
Default parameters (defined in nextflow.config):
params {
    input = 'test_expression.csv'
    min_gene_mean = 0.0
    sample_cols = []
}

🗂️ Output directory

All results are automatically copied into a folder outside work/ using the publishDir directive. You can find them here after each run: nf-cosimflow/results/

💡 Example command:

nextflow run main.nf -profile test -resume

Example output:

✅ Cosine matrix generated: results/cosine_matrix.csv

✅ Heatmap generated: results/cosine_heatmap.png

🧬 Why this pipeline?

Reproducible similarity analysis is crucial for: 

- Checking sample consistency in RNA-seq or proteomics.
- Detecting outliers or batch effects.
- Comparing replicates in large-scale projects.

nf-cosimflow brings these analyses into a modular Nextflow framework — easy to extend, easy to share, and 100% reproducible.

Citation & contact

Developed by Miguel Rosell Hidalgo.
For feedback or collaboration, contact: <miguelnorty@gmail.com>
Or open an issue on GitHub!
