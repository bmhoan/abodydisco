
# IPIAbDiscov: A Python Package, developped at IPI, for NGS-Based Antibody Discovery
IPIAbDiscov is an open-source Python package designed to streamline the analysis of Next-Generation Sequencing (NGS) data from antibody display technologies, such as phage display or yeast display libraries. It provides a command-line interface for processing raw FASTQ files from antibody selection campaigns, enabling researchers to quantify sequence abundance, track enrichment across selection rounds, and identify promising lead candidates for therapeutic antibody development.

IPIAbDiscov allow to accelerate the transition from huge amount of NGS data to validated antibody leads, reducing reliance on traditional low-throughput screening while uncovering rare, high-potential clones often missed in conventional pipelines. Ideal for academic and biotech researchers in antibody discovery engineering.

# Key Features

FASTQ Processing: Quality trimming, filtering, and annotation of antibody sequences using tools like fastp and ANARCI for accurate numbering and germline assignment of Fab or scFv fragments.
Data Aggregation: Combine results from multiple samples or rounds to generate comprehensive repertoire tables.
Lead Selection: Automated ranking and selection of top-enriched sequences based on read counts and enrichment metrics.
Repeatability Checks: Identify and quantify repeated or duplicated sequences across datasets.

The package is lightweight, dependency-managed (via requirements.txt and Bioconda tools), and configured through YAML files and sample sheets, making it suitable for MiSeq or similar NGS runs in antibody discovery projects.
Potential Future Enhancements

To further empower antibody discovery workflows, future versions could include:

Advanced Visualization: Interactive plots for repertoire diversity (e.g., Shannon entropy, clonal frequency distributions), enrichment heatmaps across selection rounds, sequence logos for CDR regions, and pairwise similarity networks.

Fold Change and Enrichment Analysis: Statistical computation of log-fold changes between pre- and post-selection rounds, with significance testing to highlight antigen-specific binders.

Machine Learning Integration: Predictive models for affinity or developability scoring, clustering of related sequences (e.g., lineage grouping), or epitope binning using sequence features.

Developability Assessment: Built-in filters for liabilities (e.g., glycosylation sites, cysteine residues) and integration with external databases for off-target prediction.


Repository: https://github.com/bmhoan/IPIAbDiscov




# package installtion

download abodydisco from ipi githup

#install requirements
pip install -r requirements.txt

#fastp instalation
conda install -c bioconda fastp

#anarci installation
https://github.com/oxpig/ANARCI.git
conda install -c conda-forge biopython -y
conda install -c bioconda hmmer=3.3.2 -y
cd ANARCI
python setup.py install

# How To Use AbodyDisco

cd abodydisco

python __main__.py process \
  --config config.yaml \
  --sample-sheet /Users/Hoan.Nguyen/ComBio/NGS/Projects/AntibodyDiscovery/Miseq104/Fastq/Miseq104_SampleSheet.xlsx \
  --fastq-folder /Users/Hoan.Nguyen/ComBio/NGS/Projects/AntibodyDiscovery/Miseq104/Fastq

python  __main__.py combine --folder /Users/Hoan.Nguyen/ComBio/NGS/Projects/AntibodyDiscovery/Miseq104/results
python  __main__.py  pick-leads --folder /Users/Hoan.Nguyen/ComBio/NGS/Projects/AntibodyDiscovery/Miseq104/results
python  __main__.py check-repeats --folder /Users/Hoan.Nguyen/ComBio/NGS/Projects/AntibodyDiscovery/Miseq104/results


# contact {Hoan.Nguyen, Andre.Teixeira}@proteininnovation.org}
