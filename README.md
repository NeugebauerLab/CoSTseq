# CoSTseq (Co-transcriptional structure tracking)
Data processing pipeline and analysis code for handling CoSTseq and DMS-MaPseq sequencing data sets.

## Installing the CoSTseq package

CoSTseq uses Snakemake and a custom conda environment. In addition, the following software packages need to be installed and accessible for Snakemake: `fastp`, `STAR`, `UMICollapse`, `RNAstructure`, `samtools`. Set the path to each software package in the `config/snake_CoST.yaml` file.

```bash
# clone repository
git clone https://github.com/NeugebauerLab/CoSTseq.git
cd CoSTseq

# create conda environment
conda env create --name CoST --file=CoST_env.yml # create new environment from template
conda activate CoST # activate

# install python package
pip install -e .
```

## Running the data processing pipeline

After completing the configuration process, the pipeline can be executed using the following command:

```bash
conda activate CoST
snakemake -c 16 -d smk_rundir/test --configfile config/snake_CoST.yaml --resources mem_mb=32000 --rerun-incomplete --use-conda
```

## Re-generating analyses and figures from "Rapid folding of nascent RNA regulates eukaryotic RNA biogenesis"

To re-generate the analyses, first run the analysis pipeline on deposited raw data, or download processed files from GEO accession number [GSE254264](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254264). Save the processed data files in the `data` directory. Analysis code for each individual figure is available in jupyter notebooks in the `notebooks` directory and can be executed directly.
