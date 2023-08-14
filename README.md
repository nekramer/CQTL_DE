# CQTL_DE

The pipeline used for QC'ing and processing RNA-seq data for differential expression analysis
in the CQTL project.

## Workflow
1. Clone this entire repo into working directory:

    ```bash
    git clone https://github.com/nekramer/CQTL_DE.git
    ```

2. Edit the comma-separated `samplesheet.csv` to include sequencing samples of interest. The helper script
`makeSamplesheet.R` was used for pulling various subsets of datasets from our internal sequencing
samplesheets.

3. Edit `config/config.yaml` for parameters specific to analysis.

4. Submit workflow with `sbatch`:

    ```bash
    sbatch runRNAproc
    ```

To relaunch this workflow, unlock the directory with:

    ```bash
    module load python/3.6.6
    ./unlock.sh RNAproc #or
    ./unlock.sh runRNAproc
    ```