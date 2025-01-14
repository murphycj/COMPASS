The script `preprocess.py` can be used to process a loom file (generated by the MissionBio pipeline) into the input expected by COMPASS.

## Dependencies
The dependencies required to run the script can be installed by running:

`pip install numpy pandas loompy varcode`

`pyensembl install --release 75 76`

## Usage

`python preprocess.py -i [input] -o [output_dir]`

The input can either be a single loom file, or a directory containing several loom files.
Optional arguments:

| Argument      | Description |
| ----------- | ----------- |
| --whitelist   | Path to csv file containing the mutations to include (see `mutations.csv` as an example)        |
| --SNP     | Path to a tsv file containing the population frequency of variants (see below)       |
| --region   | Region to use for inferring CNVs; must be either 'gene' or 'amplicon' (default: gene).       |
| --panel   | Path to a csv file describing the amplicons. This is sometimes useful to get the correct name of the amplicons.  |

## Whitelist

In case, you already know which mutations to include in the analysis, the list of mutations can be passed as the whitelist argument (see `mutations.csv` as an example). Otherwise, the script will attempt to select somatic variants present in some cells but not others, as well as germline variants that might be affected by LOH.

## Population Frequency

Optionally, the preprocessing script can take as input a tsv file containing the population frequency of variants. This is used in the preprocessing to remove germline variants (unless they appear to be affected by LOH in some cells) and, in COMPASS, to penalize variants with a high population frequency which are not placed at the root (since they are likely to be germline variants).
The file that we used was generated using the script `download_1000G.sh`, which was adapted from [this script](https://github.com/single-cell-genetics/cellSNP/blob/master/SNPlist_1Kgenome.sh).
