# tcr-germline-paper
### JH @ MGH, 2025

The contents of this repo produce the figures in the [currently in submission] manuscript from Heather *et al.*, in which we use various publicly available data resources to illustrate why the field of TCR research perhaps needs to reconsider how they handle germline references.

## Directions for use

The following instructions expect Python scripts to be run in a UNIX-like environment. They were tested with Python 3.12.0 running on Mac OSX.

### 1) Clone this repo
  
```bash
git clone https://github.com/JamieHeather/tcr-germline-paper.git
cd tcr-germline-paper
```

### 2) Get the GitHub data

* This script uses the contents of two other repositories in which I've collected various TCR-related datasets:
  * https://github.com/JamieHeather/genedb-releases (in which I've banked a number of IMGT/GENE-DB releases since 2013)
  * https://github.com/JamieHeather/novel-tcr-alleles/ (where I've collated a number of reported novel TCR alleles from the literature)

```bash
git clone https://github.com/JamieHeather/genedb-releases.git
git clone https://github.com/JamieHeather/novel-tcr-alleles.git
```

### 3) Install the necessary dependencies

The Python scripts require several packages which are not included in the Standard Library, which can be installed with `pip`:

```bash
pip install numpy pandas pyalluvial matplotlib seaborn
```

It also requires the installation of William Lees' [receptor_utils](https://williamdlees.github.io/receptor_utils/_build/html/introduction.html) package, which can also be installed via `pip` (`pip install receptor_utils`). The package is actually run as in the terminal by the scripts through the `subprocess` module, so the installed package needs to be available in the PATH during runtime. 

### 4) Get and pre-process the TCRseq data

* The data from the [Mikelov *et al.* 2024 Genome Research paper](doi.org/10.1101/gr.278775.123) needs to be downloaded to the `mikelov-data/` directory 

  * This is accessible from [EBI with the accession 'E-MTAB-13593'](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13593)
  * There are many ways to access this data, e.g. via FTP
  * Whichever route used, users need to ensure they end up with the right amount of data in the same format 
    * There should be two files for each sample, read one and read two, named `ERR12340xxx_[1/2].fastq` (which might be gzipped)

* Downloaded TCRseq reads need to be merged
  * Given the read length and the nature of the TCR sequencing protocol, this allows productive of complete variable domain sequences where available
  * I achieved this with [FLASH](https://doi.org/10.1093/bioinformatics/btr507), e.g. using the following loop, but other tools should behave comparably:

```bash
for r1 in ERR*1.fastq.gz
   do
   nam=$(echo $r1 | cut -d '_' -f 1)
   r2=$(echo $r1 | sed 's/_1/_2/')
   echo $nam $r2
   flash $r1 $r2 -o flash-$nam -O -M 350
done
```

* Merged reads then need to have their TCRs annotated
  * This is achieved using [`autoDCR` version 0.2.7](https://github.com/JamieHeather/autoDCR), and the `collapseTCRs.py` script located in the `mikelov-data` directory 
  * This will identify TCR rearrangements across all merged read, and then save a reduced format in which each V/J/CDR3 occurs only one, respectively
  * This can similarly be executed across all relevant merged read files:

```bash
for merge in flash*.fastq.gz
	do 
		nam=$(echo $merge | cut -d '-' -f 2 | cut -d '.' -f 1)
		echo $nam 
		python autoDCR.py -fq $merge -o dcr-$nam
	    python quick-collapse.py -in dcr-$nam
	done
```

### 5) Run the analysis scripts

In order to run the analysis, users need to navigate to the scripts directory and run the plotting scripts in turn. Note that the use of relative paths means they have to be run from this folder:

```bash
cd scripts
python plot_genedb_alleles.py 
python plot_gene_usage.py
python plot_allele_confidence.py
python plot_germline_referencing.py
```

If executing properly, these will each generate date-prefixed output subdirectories in the `plots/` folder, containing the individual plots (and occasional intermediate tables) used in the manuscript.

Note that the TCR analysis section also deposits some of its intermediate tables in the `reference-data/` directory, so that the plotting analysis can be repeated without having to re-process all of the individual repertoire files. A version of those files has been included in this repo so that the plotting code can be re-run without the need to download and process the large dataset: if users wish to repeat the entire analysis these files need to be removed before running the `plot_gene_usage.py` script.

