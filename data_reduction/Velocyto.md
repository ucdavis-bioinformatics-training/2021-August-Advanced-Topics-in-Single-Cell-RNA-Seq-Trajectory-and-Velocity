
Login to tadpole and navigate to your directory on the share space.

```bash
srun -t 1-00:00:00 -c 16 -n 1 --mem-per-cpu 2000 --partition production --account workshop --reservation workshop  --pty /bin/bash
```

We're going to need access to our home directories from the interactive session, so we need to run the aklog command.

```bash
aklog
```

#  RNA Velocity measurement using Velocyto

RNA velocity is the time derivative of the gene expression state,
 [(La Manno et al., 2018)](https://www.nature.com/articles/s41586-018-0414-6) allows for the inference of the dynamic patterns in scRNA-seq data sets, by looking at the abundance of unspliced and spliced mRNA RNA in each cell, and modelling using a system of ordinary differential equations.

From a quantification point of view, RNA velocity analysis requires the generation of two count matrices, representing the spliced/processed (exonic) and unspliced/unprocessed RNA (intronic).

RNA velocity is a high-dimensional vector that predicts the future state of a cell on a timescale of hours.


### Example dataset

The data we're using here is unpublished human immune cell data that has been subsetted and manipulated. We'll be working with the output of cellranger multi. Specifically, velocyto.py uses the barcodes.tsv.gz and sample_alignments.bam files, along with the reference GTF used by cellranger.

Let's set up the project directory

```bash
mkdir -p /share/workshop/adv_scrnaseq/$USER
cd /share/workshop/adv_scrnaseq/$USER
ln -s /share/biocore/workshops/2021_08_Trajectory_Velocity/01-Cellranger .
ln -s /share/biocore/workshops/2021_08_Trajectory_Velocity/references .
```

### Next lets install the software [velocyto](https://velocyto.org/)

To install velocyto (a python application) we are going to use conda and a virtual environment.

```bash
cd /share/workshop/adv_scrnaseq/$USER
module load anaconda3
conda create -p /share/workshop/adv_scrnaseq/$USER/velocyto
conda activate /share/workshop/adv_scrnaseq/$USER/velocyto
```

If the environment 'activated' properly, than your prompt should look something like this:

```
(/share/workshop/adv_scrnaseq/hslyman/velocyto) hslyman@tadpole:/share/workshop/adv_scrnaseq/hslyman$
```

If not, you may have gotten an error message:
```
CommandNotFoundError: Your shell has not been properly configured to use 'conda activate'.
To initialize your shell, run

    $ conda init <SHELL_NAME>

Currently supported shells are:
  - bash
  - fish
  - tcsh
  - xonsh
  - zsh
  - powershell

See 'conda init --help' for more information and options.

IMPORTANT: You may need to close and restart your shell after running 'conda init'.
```

In order to initialize our shell while in an interactive session, we need to run a few extra lines of code.

```bash
aklog
conda init bash
source ~/.bashrc
conda activate /share/workshop/adv_scrnaseq/$USER/velocyto
```

Once your conda environment is properly activated you can then install the software.
```bash
# install prerequisites
conda install numpy scipy cython numba matplotlib scikit-learn h5py click
# install velocyto
pip install velocyto
```

To test the installation, print Velocyto's help statement.
```bash
velocyto --help
```

It is also recommended to use a repetative sequence mask file, this is easiest to obtain from the [UCSC genome browser](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=611454127_NtvlaW6xBSIRYJEBI0iRDEWisITa&clade=mammal&org=Human&db=0&hgta_group=allTracks&hgta_track=rmsk&hgta_table=rmsk&hgta_regionType=genome&position=&hgta_outputType=gff&hgta_outFileName=GRCh38_rmsk.gtf
).

You want to download the mask in GTF format. Unfortunately, Ensembl doesn't provide a GTF file for the repeat sequences and we'd have to generate one ourselves (write a script to do so) which is beyond the scope of this workshop.

###  Running Velocyto on Cellranger output

Velocyto has a helpful run10x function, which is a wrapper around the run function with some preset parameters that allow you to get away with typing less on the command line. However, the output folders produced by cellranger multi differ slightly from those produced by cellranger count, so we will be using the more general "run" function. Let's first take a look at the help documentation (also available [online](https://velocyto.org/velocyto.py/tutorial/cli.html#run-run-on-any-technique-advanced-use)).

```bash
velocyto run --help
```
We're ready to run velocyto on our first sample!

```bash
cd /share/workshop/adv_scrnaseq/$USER/
velocyto run --bcfile /share/workshop/adv_scrnaseq/$USER/01-Cellranger/Sample1/outs/multi/count/raw_feature_bc_matrix/barcodes.tsv.gz \
             --mask /share/workshop/adv_scrnaseq/$USER/references/GRCh38_rmsk.gtf \
             --outputfolder /share/workshop/adv_scrnaseq/$USER/02-Velocyto \
             --samtools-threads 24 \
             --samtools-memory 2000 \
             /share/workshop/adv_scrnaseq/$USER/01-Cellranger/Sample1/outs/per_sample_outs/Sample1/count/sample_alignments.bam \
             /share/workshop/adv_scrnaseq/$USER/references/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
             2> velocyto_sample1.err > velocyto_sample1.out
```
Even with only 10,000,000 reads, this takes a while. You can observe the progress of your velocyto job by viewing the tail of velocyto_sample1.out. We can progress to the R analysis with output from all three samples made earlier this week.

### Output

Velocyto produces a single [loom](https://linnarssonlab.org/loompy/format/index.html) file containing the needed matrices for the analysis.

To download all three output loom files and the cellranger count data, run the following from the terminal on your local computer, making sure to replace "user" with your username.

```bash
scp user@tadpole.genomecenter.ucdavis.edu:/share/biocore/workshops/2021_08_Trajectory_Velocity/data_download.tar.gz .
tar -xzv data_download.tar.gz
```
We will be using these files in R, so move them to your R project directory.

### More reading

A detailed study of the impact of quantification on RNA velocity estimates and interpretation, see [Soneson et al., 2020](https://www.biorxiv.org/content/10.1101/2020.03.13.990069v1).
