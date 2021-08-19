# fasterqParseR

**Declaration**: *This tool was designed for the sole purpose of parsing fasterq/fastq-dump outputs for single cell RNAseq experiments, in particular, 10X genomics datasets.*

This tool is a quick and dirty method for checking filenames of fasterq/fastq-dump outputs. It requires the outputs obtained with the following option arguments:

fasterq-dump --include-technical --split-files</br>
or</br>
fastq-dump --split-files</br>

The tool will then take an input directory and iteratively sample the top 250 reads and assign them to R1, R2 or I1 using the following paramaters and generate correct names based on the following parameters:

- If the mean lengths of the sampled reads are <=10bp, they will be assigned as "I"
- If the reads are over 10bp, they will be automatically matched against 10X whitelists for v1, v2 and v3 chemistries. If >1000 barcodes are detected, they will be assigned to "R1" 
- When matching against whitelists, 10X chemistries are also returned as the chemistry with the highest number of matches to the respective whitelist. 
- Lastly, if <1000 barcodes are detected it will be assigned to "R2"

The motivation for writing this was born out of repeated misordering of fasterq/fastq-dump outputs based on the order of upload of the depositing researchers. I hope it makes some of your lives that little bit easier! All feedback is welcome and if there are any unstable or unwanted behaviours, please do report these in the issues section. 

# Instructions for use

```
library(devtools)
install_github('Nusob888/fasterqParseR')
```

To avoid unwanted behaviours, always input absolute paths. I have used relative paths in this example for brevity. 
If outdir is left empty, the function *assignSRAreads()* will make a "read_lengths" directory in the working_dir. 

The functions *assignSRAreads()* and *renameAll()* will always output a csv file named "assigned_SRAreads.csv" into outdir. This will contain the updated dataframes which will index original SRA run IDs, fastq output names and the proposed new names. This way one may change the names back easily in place. 

**Important considerations** 
1) Always check the new assignments. In some cases, SRA deposits contain duplicated fastqs representing the same reads despite being labelled as _1 or _2. Sadly, this is something that needs to be taken up with SRA and the authors. 

2) fasterqParseR currently does not support scATAC. Unfortunately there is no way of identifying R1 and R3 as both are 50bp. 

3) parallelisation is supported via the parallel package. This is currently turned off as default to avoid accidently hogging resources on clusters. It can be activated by: @param parallel = TRUE

```
##Assign directories
working_dir = "~/Project_X/"
input_dir = "~/Project_X/fasterq_output/"
outdir="~/Project_X/read_lengths/"

##Assign SRA reads. This will output a datatable containing columns: "SRR_ID", "assigned_read", "new_names", "cellranger_names"
assigned_files <- assignSRAreads(working_dir = working_directory, input_dir = input_dir, outdir =outdir, parallel=TRUE)

##Correct names. Here format can be assigned as "read_correct" or "cellranger". read_correct simple corrects the _1/_2/_3 suffixes to the correct assignments. cellranger, will rename all files to cellranger compatible formats SRRXXXX_S1_L001_RX_001.fastq

renameAll(assigned_SRA= assigned_files, input_dir=input_dir, format="cellranger")
```

I hope you all enjoy using this small package. 

**Happy public dataset hunting!**
