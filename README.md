# fasterqParseR

**Declaration**: *This tool was designed for the sole purpose of parsing fasterq/fastq-dump outputs for single cell RNAseq experiments, in particular, 10X genomics datasets.*

This tool is a quick and dirty method for checking filenames of fasterq/fastq-dump outputs. It requires the outputs obtained with the following option arguments:

fasterq-dump --include-technical --split-files</br>
or</br>
fastq-dump --split-files</br>

The tool will then take an input directory and iteratively sample the top 250 reads and assign them to R1, R2 or I1 using the following paramaters and generate correct names based on these assignments or cellranger compatible names based on these assignments. 

```
assignRead <- function(mean_length=mean_length){
  if(mean_length %in% seq(20, 30)){
    assign <- "R1"
  }else if(mean_length %in% seq(80, 100)){
    assign <- "R2"
  }else if(mean_length %in% seq(5,10)){
    assign <- "I3"
  }
  return(assign)
}
```
The motivation for writing this was born out of repeated misordering of fasterq/fastq-dump outputs based on the order of upload of the depositing researchers. I hope it makes some of your lives that little bit easier! All feedback is welcome and if there are any unstable or unwanted behaviours, please do report these in the issues section. 

# Instructions for use

```
library(devtools)
install_github('Nusob888/fasterqParseR')
```

To avoid unwanted behaviours, always input absolute paths. I have used relative paths in this example for brevity. 
If outdir is left empty, the function *assignSRAreads()* will make a "read_lengths" directory in the working_dir. 

The functions *assignSRAreads()* and *renameAll()* will always output a csv file named "assigned_SRAreads.csv" into outdir. This will contain the updated dataframes which will index original SRA run IDs, fastq output names and the proposed new names. This way one may change the names back easily in place. 

```
##Assign directories
working_dir = "~/Project_X/"
input_dir = "~/Project_X/fasterq_output/"
outdir="~/Project_X/read_lengths/"

##Assign SRA reads. This will output a datatable containing columns: "SRR_ID", "assigned_read", "new_names", "cellranger_names"
assigned_files <- assignSRAreads(working_dir = working_directory, input_dir = input_dir, outdir =outdir)

##Correct names. Here format can be assigned as "read_correc" or "cellranger". read_correct simple corrects the _1/_2/_3 suffixes to the correct assignments. cellranger, will rename all files to cellranger compatible formats SRRXXXX_S1_L001_RX_001.fastq

renameAll(assigned_SRA= assigned_files, input_dir=input_dir, format="cellranger")
```

## Get correct chemistry versions
Another issue with SRA deposits and indeed papers, is when 10X chemistries are incorrectly stated, or when chemistries are mixed without metadata to index. This is not an issue for cellranger, however does become an issue for other tools requiring the user to state an input chemistry. 

My personal advice for kallisto users, is to always recover chemistries this way and not rely on papers and SRA meta to be correct. 

To obtain correct chemistries, use option *get_chemistry* within the *assignSRAreads()* function. This is defaulted to *FALSE*. </br>
If *TRUE* the first 16bp of the top 10,000 reads from the corrected R1 assignments will be extracted and matched against 10X genomics [whitelists](https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist-). 



```
##If chemistries need to be cross checked, set get_chemistry to TRUE
assigned_files <- assignSRAreads(working_dir = working_directory, input_dir = input_dir, outdir =outdir, get_chemistry=TRUE)

versions <- assigned_files[!is.na(assigned_files$chemistry),c("SRR_ID", "chemistry")]

```

I hope you all enjoy using this small package. 

**Happy public dataset hunting!**
