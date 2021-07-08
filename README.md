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

To avoid unwanted behaviours, always input absolute paths. I have used relative paths in this example for brevity. 
If outdir is left empty, the function *assignSRAreads()* will make a "read_lengths" directory in the working_dir. 

The functions *assignSRAreads()* and *renameAll()* will always output a csv file named "assigned_SRAreads.csv" into outdir. This will contain the updated dataframes which will index original SRA run IDs, fastq output names and the proposed new names. This way one may change the names back easily in place. 

```
##Assign directories
working_dir = "~/Project_X/"
input_dir = "~/Project_X/fasterq_output/"
outdir="/well/combat/users/vkh192/public/CSF_data/read_lengths/"

##Assign SRA reads. This will output a datatable containing columns: "SRR_ID", "assigned_read", "new_names", "cellranger_names"
assigned_files <- assignSRAreads(working_dir = working_directory, input_dir = input_dir, outdir =outdir)

##Correct names. Here format can be assigned as "read_correc" or "cellranger". read_correct simple corrects the _1/_2/_3 suffixes to the correct assignments. cellranger, will rename all files to cellranger compatible formats SRRXXXX_S1_L001_RX_001.fastq

renameAll(assigned_SRA= assigned_files, input_dir=input_dir, format="cellranger")
```

I hope you all enjoy using this. When I have time, I will write functionality for autodetecting 10X chemistries as I have noticed these are often inaccurately assigned in SRA and papers. Whilst this is fine for cellranger which auto-detects against whitelists, for those using kallisto bustools, it can be very frustrating. 

**Happy public dataset hunting!**
