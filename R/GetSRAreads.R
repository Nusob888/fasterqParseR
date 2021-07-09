#'@importFrom data.table rbindlist


##assign reads to likely R1, R2, I3
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

##correct names
correctNames <- function(input = assigned_SRA){
  if(is.null(input)){
    stop("Please supply the output of GetSRAreads")
  }else if(any(grepl("fastq.gz", input$orig_names))){
    new_names <- gsub("SRR.*", "", input$orig_names)
    input$new_names <- paste0(new_names, gsub("\\_.*","",input$SRR_ID), gsub("R|I", "_",input$assigned_read), ".fastq.gz")
    input$cellranger_names <- paste0(new_names, gsub("\\_.*","",input$SRR_ID),"_S1_L001_" ,input$assigned_read, "_001.fastq.gz")
  }else if(any(grepl("fastq", input$orig_names))){
    new_names <- gsub("SRR.*", "", input$orig_names)
    input$new_names <- paste0(new_names, gsub("\\_.*","",input$SRR_ID), gsub("R|I", "_",input$assigned_read), ".fastq")
    input$cellranger_names <- paste0(new_names, gsub("\\_.*","",input$SRR_ID),"_S1_L001_" ,input$assigned_read, "_001.fastq")
  }else{
    stop("Oops... It appears you have no assigned fastq paths")
  }
  return(input)
}


#' Assign SRA reads for desposited 10X single cell RNAseq datasets
#'
#' This function takes an input directory containing fasterq-dump or fastq-dump files.
#' Fasterq and fastq outputs must be generated from using the split option.
#' Example:
#' fasterq-dump --include-technical --split-files
#' or
#' fastq-dump --split-files
#'**PLEASE NOTE** This will not take the output of fasterq/fastq-dump if you specified the --type argument
#'
#' 250 reads from all fastq files are then sampled from the listed fastqs and assigned to R1, R2 or I3.
#' 10,000 reads will also be sampled if the reads are matched for R1
#'
#' @param working_dir Set working directory. Provide absolute paths
#' @param input_dir Set to directory containing fasterq-dump or fastq-dump files
#' @param outdir Directory to output the read legnths as SRRXXXXX_X.readlength.txt
#' @param get_chemistry If *TRUE* will take first 10,000 reads and match to 10X whitelists stored within
#' the package.  *default = FALSE*
#' @return A dataframe containing SRR_ID, assigned_read, orig_names, new_names and cellranger_names
#' A csv will be written into the outdir "assigned_SRAreads.csv".
#' @export
assignSRAreads <- function(working_dir=NULL, input_dir=NULL, outdir=NULL, get_chemistry=FALSE) {

  if(is.null(working_dir)){
    stop("please provide a working directory")
  }else if(is.null(input_dir)){
    stop("please provide input directory containing all fasterq-dump outputs")
  }else if(is.null(outdir)){
    system(paste0("mkdir ", "read_lengths"))
  }else{
    print("First n=250 read lengths will be output into outdir")
    ##set working directory
    setwd(working_dir)

    if(isTRUE(get_chemistry)){
    #Extract whitelists from package
    data(tenXv1, tenXv2, tenXv3)
    whitelists <- list(tenXv1, tenXv2, tenXv3)
    names(whitelists) <- c("10xv1", "10xv2", "10xv3")
    }

    ##Get fastqs and create dataframe of references
    fasterq_list <- list.files(input_dir, pattern="fastq.gz$|fastq$", recursive=TRUE,full.names = TRUE)
    orig_names <- gsub(".*SRR", "SRR", fasterq_list)
    orig_names <- gsub(".fastq.*", "", orig_names)
    toParse <- data.frame(fasterq_list=as.character(fasterq_list), orig_names=as.character(orig_names))

    ##Sample first 250 reads and assign to either R1, R2 or I1 dependenton read lengths
    assigned_SRA<- lapply(toParse$fasterq_list, function(x){
      system(paste0("zcat ", x, " | head -1000 | awk '{if(NR%4==2) print length($1)}' > ", outdir,toParse[x,"orig_names"],".readslength.txt"))
      mean_length <- mean(scan(paste0(outdir,toParse[x,"orig_names"],".readslength.txt"), numeric(), quiet = TRUE))
     assigned_read<- assignRead(mean_length)

    if(isTRUE(get_chemistry)){
      if(assigned_read == "R1"){
        system(paste0("zcat ", x, " | head -40000 | awk '{if(NR%4==2) print /^@/ ? $1 : substr($0,1,16)}' > ", outdir,toParse[x,"fastq_names"],".seqs.txt"))
        seqs <- scan(paste0(outdir,toParse[x,"fastq_names"],".seqs.txt"),character(), quiet = TRUE)

        #Sum matched seqs to whitelists and take version of greatest matches
        whitelist_counts <- cbind(lapply(whitelists, function(x){
          sum(seqs %in% x)
        }))
        version <- names(whitelist_counts[which.max(whitelist_counts),])
      }else{
        version <- NA
      }

      SRR_ID <- toParse[x,"orig_names"]
      assigned <- data.frame(SRR_ID, assigned_read, version)
      return(assigned)

    }else{

      SRR_ID <- toParse[x,"orig_names"]
      assigned <- data.frame(SRR_ID, assigned_read)
      return(assigned)

    }
    })

    ##rbind output
    bound <- data.table::rbindlist(assigned_SRA)
    bound$orig_names <- fasterq_list
    assigned_SRA <- correctNames(bound)

    ##write output to outdir
    write.csv(assigned_SRA, paste0(outdir, "assigned_SRAreads.csv"))
  }
  return(assigned_SRA)
}


#' Rename all fastas to the correct filenames
#'
#' This function takes the output of *assignSRAreads* and directly changes the filename in place.
#'
#' @param assigned_SRA Input the output of *assignSRAreads*
#' @param input_dir Input the same input_dir argument as used for *assignSRAreads*
#' @return The "assigned_SRAreads.csv" file will be updated with a name_check column containing
#' the boolean outcome of listed fastq names matched to the intended new_names column.
#' @export
renameAll <- function(assigned_SRA = NULL, input_dir=NULL, format = NULL){
  if(is.null(assigned_SRA)){
  stop("Please input the output of assignSRAreads")
  }else if(is.null(input_dir)){
    stop("Please fill input_dir with same input file used for assignSRAreads")
  } else if(any(!grepl("read_correct|cellranger", format))){
    stop("Please decide on a renaming format from 'read_correct' or 'cellranger'")
  }else if(format == "read_correct"){
    for (i in 1:nrow(assigned_SRA)){
      system(paste0( "mv ", assigned_SRA$orig_name[i], " ", assigned_SRA$new_name[i]))
    }
  }else if(format == "cellranger"){
    for (i in 1:nrow(assigned_SRA)){
      system(paste0( "mv ", assigned_SRA$orig_name[i], " ", assigned_SRA$cellranger_names[i]))
    }
  }
  name_check <- list.files(input_dir, pattern="fastq.gz$|fastq$", recursive=TRUE,full.names = TRUE)
  assigned_SRA$name_check <- name_check %in% assigned_SRA$new_name
  ##write output to outdir
  write.csv(assigned_SRA, paste0(outdir, "assigned_SRAreads.csv"))
  print("All done! Please check in assigned_SRAreads.csv for new column name_check and make sure you are happy!")
}


