#'@importFrom data.table rbindlist
#'@importFrom parallel mclapply
#'@importFrom parallel detectCores


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
#' @param parallel if *TRUE*, uses mclapply from parallel package. *default=FALSE*
#' @return A dataframe containing SRR_ID, assigned_read, orig_names, new_names and cellranger_names
#' A csv will be written into the outdir "assigned_SRAreads.csv".
#' @export
assignSRAreads <- function(working_dir=NULL, input_dir=NULL, outdir=NULL, parallel = FALSE) {
  ##set working directory
  setwd(working_dir)

  if(is.null(outdir)){
    system(paste0("rm -r ",working_dir,"/", "read_lengths"))
    system(paste0("mkdir ",working_dir,"/", "read_lengths"))
    outdir <- paste0(working_dir,"/", "read_lengths/")
  }

  if(is.null(working_dir)){
    stop("please provide a working directory")
  }else if(is.null(input_dir)){
    stop("please provide input directory containing all fasterq-dump outputs")
  }else{
    print("First n=250 read lengths will be output into outdir")

    #detect cores
    numCores <- detectCores()-1

    #Extract whitelists from package
    data(tenXv1, tenXv2, tenXv3)
    whitelists <- list(tenXv1, tenXv2, tenXv3)
    names(whitelists) <- c("10xv1", "10xv2", "10xv3")


    ##Get fastqs and create dataframe of references
    fasterq_list <- list.files(input_dir, pattern="fastq.gz$|fastq$", recursive=TRUE,full.names = TRUE)
    orig_names <- gsub(".*SRR", "SRR", fasterq_list)
    orig_names <- gsub(".fastq.*", "", orig_names)
    toParse <- data.frame(fasterq_list=as.character(fasterq_list), orig_names=as.character(orig_names))
    toParse$cat <- ifelse(grepl(".gz", toParse$fasterq_list), "zcat", "cat")

    ##Sample first 250 reads and assign to either R1, R2 or I1 dependent on read lengths
    if (parallel==TRUE){
    assigned_SRA<- mclapply(toParse$fasterq_list, function(x){
      system(paste0(toParse[grepl(x, toParse$fasterq_list),"cat"], " ", x, " | head -1000 | awk '{if(NR%4==2) print length($1)}' > ", outdir,toParse[grepl(x, toParse$fasterq_list),"orig_names"],".readslength.txt"))
      mean_length <- mean(scan(paste0(outdir,toParse[grepl(x, toParse$fasterq_list),"orig_names"],".readslength.txt"), numeric(), quiet = TRUE))

      if(mean_length %in% seq(5,10)){
        assigned_read <- "I3"
        chemistry <- NA
      }else if(mean_length > 10){
        system(paste0(toParse[grepl(x, toParse$fasterq_list),"cat"], " ", x, " | head -40000 | awk '{if(NR%4==2) print /^@/ ? $1 : substr($0,1,16)}' > ", outdir,toParse[grepl(x, toParse$fasterq_list),"orig_names"],".seqs.txt"))
        seqs <- scan(paste0(outdir,toParse[grepl(x, toParse$fasterq_list),"orig_names"],".seqs.txt"),character(), quiet = TRUE)

        whitelist_counts <- cbind(lapply(whitelists, function(x){
          sum(seqs %in% x)
        }))

        if(sum(unlist(whitelist_counts)) > 1000){
          assigned_read <- "R1"
          chemistry <- names(whitelist_counts[which.max(whitelist_counts),])
        }else{
          assigned_read <- "R2"
          chemistry <- NA
        }
      }

      SRR_ID <- toParse[grepl(x, toParse$fasterq_list),"orig_names"]
      assigned <- data.frame(SRR_ID, mean_length,assigned_read, chemistry)
      return(assigned)
    }, mc.cores=numCores)
    }else{
      assigned_SRA<- lapply(toParse$fasterq_list, function(x){
        system(paste0(toParse[grepl(x, toParse$fasterq_list),"cat"], " ", x, " | head -1000 | awk '{if(NR%4==2) print length($1)}' > ", outdir,toParse[grepl(x, toParse$fasterq_list),"orig_names"],".readslength.txt"))
        mean_length <- mean(scan(paste0(outdir,toParse[grepl(x, toParse$fasterq_list),"orig_names"],".readslength.txt"), numeric(), quiet = TRUE))

        if(mean_length %in% seq(5,10)){
          assigned_read <- "I3"
          chemistry <- NA
        }else if(mean_length > 10){
          system(paste0(toParse[grepl(x, toParse$fasterq_list),"cat"], " ", x, " | head -40000 | awk '{if(NR%4==2) print /^@/ ? $1 : substr($0,1,16)}' > ", outdir,toParse[grepl(x, toParse$fasterq_list),"orig_names"],".seqs.txt"))
          seqs <- scan(paste0(outdir,toParse[grepl(x, toParse$fasterq_list),"orig_names"],".seqs.txt"),character(), quiet = TRUE)

          whitelist_counts <- cbind(lapply(whitelists, function(x){
            sum(seqs %in% x)
          }))

          if(sum(unlist(whitelist_counts)) > 1000){
            assigned_read <- "R1"
            chemistry <- names(whitelist_counts[which.max(whitelist_counts),])
          }else{
            assigned_read <- "R2"
            chemistry <- NA
          }
        }

        SRR_ID <- toParse[grepl(x, toParse$fasterq_list),"orig_names"]
        assigned <- data.frame(SRR_ID, mean_length,assigned_read, chemistry)
        return(assigned)
      })
    }

    ##rbind output
    bound <- data.table::rbindlist(assigned_SRA)
    bound$orig_names <- fasterq_list
    bound <- correctNames(bound)
    ##write output to outdir
    write.csv(assigned_SRA, paste0(outdir, "assigned_SRAreads.csv"))
  }
  return(bound)
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
      ##create dummy name to avoid overwriting duplicate names
      assigned_SRA$parent_dir <- gsub("/[^/]*$", "", assigned_SRA$orig_names)
      assigned_SRA$dummy_name <- paste0(assigned_SRA$parent_dir, "/","dummy", gsub(".*SRR","SRR",assigned_SRA$new_name))

      ##assign to dummy name
      for (i in 1:nrow(assigned_SRA)){
        system(paste0( "mv ", assigned_SRA$orig_name[i], " ", assigned_SRA$dummy_name[i]))
      }
      ##assign to new name
      for (i in 1:nrow(assigned_SRA)){
        system(paste0( "mv ", assigned_SRA$dummy_name[i], " ", assigned_SRA$new_name[i]))
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

