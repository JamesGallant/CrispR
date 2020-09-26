if(!require(optparse)){
  install.packages("optparse")
}

option_list = list(
  make_option("--method", type = "character", default = NULL, action = "store",
              help = "Function to run, use 'list_options' to see methods"),
  make_option("--config_file_path", type = "character", default = NULL, action = "store",
              help = "path to configuration file for CrisprSeek settings"), 
  make_option("--out_dir", type = "character", default = NULL, action = "store",
              help = "path to a output directory")

)

options_parser = OptionParser(option_list = option_list)
command_args = parse_args(options_parser)

if (is.null(command_args$method)) {
  stop("No arguments given, use --help for options")
}

#functions------------------------------------------------------------------------>
avail_methods <- function(){
  #methods for command_args$method
  return(c("list_options", "generate_config", "create_gRNA_library"))
}


generate_config <- function(out_file){
  df <- data.frame(organism = "Organism name",
                   circular_chromosome = "logical TRUE/FALSE",
                   input_file = "path to input fasta file",
                   gff_file = "Path to gff file",
                   find_gRNA_with_cutsites = "logical TRUE/FALSE",
                   find_paired_gRNA = "logical TRUE/FALSE",
                   BSgenome = "Name of the BSgenome object (BSgenome.Mtuberculosis.VUmc.H37Rv for H37Rv)",
                   chromosomes_to_search = "all",
                   min_gap = 0,
                   max_gap = 20, 
                   max_mismatch_gRNA = 4,
                   PAM_sequence = "AGAA",
                   gRNA_size = 20,
                   scoring_method = "CFDscore/Hsu-Zhang")
  
  write.table(df, file = paste0(out_file, "/config.txt"), sep = "\t", row.names = F)
}


create_txdb <- function(path_to_gff, organism, circular_chroms, bsgenome_info){
 
  #Path_to_gff can be a url or absolute path
  #orgranism is the name of the organisms
  #Circular_chroms can be true or false
  
  library(GenomicFeatures)
  library(bsgenome_info, character.only = TRUE)
  
  bsgenome_data <- eval(as.name(paste(bsgenome_info)))
  chromosomes <- seqnames(bsgenome_data)
  chrom_len <- sapply(chromosomes, function(x){ return(length(bsgenome_data[[x]]))})
  if (circular_chroms) {
    circ <- sapply(1:length(chromosomes), function(x){return(TRUE)})
  } else {
    circ <- sapply(1:length(chromosomes), function(x){return(FALSE)})
  }
  
  return(makeTxDbFromGFF(file = path_to_gff,
                         format = "auto",
                         dataSource = "User GFF file",
                         organism = organism,
                         chrominfo = data.frame(chrom = chromosomes,
                                                length = chrom_len,
                                                is_circular = circ)))
}


create_gRNA_library <- function(config_file){
  #Requires a configuration file with offTragetAnalysis settings, tab separated
  #takes a single line now, can iterate in a loop once we have more handling
  
  library(CRISPRseek)
  
  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  script.dirname <- dirname(script.name)
  
  
  REpatternFile <- paste0(script.dirname, '\\', "metadata", "\\", "NEBenzymes.fa")
  
  config <- read.delim(file = config_file, sep = "\t", header = TRUE)
  
  txdb_file <- create_txdb(path_to_gff = config[["gff_file"]],
                           organism = config[["organism"]],
                           circular_chroms = config[["circular_chromosome"]],
                           bsgenome_info =  config[["BSgenome"]])
  
  
  offTargetAnalysis(inputFilePath = config[["input_file"]], 
                    findgRNAsWithREcutOnly = config[["find_gRNA_with_cutsites"]], 
                    REpatternFile = REpatternFile, 
                    findPairedgRNAOnly = config[["find_paired_gRNA"]],
                    BSgenomeName = eval(as.name(paste(config[["BSgenome"]]))), 
                    chromToSearch = config[["chromosomes_to_search"]],
                    min.gap = config[["min_gap"]], 
                    max.gap = config[["max_gap"]],
                    gRNA.size = config[["gRNA_size"]],
                    txdb = txdb_file,
                    max.mismatch = config[["max_mismatch_gRNA"]],
                    outputDir = command_args$out_dir, overwrite = T, 
                    PAM = config[["PAM_sequence"]],
                    PAM.size = config[["PAM_length"]],
                    scoring.method = config[["scoring_method"]])
}


#Logic---------------------------------------------------------------------------->
switch (command_args$method,
  'list_options' = {avail_methods()},
  'generate_config' = {generate_config(out_file = command_args$out_dir)},
  'create_gRNA_library' = {create_gRNA_library(config_file = command_args$config_file_path)}
)
