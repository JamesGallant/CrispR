if(!require(optparse)){
  install.packages("optparse")
}

source("../R/install_bioconductor.R")

option_list = list(
  make_option("--method", type = "character", default = NULL, action = "store",
              help = "Function to run, use 'list_options' to see methods"),
  make_option("--path", type = "character", default = NULL, action = "store",
              help = "path to file"),
  make_option("--out_dir", type = "character", default = NULL, action = "store",
              help = "output directory"),
  make_option("--package", type = "character", default = NULL, action = "store",
              help = "Name of package to remove")
)

options_parser = OptionParser(option_list = option_list)
command_args = parse_args(options_parser)


if (is.null(command_args$method)) {
  stop("No arguments given, use --help for options")
}

#Functions------------------------------------------------------------------->
avail_methods <- function(){
  #methods for command_args$method
  
  return(c("list_options", "Forge", "Detect", 
           "Build", "Install", "Remove"))
}

detect_package <- function(PATH) {
  if (file.exists(PATH)) {
    print(paste0("file at location:", " ", PATH, " ", "exists"))
  } else {
    return("file does not exist")
  }
}

forge_genome_package <- function(PATH, OUT_DIR) {
  #creates package from seed file
  #requires seed file (.DCM extention)

  if (is.null(PATH)) {
    stop("Requires a seed file (.DCM extention)")
  } else {
    require(GenomicFeatures)
    require(BSgenome)
    
    forgeBSgenomeDataPkg(PATH, destdir = OUT_DIR)
  }
  
}

build_genome_package <- function(PATH, OUT_DIR){
  #Build the package and tars it for installation
  
  if (is.null(PATH)) {
    stop("Requires a folder built with forge_genome_package")
  } else {
    
    if(!require(devtools)){
      install.packages("devtools")
    }
    
    devtools::build(pkg = PATH, path = OUT_DIR)
  }
}

install_genome_package <- function(PATH) {
  #installs the package build with devtools
  if (is.null(PATH)) {
    stop("Requires a tared file built with build_genome_package")
  } else {
    install.packages(PATH, repos = NULL, type = "source")
  }
  
}

remove_genome_package <- function(PACKAGE) {
  if (is.null(PACKAGE)) {
    stop("No package name provider, cant remove")
  } else {
    remove.packages(PACKAGE)
  }
}

#Logic----------------------------------------------------------------------->
if (is.na(match(command_args$method, avail_methods()))) {
  stop("--method arguments don't match availible methods, use list options to see methods")
}

switch (command_args$method,
  "list_options" = {avail_methods()},
  "Detect" = {detect_package(PATH = command_args$path)},
  "Forge" = {forge_genome_package(PATH = command_args$path, OUT_DIR = command_args$out_dir)},
  "Build" = {build_genome_package(PATH = command_args$path, OUT_DIR = command_args$out_dir)},
  "Install" = {install_genome_package(PATH = command_args$path)},
  "Remove" = {remove_genome_package(PACKAGE = command_args$package)}
)