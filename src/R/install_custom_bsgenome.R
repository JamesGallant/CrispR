# Title     : Install user bsgenome
# Objective : installer
# Created by: James Gallant
# Created on: 2020/10/12
if(!require(optparse)){
  install.packages("optparse")
}

option_list = list(
  make_option("--method", type = "character", default = NULL, action = "store",
              help = "Method"),

 make_option("--filepath", type = "character", default = NULL, action = "store",
              help = "Filepath to bsgenome package"))

options_parser = OptionParser(option_list = option_list)
command_args = parse_args(options_parser)

if (is.null(command_args$method)) {
  stop("No arguments given, use --help for options")
}

installer <- function (filepath) {
  install.packages(filepath, repos = NULL, type = "source")
}
print(command_args$filepath)
switch(command_args$method,
       'install': installer(filepath = command_args$filepath))

