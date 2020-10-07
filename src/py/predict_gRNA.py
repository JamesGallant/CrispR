import os
import subprocess

class Get_gRNA:
	def __init__(self):
		self.root = os.path.dirname(os.path.abspath(__name__))
		self.rscript_path = os.path.join(self.root, "src\\R\\find_gRNA.R")

	def show_class_methods(self):
		methods = ["check_help()", "check_options()",
		"generate_config(out_dir)",
		"create_gRNA_library(config_file, out_dir, gff_file)"]

		for items in methods:
			print(items)


	def check_help(self):
		subprocess.call(["Rscript", "--vanilla", self.rscript_path,
			"--help"], shell=True)

	def check_options(self):
		subprocess.call(["Rscript", "--vanilla", self.rscript_path,
		 "--method", "list_options"], shell=True)

	def generate_config(self, out_dir):
		if out_dir == "":
			assert ValueError("Output directory required")

		subprocess.call(["Rscript", "--vanilla", self.rscript_path,
			"--method", "generate_config",
			"--out_dir", out_dir], shell=True)

	def create_gRNA_library(self, config_file, out_dir):
		subprocess.call(["Rscript", "--vanilla", self.rscript_path,
			"--method", "create_gRNA_library",
			"--config_file", config_file,
			"--out_dir", out_dir], shell=True)