import os 
import shutil
import subprocess


class BSgenome:
	def __init__(self):
		self.root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
		self.genome_build_funcs = os.path.join(self.root, "R/scripts/BSgenome_builder.R")

		if "temp" not in os.listdir(self.root):
			os.mkdir(os.path.join(self.root, "temp"))

		self.tempdir = os.path.join(self.root, "temp")
		self.package_dir = os.path.join(self.root, "R/local_packages")


	def _get_version_number(self, seedfile):
		with open(seedfile, 'r') as dcf_file:
			for line in dcf_file:
				if "Version:" in line:
					return line.split()[1]
				

	def _package_exists(self):
		pass

	def list_help(self):
		pass

	def list_methods(self):
		""""
		Lists the availible options in the Rscript
		"""
		subprocess.call(['Rscript','--vanilla', self.genome_build_funcs, '--method', 'list_options'], shell=True)


	def bsgenome_from_seed(self, seedfile, package_name):
		"""
		requires a seed file (SEED). 
		see https://www.bioconductor.org/packages//2.7/bioc/vignettes/BSgenome/inst/doc/BSgenomeForge.pdf on how to make one
		The seed will point to the fasta file
		"""
		subprocess.call(['Rscript', 
			'--vanilla', 
			self.genome_build_funcs, 
			'--method', 'Forge',
			'--path', seedfile,
			'--out_dir', self.tempdir], shell=True)

		subprocess.call(['Rscript', 
			'--vanilla', 
			self.genome_build_funcs, 
			'--method', 'Build',
			'--path', os.path.join(self.root, "temp", package_name),
			'--out_dir', self.package_dir], shell=True)

		subprocess.call(['Rscript', 
			'--vanilla', 
			self.genome_build_funcs, 
			'--method', 'Install',
			'--path', os.path.join(self.package_dir, f"{package_name}_{self._get_version_number(seedfile)}.tar.gz")], shell=True)

		shutil.rmtree(os.path.join(self.root, "temp"), ignore_errors=True)








	

