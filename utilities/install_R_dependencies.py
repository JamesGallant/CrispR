import os
import subprocess

class Installer:
	def __init__(self):
		self.root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
		self.r_path = os.path.join(self.root, "R/scripts/install_bioconductor.R")

	def install(self):
		subprocess.call(['Rscript', '--vanilla',
			self.r_path], shell=True)
