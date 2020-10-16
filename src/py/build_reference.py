import os
import shutil
import subprocess

class BSgenome:
    def __init__(self):
        self.root = os.path.dirname(os.path.abspath(__name__))
        self.genome_build_funcs = os.path.join(self.root, "src/r/BSgenome_builder.R")

        if "temp" not in os.listdir(self.root):
            os.mkdir(os.path.join(self.root, "temp"))

        self.tempdir = os.path.join(self.root, "temp")
        self.package_dir = os.path.join(self.root, "src/local_packages")

    def _dcf_metadata(self, seedfile, extract):
        """
		Possible error handling needed: Is the extract valid?
		Requires a dcf file (seed file, full path) and a keyword to extract
		Returns string
		"""

        if os.path.splitext(seedfile)[1] != ".dcf":
            raise ValueError("dcf extention not detected, is this a dcf file?")

        with open(seedfile, 'r') as dcf_file:
            for line in dcf_file:
                if extract in line:
                    return line.split()[1]

    def package_exists(self, path_to_package):
        subprocess.call(['Rscript', '--vanilla', self.genome_build_funcs, '--method', 'Detect',
                         '--path', path_to_package], shell=True)

        shutil.rmtree(os.path.join(self.root, "temp"), ignore_errors=True)

    def list_help(self):
        """
		List th help in R
		"""
        subprocess.call(['Rscript', '--vanilla', self.genome_build_funcs, '--help'], shell=True)
        shutil.rmtree(os.path.join(self.root, "temp"), ignore_errors=True)

    def list_methods(self):
        """
		Lists the availible options in the Rscript
		"""
        subprocess.call(['Rscript', '--vanilla', self.genome_build_funcs, '--method', 'list_options'], shell=True)
        shutil.rmtree(os.path.join(self.root, "temp"), ignore_errors=True)

    def remove_package(self, package):
        """
		Requires the package name only, handle paths internally
		"""
        subprocess.call(['Rscript',
                         '--vanilla',
                         self.genome_build_funcs,
                         '--method', 'Remove',
                         '--package', package])

        shutil.rmtree(os.path.join(self.root, "temp"), ignore_errors=True)

    def bsgenome_from_seed(self, seedfile):
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
                         '--out_dir', self.tempdir])

        subprocess.call(['Rscript',
                         '--vanilla',
                         self.genome_build_funcs,
                         '--method', 'Build',
                         '--path',
                         os.path.join(self.root, "temp", self._dcf_metadata(seedfile=seedfile, extract='Package')),
                         '--out_dir', self.package_dir])

        subprocess.call(['Rscript',
                         '--vanilla',
                         self.genome_build_funcs,
                         '--method', 'Install',
                         '--path', os.path.join(self.package_dir,
                                                f"{self._dcf_metadata(seedfile=seedfile, extract='Package')}_{self._dcf_metadata(seedfile=seedfile, extract='Version')}.tar.gz")])

        shutil.rmtree(os.path.join(self.root, "temp"), ignore_errors=True)
