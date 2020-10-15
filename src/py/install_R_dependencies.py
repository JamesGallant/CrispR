import os
import subprocess

class Installer:
    def __init__(self, script):
        self.root = os.path.dirname(os.path.abspath(__name__))
        self.rpath = os.path.join(self.root, "src", "R", f"{script}")
        
    def install_external_bsgenome(self, bsgenome_path):
        subprocess.call(['Rscript', '--vanilla',
                         self.rpath, '--method', 'install',
                         '--filepath', bsgenome_path], shell=False)