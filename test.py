
from utilities.build_reference import BSgenome
from utilities.install_R_dependencies import Installer

bsgenome = BSgenome()
ins_dep = Installer()

#print(bsgenome._get_version_number(seedfile="D:\\Eigenaara\\Documents\\Python_apps\\CrispR\\references\\dcf_files\\BSgenome_MtbH37Rv-seed.dcf"))
bsgenome.bsgenome_from_seed(seedfile="D:\\Eigenaara\\Documents\\Python_apps\\CrispR\\references\\dcf_files\\test.dcf")

#bsgenome.package_exists(path_to_package="D:\\Eigenaara\\Documents\\Python_apps\\CrispR\\R\\local_packages\\BSgenome.Mtuberculosis.VUmc.H37Rv_20.07.tar.gz")
#bsgenome.list_help()
#bsgenome.list_methods()
#bsgenome.remove_package(seedfile="D:\\Eigenaara\\Documents\\Python_apps\\CrispR\\references\\dcf_files\\BSgenome_MtbH37Rv-seed.dcf")
#print(bsgenome._dcf_metadata(seedfile="D:\\Eigenaara\\Documents\\Python_apps\\CrispR\\references\\dcf_files\\BSgenome_MtbH37Rv-seed.dcf", extract="Package"))
#print(bsgenome.list_avail_methods)

#ins_dep.install()