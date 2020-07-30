import unittest
import os

from utilities.build_reference import BSgenome

bsgenome = BSgenome()

def build_bsgenome():
	bsgenome.bsgenome_from_seed(seedfile="D:\\Eigenaara\\Documents\\Python_apps\\CrispR\\tests\\testing_resources\\test.dcf")

def remove_test_package():
	bsgenome.remove_package(seedfile="D:\\Eigenaara\\Documents\\Python_apps\\CrispR\\tests\\testing_resources\\test.dcf")

def file_creation():
	root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
	build_bsgenome()
	check_file = os.listdir(os.path.join(root, "R/local_packages"))

	if "BSgenome.testorg.test_0.0.tar.gz" in check_file:
		out = True
	else:
		out = False

	remove_test_package()
	os.remove(os.path.join(root, "R/local_packages/BSgenome.testorg.test_0.0.tar.gz"))

	return out


class Bsgenome_test(unittest.TestCase):

	def test_file_creation(self):
		self.assertEqual(True, file_creation())


if __name__ == '__main__':
    unittest.main()