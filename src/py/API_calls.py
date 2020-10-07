import requests

class TaxonomyAPI:
    def __init__(self, genus, species):

        self.url = f"https://www.ebi.ac.uk/ena/taxonomy/rest/any-name/{genus}%20{species}"

    def convert(self, lst):
        res_dct = {lst[i]: lst[i + 1] for i in range(0, len(lst), 2)}
        return res_dct

    def get_taxon(self):
        call = requests.get(self.url).json()
        call_single = call[0]
        return repr(call_single)