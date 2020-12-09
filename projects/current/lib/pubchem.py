# Import the library necessary for making a web sevice request.
import urllib.request

from util import *
from pathlib import Path

new_dir = '..' + SEP + 'misc_files' + SEP + 'pubchem' + SEP
results = new_dir + 'results.csv'
Path(new_dir).mkdir(parents=True, exist_ok=True)

url_prefix = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug'
cid_str = '/compound/cid/297,6324,6334,11327430,1234567891'
props = '/property/MolecularFormula,MolecularWeight,XLogP,InChi,CanonicalSMILES,IsomericSMILES'
output = '/CSV'
url = url_prefix + cid_str + props + output

request = urllib.request.urlopen(url)
content = request.read()
with open(results, 'w')as file:
    file.write(content)
