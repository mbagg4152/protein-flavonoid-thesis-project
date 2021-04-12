import os
from parseutil import *

SEP = os.sep

flav_names = 'pajson' + SEP + 'flav_names.json'
flav_rel = 'pajson' + SEP + 'flav_related.json'
flav_syns = 'pajson' + SEP + 'flav_syns.json'
ncomp_codes = 'pajson' + SEP + 'chem_codes.json'
ncomp_dict = 'pajson' + SEP + 'chem_codes_names.json'
nplant_codes = 'pajson' + SEP + 'organism_codes.json'
nplant_dict = 'pajson' + SEP + 'organism_codes_names.json'
plant_names_codes = 'pajson' + SEP + 'plant_names_codes.json'
lig_id_dict = 'pajson' + SEP + 'lig_identifiers.json'
lig_id_part = 'pajson' + SEP + 'lig_identifiers_part.json'
pdb_knap_ids = 'pajson' + SEP + 'pdb_knap_ids.json'

flav_list = get_json_data(flav_names)
flav_relatives = get_json_data(flav_rel)
flav_synonyms = get_json_data(flav_syns)
plant_dict = get_json_data(plant_names_codes)

cmp_codes = '..' + SEP + 'misc_files' + SEP + 'npass' + SEP + 'comp_codes.txt'
comp_dict = '..' + SEP + 'misc_files' + SEP + 'npass' + SEP + 'comp_dict.txt'
comp_in = '..' + SEP + 'misc_files' + SEP + 'npass' + SEP + 'comp_info.txt'
org_codes = '..' + SEP + 'misc_files' + SEP + 'npass' + SEP + 'plant_codes.txt'
org_in = '..' + SEP + 'misc_files' + SEP + 'npass' + SEP + 'npass_species.txt'
org_out = '..' + SEP + 'misc_files' + SEP + 'npass' + SEP + 'parsed_species.txt'
pairs_in = '..' + SEP + 'misc_files' + SEP + 'npass' + SEP + 'plants_compounds.txt'
selected = '..' + SEP + 'misc_files' + SEP + 'npass' + SEP + 'selected.txt'

re_aroma = r'<tr.id=.c[a-z]+Ar.+A.+C.+t.>\s*<th>A.+\sB.+\sC.+t</th>\s*<td>(\d+)</td>\s*</tr>'
re_atom_count = r'\s*<tr id=.ch.*A.*C.*t.>\s*<th>A.*C.*t</th>\s*<td>(.*)</td>\s*</tr>\s*<tr id=.ch.*Ch.*A.*C.*.>\s*<th>'
re_bond = r'<tr.id=.c[a-z]+B.+C.+t.>\s*<th>B.+\sC.+t</th>\s*<td>(\d*)</td>\s*</tr>\s*<tr.id=.ch.+A.+A.+C.+t.>'
re_charge = r'<tr id=.ch.+Fo.+C.+e.>\s*<th>Fo.+\sC.*e</th>\s*<td>([\w\.\-\+]+)</td>'
re_chebi = r">CHEBI:([0-9]*)"
re_chembl = r">CHEMBL([0-9]*)"
re_chiral = r'<tr.id=.c[a-z]+C.+A.+t.>\s*<th>C[a-z]+\sAt.+Co.+t</th>\s*<td>(\d+)</td>\s*</tr>\s*<tr.id=.c.+B\w+C\w+t.>'
re_comp_id = r'(NPC\S*)'
re_comp_name = r'NPC\S*\s*([^\t\r\n]*)'
re_form = r'.*<tr id=.ch.*F.*a.>\s*<th>F.*a</th>\s*<td>(.*)</td>\s*</tr>\s*<tr id=.ch.*Mo.*W.*t\">'
re_has_cas = r'<tr><td class="d1"><a href=information\.php\?word=C[0-9]* target="_blank">(C[0-9]*)<\/a><\/td><td ' \
             r'class="d1">([0-9]*-[0-9]*-[0-9]*)<\/td><td class="d1">([^<]*)<\/td><td class="d1"> '
re_inchi = r"<tr id=\"chemicalInChI\"><th>InChI<\/th><td style=\"word-wrap: break-word\">(InChI=.*)</td></tr><tr " \
           r"id=\"chemicalInChIKey\"> "
re_inchi_key = r"<th>InChIKey<\/th><td>(.*)</td></tr></table></div><div class=\"col-md-4"
re_knap_entry = r'<tr>[.|\s]*<td\sc.*d1\">[.|\s]*<a.*blank\">(C[0-9]*)</a>[.|\s]*</td>[.|\s]*<td.*d1\">([ ' \
                r'0-9|-]*)</td>\s*<td\s.*d1\">(.*)</td>[.|\s]*<td\s.*d1\">([0-9|A-Z]+)</td>[.|\s]*<td\s.*d1\">([' \
                r'\d|\.]+).*(<[|/]td.*<[|/]tr>) '
re_knap_name = r'<th class=\"inf\">Name</th>\s*<td colspan=\"4\" class=\"inf\">(.*)</td>\s*</tr>'
re_knap_org = r'.*</td><td class=\"?org2\"?>([a-zA-Z0-9!@#$&()\\-`.+,/\"\s]*)'
re_long = r'id=.ch[a-z]+I[a-z]+s.>\s*<th>Id.+rs</th>\s*<td.s.+=.word.+word.>(.*)</td></tr><tr.id=.ch[a-z]+Formula'
re_name = r'</ul>\s*</div>\s*<h1.id=.m.+eId.>.{3}</h1>\s*<h4 id=.m.+e.>\s*([^><!/].*)</h4>.*<div.class=.form-group'
re_no_cas = r'<tr><td class="d1"><a h.*n\.php\?word=C[0-9]* target="_blank">(C[0-9]*)<\/a><\/td><td ' \
            r'class="d1"><\/td><td class="d1">([^<]*)<\/td><td class="d1"> '
re_num_knap_results = r'Number of matched data :([0-9]*).*<br>'
re_org_comp_id = r'-(NPC\S*)'
re_org_id = r'(NPO\S*)'
re_org_name = r'NPO[^ \t\r]*\s*(\S*\s*\S*)\s*\bSpecies'
re_org_name_pairs = r'(NPO\S*)\b-'
re_pubchem = r"https:\/\/pubchem\.ncbi\.nlm\.nih\.gov\/compound\/([0-9]*)"
re_smile = r'id=.ch.*Is.*c.>\s*<th>Is.*S.*S</th>\s*<td.*\">(.*)</td>\s*</tr>\s*<tr.id=.ch.*I.*I.>\s*'
re_weight = r'.*<tr id=.ch.*Mo.*W.*t.>\s*<th>Mo.*W.*t</th>\s*<td>(.*)</td>\s*</tr>\s*<tr id=\"ch.*lT.*e\">'

URL_KNAP = "'http://www.knapsackfamily.com/knapsack_core/result.php?sname=organism&word="  # knapsack partial URL
PDB_URL = "http://www.rcsb.org/ligand/"
URL_DBGET = 'https://www.kegg.jp/dbget-bin/www_bget?-f+-n+n+'
K_CHEB = 'CHEBI'
K_CHEM = 'CHEMBL'
K_FORM = 'FORMULA'
K_IKEY = 'INCHI_KEY'
K_INCHI = 'INCHI'
K_ISO = 'ISO_SMILES'
K_LONG = 'LONG'
K_NAME = 'COMMON'
K_PUB = 'PUBCHEM'
K_SMILE = 'SMILES'
S_IKEY = 'INCHIKEY'
S_INCHI = 'INCHI_CD'
S_METAB = 'metabolite'
