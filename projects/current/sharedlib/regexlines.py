RE_ALPH = r'[^a-zA-Z]+'
RE_AROMA = r'<tr.id=.c[a-z]+Ar.+A.+C.+t.>\s*<th>A.+\sB.+\sC.+t</th>\s*<td>(\d+)</td>\s*</tr>'
RE_ATOM_COUNT = r'\s*<tr id=.ch.*A.*C.*t.>\s*<th>A.*C.*t</th>\s*<td>(.*)</td>\s*</tr>\s*<tr id=.ch.*Ch.*A.*C.*.>\s*<th>'
RE_BOND = r'<tr.id=.c[a-z]+B.+C.+t.>\s*<th>B.+\sC.+t</th>\s*<td>(\d*)</td>\s*</tr>\s*<tr.id=.ch.+A.+A.+C.+t.>'
RE_CAS_ID_FORM = r"([0-9]*-[0-9]*-[0-9]*)"  # regular expression for CAS ID
RE_CHARGE = r'<tr id=.ch.+Fo.+C.+e.>\s*<th>Fo.+\sC.*e</th>\s*<td>([\w\.\-\+]+)</td>'
RE_CHEBI = r">CHEBI:([0-9]*)"
RE_CHEMBL = r">CHEMBL([0-9]*)"
RE_CHIRAL = r'<tr.id=.c[a-z]+C.+A.+t.>\s*<th>C[a-z]+\sAt.+Co.+t</th>\s*<td>(\d+)</td>\s*</tr>\s*<tr.id=.c.+B\w+C\w+t.>'
RE_EC = r'([ |:][0-9]*\.[0-9|-]*\.[0-9|-]*\.[0-9|-]*)'
RE_FORM = r'.*<tr id=.ch.*F.*a.>\s*<th>F.*a</th>\s*<td>(.*)</td>\s*</tr>\s*<tr id=.ch.*Mo.*W.*t\">'
RE_INCHI = r"<tr id=\"chemicalInChI\"><th>InChI<\/th><td style=\"word-wrap: break-word\">(InChI=.*)</td></tr><tr " \
           r"id=\"chemicalInChIKey\"> "
RE_INCHI_KEY = r"<th>InChIKey<\/th><td>(.*)</td></tr></table></div><div class=\"col-md-4"
RE_KNAP_ENTRY = r'<tr>[.|\s]*<td\sc.*d1\">[.|\s]*<a.*blank\">(C[0-9]*)</a>[.|\s]*</td>[.|\s]*<td.*d1\">([ ' \
                r'0-9|-]*)</td>\s*<td\s.*d1\">(.*)</td>[.|\s]*<td\s.*d1\">([0-9|A-Z]+)</td>[.|\s]*<td\s.*d1\">([' \
                r'\d|\.]+).*(<[|/]td.*<[|/]tr>) '
RE_KNAP_NAME = r'<th class=\"inf\">Name</th>\s*<td colspan=\"4\" class=\"inf\">(.*)</td>\s*</tr>'
RE_KNAP_ORG = r'.*</td><td class=\"?org2\"?>([a-zA-Z0-9!@#$&()\\-`.+,/\"\s]*)'
RE_KO = r'(\[KO:K[^\]]*\])'
RE_LONG = r'id=.ch[a-z]+I[a-z]+s.>\s*<th>Id.+rs</th>\s*<td.s.+=.word.+word.>(.*)</td></tr><tr.id=.ch[a-z]+Formula'
RE_NAME = r'</ul>\s*</div>\s*<h1.id=.m.+eId.>.{3}</h1>\s*<h4 id=.m.+e.>\s*([^><!/].*)</h4>.*<div.class=.form-group'
RE_NT_HEAD = r'(&gt;.*)'
RE_NT_SEQ = r'([acgt]{10,})'
RE_NUM_KNAP_RESULTS = r'Number of matched data :([0-9]*).*<br>'
RE_PUBCHEM = r"https:\/\/pubchem\.ncbi\.nlm\.nih\.gov\/compound\/([0-9]*)"
RE_SEARCH_TYPE = r'.*input?\stype?\s=?\s<font?\sclass=\"iw\">\s*(\S*)?\s,?\s</font>'
RE_SMILE = r'id=.ch.*Is.*c.>\s*<th>Is.*S.*S</th>\s*<td.*\">(.*)</td>\s*</tr>\s*<tr.id=.ch.*I.*I.>\s*'
RE_SQ_BRACKETS = r"(\[.*\])"
RE_WEIGHT = r'.*<tr id=.ch.*Mo.*W.*t.>\s*<th>Mo.*W.*t</th>\s*<td>(.*)</td>\s*</tr>\s*<tr id=\"ch.*lT.*e\">'
