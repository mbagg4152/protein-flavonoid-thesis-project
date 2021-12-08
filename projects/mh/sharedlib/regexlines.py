RE_ALPH = r'[^a-zA-Z]+'
RE_CAS_ID_FORM = r"([0-9]*-[0-9]*-[0-9]*)"  # regular expression for CAS ID
RE_EC = r'([ |:][0-9]*\.[0-9|-]*\.[0-9|-]*\.[0-9|-]*)'


RE_KO = r'(\[KO:K[^\]]*\])'
RE_NT_HEAD = r'(&gt;.*)'
RE_NT_SEQ = r'([acgt]{10,})'
RE_SEARCH_TYPE = r'.*input?\stype?\s=?\s<font?\sclass=\"iw\">\s*(\S*)?\s,?\s</font>'
RE_SQ_BRACKETS = r"(\[.*\])"
