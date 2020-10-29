# misc. string values
CHUNK_SIZE, RETRY = 30, 20
# misc strings
NL, SP, NIX, CSV, GDATA = '\n', ' ', '', '.csv', 'Gene_data_'

URL_DBGET = 'https://www.kegg.jp/dbget-bin/www_bget?-f+-n+n+'
URL_KNAP = "'http://www.knapsackfamily.com/knapsack_core/result.php?sname=organism&word="  # knapsack partial URL

# keys for accessing dictionaries
O_KEY, E_KEY, N_KEY, P_KEY, G_KEY = 'ORTHOLOGY', 'EC', 'NTSEQ', 'PLANT', 'GENE'

# regular expressions
RE_KO, RE_EC, RE_ALPH, RE_NT_HEAD, RE_NT_SEQ, RE_SQ_BRACKETS = r'(\[KO:K[^\]]*\])', \
                                                               r'([ |:][0-9]*\.[0-9|-]*\.[0-9|-]*\.[0-9|-]*)', \
                                                               r'[^a-zA-Z]+', r'(&gt;.*)', r'([acgt]{10,})', r"(\[.*\])"
# sections of HTML to remove from knapsack
TO_REMOVE = ['<tr>', '</tr>', '</font>', "\n", "\\n", ">", '</a>', '</td>',
             "target=\"_blank\">", "target=\"_blank\"", "<td class=\"d1\">",
             '<font color=#FF00BF>', '<a href=information.php?word=']

# list of HTTP errors
HTTP_ERRS = ['400', '401', '402', '403', '404', '405', '406', '407', '408', '409', '410', '411',
             '412', '413', '414', '415', '416', '417', '500', '501', '502', '503', '504', '505']

# labels for each of the compounds
AGI, BUN, CAQ, EC, EGT, ERD, GC, GEN, HCC, HWB, KMP, KXN, LU2, MYC, NAR, PYG, QUE, RCO, STL = \
    'Apigenin', 'Butein', 'Catechol', 'Epicatechin', 'Epigallocatechin', 'Eriodictyol', \
    'Gallocatechin', 'Genistein', 'Isoliquiritigenin', 'Cyanidin', 'Kaempferol', 'Catechin', \
    'Luteolin', 'Myricetin', 'Naringenin', 'Pyrogallol', 'Quercetin', 'Resorcinol', 'Resveratrol'
