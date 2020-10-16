# misc. string values
CHUNK_SIZE = 30
RETRY = 20
CSV = '.csv'
GDATA = 'Gene_data_'
KS_URL = "'http://www.knapsackfamily.com/knapsack_core/result.php?sname=organism&word="  # knapsack partial URL
NIX = ""
NL = '\n'
SP = ' '
RE_SQ_BRACKETS = "(\[.*\])"
O_KEY = 'ORTHOLOGY'
E_KEY = 'EC'
N_KEY = 'NTSEQ'
P_KEY = 'PLANT'

to_remove = ['<tr>',
             '</tr>',
             '</font>',
             "\n",
             "\\n",
             ">",
             '</a>',
             '</td>',
             "target=\"_blank\">",
             "target=\"_blank\"",
             "<td class=\"d1\">",
             '<font color=#FF00BF>',
             '<a href=information.php?word='
             ]
http_errs = ['400',
             '401',
             '402',
             '403',
             '404',
             '405',
             '406',
             '407',
             '408',
             '409',
             '410',
             '411',
             '412',
             '413',
             '414',
             '415',
             '416',
             '417',
             '500',
             '501',
             '502',
             '503',
             '504',
             '505'
             ]

AGI = 'Apigenin'
BUN = 'Butein'
CAQ = 'Catechol'
EC = 'Epicatechin'
EGT = 'Epigallocatechin'
ERD = 'Eriodictyol'
GC = 'Gallocatechin'
GEN = 'Genistein'
HCC = 'Isoliquiritigenin'
HWB = 'Cyanidin'
KMP = 'Kaempferol'
KXN = 'Catechin'
LU2 = 'Luteolin'
MYC = 'Myricetin'
NAR = 'Naringenin'
PYG = 'Pyrogallol'
QUE = 'Quercetin'
RCO = 'Resorcinol'
STL = 'Resveratrol'
