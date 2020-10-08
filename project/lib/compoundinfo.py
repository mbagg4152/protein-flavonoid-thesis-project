ENA = 'EC:4.3.1.24'
ENB = 'EC:4.3.1.25'
ENC = 'EC:1.14.14.91'
END = 'EC:6.2.1.12'
ENE = 'EC:2.3.1.170'
ENF = 'EC:2.3.1.133'
ENG = 'EC:1.14.14.96'
ENH = 'EC:1.14.13.-'
ENI = 'EC:2.3.1.74'
ENJ = 'EC:5.5.1.6'
ENK = 'EC:1.14.20.5'
ENL = 'EC:1.14.19.76'
ENM = 'EC:1.14.14.81'
ENN = 'EC:1.14.14.82'
ENO = 'EC:1.14.11.9'
ENP = 'EC:1.14.20.6'
ENQ = 'EC:1.1.1.219'
ENR = 'EC:1.14.20.4'
ENS = 'EC:1.3.1.77'
ENT = 'EC:1.17.1.3'
ENU = 'EC:1.14.14.87'
ENV = 'EC:4.2.1.105'
ENW = 'EC:2.4.1.74'
ENX = 'EC:2.3.1.70'
ENY = 'EC:2.3.1.30'
ENZ = 'EC:2.4.1.136'

# logical combination of pathways for each chemical
def start(vals): return (ENA or ENB) and ENC and END in vals
def has_nari(vals): return ENI and ENJ in vals and start(vals)
def has_erio(vals): return (ENM or ENN in vals and has_nari(vals)) or (ENH or (ENF and ENG) and ENI in vals)
def has_kaem(vals): return has_nari(vals) and (ENO and ENP in vals)
def has_quer(vals): return (has_erio(vals) and (ENO and ENP in vals)) or (has_kaem(vals) and (ENM or ENN in vals))
def has_cate(vals): return has_erio(vals) and (ENO and ENQ and ENT in vals)
def has_gall(vals): return has_erio and (ENO and ENM and ENQ and ENT in vals)
def has_epig(vals): return has_erio(vals) and (ENO and ENM and ENQ and ENR and ENS in vals)
def has_cyan(vals): return has_erio(vals) and (ENO and ENQ and ENR in vals)
def has_ecat(vals): return has_cyan(vals) and (ENS in vals)
def has_myri(vals): return (has_quer(vals) and (ENM in vals)) or \
                           (has_erio(vals) and (ENM and ENP in vals)) or \
                           (has_nari(vals) and (ENO and (ENM or ENN) and ENM and ENP in vals))
def has_isoflav(vals): return start(vals) and (ENE and ENJ in vals)
def has_apig(vals): return ((ENK or ENL in vals) and has_nari(vals)) or (has_isoflav(vals) and (ENK or ENL in vals))
def has_lute(vals): return has_apig(vals) and (ENM or ENN in vals)
def has_geni(vals): return has_isoflav(vals) and (ENU and ENV in vals)
def has_isol(vals): return has_isoflav(vals) and (ENJ in vals)
def has_bute(vals): return start(vals) and (ENE or ('KO:K08243' and 'KO:00600') in vals)
def has_w(vals): return ENW in vals
def has_x(vals): return ENX in vals
def has_y(vals): return ENY in vals
def has_z(vals): return ENZ in vals

APIG = 'Apigenin'
BUTE = 'Butein'
CATE = 'Catechin'
CYAN = 'Cyanidin'
ECAT = 'Epicatechin'
EPIG = 'Epigallocatechin'
ERIO = 'Eriodictyol'
GALL = 'Gallocatechin'
GENI = 'Genistein'
KAEM = 'Kaempferol'
LUTE = 'Luteolin'
MYRI = 'Myricetin'
NARI = 'Naringenin'
QUER = 'Quercetin'
