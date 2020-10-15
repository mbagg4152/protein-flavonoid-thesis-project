from lib.miscvals import *

E1 = 'EC:4.3.1.24'
E2 = 'EC:4.3.1.25'
E3 = 'EC:1.14.14.91'
E4 = 'EC:6.2.1.12'
E5 = 'EC:2.3.1.170'
E6 = 'EC:2.3.1.133'
E7 = 'EC:1.14.14.96'
E8 = 'EC:1.14.13.-'
E9 = 'EC:2.3.1.74'
E10 = 'EC:5.5.1.6'
E11 = 'EC:1.14.20.5'
E12 = 'EC:1.14.19.76'
E13 = 'EC:1.14.14.81'
E14 = 'EC:1.14.14.82'
E15 = 'EC:1.14.11.9'
E16 = 'EC:1.14.20.6'
E17 = 'EC:1.1.1.219'
E17_2 = 'EC:1.1.1.219 1.1.1.234'
E18 = 'EC:1.14.20.4'
E19 = 'EC:1.3.1.77'
E20 = 'EC:1.17.1.3'
E21 = 'EC:1.14.14.87'
E22 = 'EC:4.2.1.105'
E23 = 'EC:2.4.1.74'
E24 = 'EC:2.3.1.70'
E25 = 'EC:2.3.1.30'
E26 = 'EC:2.4.1.136'

K1 = 'KO:K08243'
K2 = 'KO:00600'
C1 = 'C00079'  # phenylalanine

def flav_check(label, ec_list):
    options = {
        AGI: agi, BU: bu, CA: ca, CYN: cyn, EC: ecat, EGCG: egcg, ERD: erd, GC: gc, GEN: gen, KMP: kmp, LU2: lu2,
        MYC: myri, NAR: nar, QUE: quer, E23: ec23, E24: ec24, E25: ec25, E26: ec26, HCC: bu
    }
    res = ''
    try: res = options[label](ec_list)
    except KeyError: res = False
    return res

# returns true if at least 1 arg is in the list
def or_in(items, *args):
    found = 0
    for a in args:
        if a in items: found += 1
    if found > 0: return True
    else: return False

# returns true only if all args are in the list
def and_in(items, *args):
    for a in args:
        if a not in items: return False
    return True

def cia(e): return or_in(e, E1, E2)  # cinnamic acid
def pca(e): return cia(e) and (E3 in e)  # p-coumaric acid
def cicoa(e): return cia(e) and (E4 in e)  # cinnamoyl-coa
def pcoa(e): return (cicoa(e) and (E3 in e)) or (pca(e) and (E4 in e))  # p-coumaroyl-coa
def narc(e): return pcoa(e) and (E9 in e)  # naringenin chalcone
def nar(e): return narc and (E10 in e)  # naringenin
def agi(e): return nar(e) and or_in(e, E11, E12)  # apigenin
def lu2(e): return agi(e) and (or_in(e, E13, E14))  # luteolin
def cacoa(e): return and_in(e, E6, E7) or (E8 in e)  # caffeoyl-coa
def erd(e): return (nar(e) and or_in(e, E13, E14)) or (cacoa(e) and (E9 in e))  # eriodictyol
def dhk(e): return nar(e) and (E15 in e)  # dihydrokaempferol
def kmp(e): return dhk(e) and (E16 in e)  # kaempferol
def dhq(e): return (dhk(e) and or_in(e, E13, E14)) or (erd(e) and (E15 in e))  # dihydroquercetin
def lc(e): return dhq(e) and or_in(e, E17, E17_2)  # leucocyanidin
def quer(e): return (kmp(e) and or_in(e, E13, E14)) or (dhq(e) and (E16 in e))  # quercetin
def ca(e): return lc(e) and (E20 in e)  # catechin
def cyn(e): return lc(e) and (E18 in e)  # cyanidin
def ecat(e): return cyn(e) and (E19 in e)  # epicatechin
def dhm(e): return (dhq(e) and (E13 in e)) or (erd(e) and and_in(e, E13, E15))  # dihydromyricetin
def myri(e): return (dhm(e) and (E16 in e)) or (quer(e) and (E13 in e))  # myricetin
def ldn(e): return dhm(e) and or_in(e, E17, E17_2)  # leucodelphinidin
def dpn(e): return ldn(e) and (E18 in e)  # delphinidin
def gc(e): return ldn(e) and (E20 in e)  # gallocatechin
def egcg(e): return dpn(e) and (E19 in e)  # epigallocatechin
def gen(e): return nar(e) and and_in(e, E21, E22)  # genistein
def bu(e): return pcoa(e) and (or_in(e, E5, E9))  # butein
def ec23(e): return E23 in e
def ec24(e): return E24 in e
def ec25(e): return E25 in e
def ec26(e): return E26 in e
