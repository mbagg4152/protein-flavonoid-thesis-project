from flib.fconstants import *

def flav_check(label, ec_list):
    # determine which function to call based on the label passed in using the label function pairs
    options = {AGI: agi, APIF: apif, AZEL: azel, BUN: bun, BUTN: butn, CHSN: chsn, DDZN: ddzn, DEC: dec, DF74: df74,
               DFV: dfv, DHGN: dhgn, DHKM: dhkm, DHMF: dhmf, DHMY: dhmy, DHQU: dhqu, DLM: dlm, EC: ec, EGT: egt,
               ERCH: erch, ERD: erd, EZEL: ezel, FSTN: fstn, G50: g50, GBZL: gbzl, GC: gc, GEN: gen, GGT: ggt,
               GLGN: glgn, HCC: hcc, HDZ6: hdz6, HESP: hesp, HWB: hwb, KMP: kmp, KXN: kxn, LDLM: ldlm, LHWB: lhwb,
               LPLR: lplr, LU2: lu2, LUTF: lutf, MYC: myc, MYF: myf, NACH: nach, NAR: nar, PBAN: pban, PCCH: pcch,
               PCEM: pcem, PLRG: plrg, QUE: que, SOA: soa, T2674: t2674, T274: t274, T674: t674, V1G: v1g}
    try: res = options[label](ec_list)
    except KeyError: res = False
    return res

# returns true if at least 1 arg is in the list
def or_in(items, *args):
    found = 0
    for a in args:
        if a in items: found += 1
    if found > 0:
        return True
    else:
        return False

# returns true only if all args are in the list
def and_in(items, *args):
    for a in args:
        if a not in items: return False
    return True

def tca(e): return or_in(e, E01, E02)  # cinnamic acid

def hc4(e): return tca(e) and (E03 in e)  # p-coumaric acid

def nca(e): return tca(e) and (E04 in e)  # cinnamoyl-coa

def pcch(e): return nca(e) and (E09 in e)  # pinocembrin chalcone

def pcem(e): return pcch(e) and (E10 in e)  # pinocembrin

def chsn(e): return pcem(e) and (E11 in e)  # chrysin

def pban(e): return pcem(e) and (E15 in e)  # pinobanksin

def glgn(e): return pban(e) and (E16 in e)  # galangin

def wca(e): return (nca(e) and (E03 in e)) or (hc4(e) and (E04 in e))  # p-coumaroyl-coa

def g50(e): return wca(e) and (and_in(e, E24, E09))  # phloretin

def nach(e): return wca(e) and (E09 in e)  # naringenin chalcone

def nar(e): return nach(e) and (E10 in e)  # naringenin

def dhgn(e): return nar(e) and (E21 in e)  # 2-Hydroxy-2,3-dihydrogenistein

def hesp(e): return nar(e)  # hesperetin

def apif(e): return nar(e) and (or_in(e, E17_1, E17_FULL, E17_2))  # apiforol

def agi(e): return nar(e) and or_in(e, E11, E12)  # apigenin

def lu2(e): return agi(e) and (or_in(e, E13, E14))  # luteolin

def fca(e): return wca(e) and (and_in(e, E06, E07) or (E08 in e))  # caffeoyl-coa

def erch(e): return (nar(e) and or_in(e, E13, E14)) or (fca(e) and (E09 in e))  # eriodictyol chalcone

def erd(e): return erch(e)  # eriodictyol. same as erch for now, one of the direct steps on map

def lutf(e): return erd(e) and (or_in(e, E17_1, E17_FULL, E17_2))  # luteoforol

def dhmf(e): return erd(e) and (E13 in e)  # dihydrotricetin

def myf(e): return lutf(e) and (E13 in e) or (dhmf(e) and (E11 in e))  # tricetin

def dhkm(e): return nar(e) and (E15 in e)  # dihydrokaempferol

def lplr(e): return dhkm(e) and (E17_1 in e)  # leucopelargonidin

def plrg(e): return lplr(e) and (E18 in e)  # pelargonidin

def ezel(e): return plrg(e) and (E19 in e)  # epiafzelechin

def azel(e): return lplr(e) and (E20 in e)  # afzelechin

def kmp(e): return dhkm(e) and (E16 in e)  # kaempferol

def dhqu(e): return (dhkm(e) and or_in(e, E13, E14)) or (erd(e) and (E15 in e))  # dihydroquercetin

def lhwb(e): return dhqu(e) and or_in(e, E17_1, E17_FULL)  # leucocyanidin

def que(e): return (kmp(e) and or_in(e, E13, E14)) or (dhqu(e) and (E16 in e))  # quercetin

def kxn(e): return lhwb(e) and (E20 in e)  # catechin

def hwb(e): return lhwb(e) and (E18 in e)  # cyanidin

def ec(e): return hwb(e) and (E19 in e)  # epicatechin

def dhmy(e): return (dhqu(e) and (E13 in e)) or (erd(e) and and_in(e, E13, E15))  # dihydromyricetin

def myc(e): return (dhmy(e) and (E16 in e)) or (que(e) and (E13 in e))  # myricetin

def ldlm(e): return dhmy(e) and or_in(e, E17_1, E17_FULL)  # leucodelphinidin

def dlm(e): return ldlm(e) and (E18 in e)  # delphinidin

def gc(e): return ldlm(e) and (E20 in e)  # gallocatechin

def egt(e): return dlm(e) and (E19 in e)  # epigallocatechin

def gen(e): return nar(e) and and_in(e, E21, E22)  # genistein

def bun(e): return wca(e) and (or_in(e, E05, E09))  # butein

def butn(e): return bun(e) and (E10 in e)  # butin

def hcc(e): return bun(e)  # isoliquiritigenin

def dfv(e): return hcc(e) and (E10 in e)  # liquiritigenin

def gbzl(e): return dfv(e) and (E15 in e)  # garbanzol

def fstn(e): return (gbzl(e) and (E14 in e)) or (butn(e) and (E15 in e))  # dihydrofisetin

def df74(e): return dfv(e) and (or_in(e, E11, E12))  # 7,4'-dihydroxyflavone

def t274(e): return dfv(e) and (E21 in e)  # 2,7,4'-Trihydroxyisoflavanone

def ddzn(e): return t274(e) and (E22 in e)  # daidzein

def t674(e): return dfv(e) and (E08 in e)  # 6,7,4'-Trihydroxyflavanone

def t2674(e): return t674(e) and (E21 in e)  # 2,6,7,4'-Tetrahydroxyisoflavanone

def hdz6(e): return t2674(e)  # 6-Hydroxydaidzein

def ggt(e): return GGT in e  # Glycosaminoglycan galactosyltransferase

def dec(e): return DEC in e  # Deleted entry

def soa(e): return SOA in e  # Serine O-acetyltransferase

def v1g(e): return V1G in e  # vanillate 1-glucosyltransferase
