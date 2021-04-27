from flib.fconstants import *

def flav_check(label, ec_list):
    # determine which function to call based on the label passed in using the label function pairs
    # options = {AGI: agi, APIF: apif, AZEL: azel, BUN: bun, BUTN: butn, CHSN: chsn, DAID: daid, DEC: dec, DF74: df74,
    #            DFV: dfv, DGEN: dgen, DKAM: dkam, DTRI: dtri, DMYR: dmyr, DQUE: dque, DLM: dlm, EC: ec, EGT: egt,
    #            ERDC: erdc, ERD: erd, EZEL: ezel, FSTN: fstn, G50: g50, GBAN: gban, GC: gc, GEN: gen, GGT: ggt,
    #            GALN: galn, HCC: hcc, HDAI: hdai, HESP: hesp, HWB: hwb, KMP: kmp, KXN: kxn, LDEL: ldel, LCYN: lcyn,
    #            LPEL: lpel, LU2: lu2, LUTF: lutf, MYC: myc, MYF: myf, NARC: narc, NAR: nar, PBAN: pban, PICH: pich,
    #            PINO: pino, PELR: pelr, QUE: que, SOA: soa, TET2: tet2, TISO: tiso, TNON: tnon, V1G: v1g}
    try: res = label(ec_list)
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

def pich(e): return nca(e) and (E09 in e)  # pinocembrin chalcone

def pino(e): return pich(e) and (E10 in e)  # pinocembrin

def chsn(e): return pino(e) and (E11 in e)  # chrysin

def pban(e): return pino(e) and (E15 in e)  # pinobanksin

def galn(e): return pban(e) and (E16 in e)  # galangin

def wca(e): return (nca(e) and (E03 in e)) or (hc4(e) and (E04 in e))  # p-coumaroyl-coa

def g50(e): return wca(e) and (and_in(e, E24, E09))  # phloretin

def narc(e): return wca(e) and (E09 in e)  # naringenin chalcone

def nar(e): return narc(e) and (E10 in e)  # naringenin

def dgen(e): return nar(e) and (E21 in e)  # 2-Hydroxy-2,3-dihydrogenistein

def hesp(e): return nar(e)  # hesperetin

def apif(e): return nar(e) and (or_in(e, E17_1, E17_FULL, E17_2))  # apiforol

def agi(e): return nar(e) and or_in(e, E11, E12)  # apigenin

def lu2(e): return agi(e) and (or_in(e, E13, E14))  # luteolin

def fca(e): return wca(e) and (and_in(e, E06, E07) or (E08 in e))  # caffeoyl-coa

def erdc(e): return (nar(e) and or_in(e, E13, E14)) or (fca(e) and (E09 in e))  # eriodictyol chalcone

def erd(e): return erdc(e)  # eriodictyol. same as erch for now, one of the direct steps on map

def lutf(e): return erd(e) and (or_in(e, E17_1, E17_FULL, E17_2))  # luteoforol

def dtri(e): return erd(e) and (E13 in e)  # dihydrotricetin

def myf(e): return lutf(e) and (E13 in e) or (dtri(e) and (E11 in e))  # tricetin

def dkam(e): return nar(e) and (E15 in e)  # dihydrokaempferol

def lpel(e): return dkam(e) and (E17_1 in e)  # leucopelargonidin

def pelr(e): return lpel(e) and (E18 in e)  # pelargonidin

def ezel(e): return pelr(e) and (E19 in e)  # epiafzelechin

def azel(e): return lpel(e) and (E20 in e)  # afzelechin

def kmp(e): return dkam(e) and (E16 in e)  # kaempferol

def dque(e): return (dkam(e) and or_in(e, E13, E14)) or (erd(e) and (E15 in e))  # dihydroquercetin

def lcyn(e): return dque(e) and or_in(e, E17_1, E17_FULL)  # leucocyanidin

def que(e): return (kmp(e) and or_in(e, E13, E14)) or (dque(e) and (E16 in e))  # quercetin

def kxn(e): return lcyn(e) and (E20 in e)  # catechin

def hwb(e): return lcyn(e) and (E18 in e)  # cyanidin

def ec(e): return hwb(e) and (E19 in e)  # epicatechin

def dmyr(e): return (dque(e) and (E13 in e)) or (erd(e) and and_in(e, E13, E15))  # dihydromyricetin

def myc(e): return (dmyr(e) and (E16 in e)) or (que(e) and (E13 in e))  # myricetin

def ldel(e): return dmyr(e) and or_in(e, E17_1, E17_FULL)  # leucodelphinidin

def dlm(e): return ldel(e) and (E18 in e)  # delphinidin

def gc(e): return ldel(e) and (E20 in e)  # gallocatechin

def egt(e): return dlm(e) and (E19 in e)  # epigallocatechin

def gen(e): return nar(e) and and_in(e, E21, E22)  # genistein

def bun(e): return wca(e) and (or_in(e, E05, E09))  # butein

def butn(e): return bun(e) and (E10 in e)  # butin

def hcc(e): return bun(e)  # isoliquiritigenin

def dfv(e): return hcc(e) and (E10 in e)  # liquiritigenin

def gban(e): return dfv(e) and (E15 in e)  # garbanzol

def fstn(e): return (gban(e) and (E14 in e)) or (butn(e) and (E15 in e))  # dihydrofisetin

def df74(e): return dfv(e) and (or_in(e, E11, E12))  # 7,4'-dihydroxyflavone

def tiso(e): return dfv(e) and (E21 in e)  # 2,7,4'-Trihydroxyisoflavanone

def daid(e): return tiso(e) and (E22 in e)  # daidzein

def tnon(e): return dfv(e) and (E08 in e)  # 6,7,4'-Trihydroxyflavanone

def tet2(e): return tnon(e) and (E21 in e)  # 2,6,7,4'-Tetrahydroxyisoflavanone

def hdai(e): return tet2(e)  # 6-Hydroxydaidzein

def ggt(e): return E_GGT in e  # Glycosaminoglycan galactosyltransferase

def dec(e): return E_DEC in e  # Deleted entry

def soa(e): return E_SOA in e  # Serine O-acetyltransferase

def v1g(e): return E_V1G in e  # vanillate 1-glucosyltransferase
