from flib.fconstants import *

def flav_check(label, ec_list):
    # determine which function to call based on the label passed in using the label function pairs
    options = {
        AGI: agi, BUN: bun, KXN: kxn, HWB: hwb, EC: ec, EGT: egt, ERD: erd, GC: gc, GEN: gen, KMP: kmp, LU2: lu2,
        MYC: myc, NAR: nar, QUE: que, GGT: ggt, DEC: dec, SOA: soa, V1G: v1g, HCC: bun, DHKM: dhkm, DHMY: dhmy,
        DHQU: dhqu, DLM: dlm, LDLM: ldlm, LHWB: lhwb, NACH: nach
    }
    try:
        res = options[label](ec_list)
    except KeyError:
        res = False
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

def wca(e): return (nca(e) and (E03 in e)) or (hc4(e) and (E04 in e))  # p-coumaroyl-coa

def nach(e): return wca(e) and (E09 in e)  # naringenin chalcone

def nar(e): return nach(e) and (E10 in e)  # naringenin

def agi(e): return nar(e) and or_in(e, E11, E12)  # apigenin

def lu2(e): return agi(e) and (or_in(e, E13, E14))  # luteolin

def fca(e): return wca(e) and (and_in(e, E06, E07) or (E08 in e))  # caffeoyl-coa

def erch(e): return (nar(e) and or_in(e, E13, E14)) or (fca(e) and (E09 in e))  # eriodictyol chalcone

def erd(e): return erch(e)  # eriodictyol. same as erch for now, one of the direct steps on map

def dhkm(e): return nar(e) and (E15 in e)  # dihydrokaempferol

def kmp(e): return dhkm(e) and (E16 in e)  # kaempferol

def dhqu(e): return (dhkm(e) and or_in(e, E13, E14)) or (erd(e) and (E15 in e))  # dihydroquercetin

def lhwb(e): return dhqu(e) and or_in(e, E17, E17_2)  # leucocyanidin

def que(e): return (kmp(e) and or_in(e, E13, E14)) or (dhqu(e) and (E16 in e))  # quercetin

def kxn(e): return lhwb(e) and (E20 in e)  # catechin

def hwb(e): return lhwb(e) and (E18 in e)  # cyanidin

def ec(e): return hwb(e) and (E19 in e)  # epicatechin

def dhmy(e): return (dhqu(e) and (E13 in e)) or (erd(e) and and_in(e, E13, E15))  # dihydromyricetin

def myc(e): return (dhmy(e) and (E16 in e)) or (que(e) and (E13 in e))  # myricetin

def ldlm(e): return dhmy(e) and or_in(e, E17, E17_2)  # leucodelphinidin

def dlm(e): return ldlm(e) and (E18 in e)  # delphinidin

def gc(e): return ldlm(e) and (E20 in e)  # gallocatechin

def egt(e): return dlm(e) and (E19 in e)  # epigallocatechin

def gen(e): return nar(e) and and_in(e, E21, E22)  # genistein

def bun(e): return wca(e) and (or_in(e, E05, E09))  # butein

def ggt(e): return GGT in e  # Glycosaminoglycan galactosyltransferase

def dec(e): return DEC in e  # Deleted entry

def soa(e): return SOA in e  # Serine O-acetyltransferase

def v1g(e): return V1G in e  # vanillate 1-glucosyltransferase
