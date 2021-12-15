# custom project library imports
from flib.fconstants import *

# take the label, which is the same as the function name & then call the function.
def flav_check(label, ec_list):
   try:
      res = label(ec_list)
   except KeyError:
      res = False
   return res

# returns true if at least 1 arg is in the list
def or_in(items, *args):
   for a in args:
      if a in items:
         return True
   return False

# returns true only if all args are in the list
def and_in(items, *args):
   for a in args:
      if a not in items:
         return False  # all values must be present
   return True

# The prediction functions. There is one for numerous compounds not being predicted
# as they are precursors to the compounds of interest, which makes the predictions
# a little more readable. Function names are either the 3 letter code assigned
# to a specific compound by the Protein Data Bank (PDB) or are alternative codes
# for compounds that either do not have a PDB code or do not have a variable-friendly
# PDB code (for example, the PDB code for Epicatechin cannot be used because it is
# 28E and the python interpreter does not allow functions to start with a number)
# Function order is determined by a compounds requirements as functions cannot
# be called before they are written.

def tca(e):  # cinnamic acid || precursor
   return or_in(e, E01, E02)

def hc4(e):  # p-coumaric acid || precursor
   return tca(e) and (E03 in e)

def nca(e):  # cinnamoyl-coa || precursor
   return tca(e) and (E04 in e)

def pich(e):  # pinocembrin chalcone || new
   return nca(e) and (E09 in e)

def pino(e):  # pinocembrin || new
   return pich(e) and (E10 in e)

def chsn(e):  # chrysin || new
   return pino(e) and (E11 in e)

def pban(e):  # pinobanksin || new
   return pino(e) and (E15 in e)

def galn(e):  # galangin || new
   return pban(e) and (E16 in e)

def wca(e):  # p-coumaroyl-coa || precursor
   return (nca(e) and (E03 in e)) or (hc4(e) and (E04 in e))

def g50(e):  # phloretin || new
   return wca(e) and (and_in(e, E24, E09))

def narc(e):  # naringenin chalcone || new
   return wca(e) and (E09 in e)

def nar(e):  # naringenin || mh thesis
   return narc(e) and (E10 in e)

def dgen(e):  # 2-Hydroxy-2,3-dihydrogenistein || new
   return nar(e) and (E21 in e)

def hesp(e):  # hesperetin || new
   return nar(e)

def apif(e):  # apiforol || new
   return nar(e) and (or_in(e, E17_1, E17_FULL, E17_2))

def agi(e):  # apigenin || mh thesis
   return nar(e) and or_in(e, E11, E12)

def lu2(e):  # luteolin || mh thesis
   return agi(e) and (or_in(e, E13, E14))

def fca(e):  # caffeoyl-coa || precursor
   return wca(e) and (and_in(e, E06, E07) or (E08 in e))

def erdc(e):  # eriodictyol chalcone || new
   return (nar(e) and or_in(e, E13, E14)) or (fca(e) and (E09 in e))

def erd(e):  # eriodictyol  || mh thesis || same as erch for now
   return erdc(e)

def lutf(e):  # luteoforol || new
   return erd(e) and (or_in(e, E17_1, E17_FULL, E17_2))

def dtri(e):  # dihydrotricetin || new
   return erd(e) and (E13 in e)

def myf(e):  # tricetin || new
   return lutf(e) and (E13 in e) or (dtri(e) and (E11 in e))

def dkam(e):  # dihydrokaempferol || new
   return nar(e) and (E15 in e)

def lpel(e):  # leucopelargonidin || new
   return dkam(e) and (E17_1 in e)

def pelr(e):  # pelargonidin || new
   return lpel(e) and (E18 in e)

def ezel(e):  # epiafzelechin || new
   return pelr(e) and (E19 in e)

def azel(e):  # afzelechin || new
   return lpel(e) and (E20 in e)

def kmp(e):  # kaempferol ||  mh thesis
   return dkam(e) and (E16 in e)

def dque(e):  # dihydroquercetin || new
   return (dkam(e) and or_in(e, E13, E14)) or (erd(e) and (E15 in e))

def lcyn(e):  # leucocyanidin || new
   return dque(e) and or_in(e, E17_1, E17_FULL)

def que(e):  # quercetin || mh thesis
   return (kmp(e) and or_in(e, E13, E14)) or (dque(e) and (E16 in e))

def kxn(e):  # catechin || mh thesis
   return lcyn(e) and (E20 in e)

def hwb(e):  # cyanidin || mh thesis
   return lcyn(e) and (E18 in e)

def ec(e):  # epicatechin || mh thesis
   return hwb(e) and (E19 in e)

def dmyr(e):  # dihydromyricetin || new
   return (dque(e) and (E13 in e)) or (erd(e) and and_in(e, E13, E15))

def myc(e):  # myricetin ||  mh thesis
   return (dmyr(e) and (E16 in e)) or (que(e) and (E13 in e))

def ldel(e):  # leucodelphinidin || new
   return dmyr(e) and or_in(e, E17_1, E17_FULL)

def dlm(e):  # delphinidin || new
   return ldel(e) and (E18 in e)

def gc(e):  # gallocatechin || mh thesis
   return ldel(e) and (E20 in e)

def egt(e):  # epigallocatechin || mh thesis
   return dlm(e) and (E19 in e)

def gen(e):  # genistein || mh thesis
   return nar(e) and and_in(e, E21, E22)

def bun(e):  # butein || mh thesis
   return wca(e) and (or_in(e, E05, E09))

def butn(e):  # butin || new
   return bun(e) and (E10 in e)

def hcc(e):  # isoliquiritigenin || mh thesis
   return bun(e)

def dfv(e):  # liquiritigenin || new
   return hcc(e) and (E10 in e)

def gban(e):  # garbanzol || new
   return dfv(e) and (E15 in e)

def fstn(e):  # dihydrofisetin || new
   return (gban(e) and (E14 in e)) or (butn(e) and (E15 in e))

def df74(e):  # 7,4'-dihydroxyflavone || new
   return dfv(e) and (or_in(e, E11, E12))

def tiso(e):  # 2,7,4'-Trihydroxyisoflavanone || new
   return dfv(e) and (E21 in e)

def daid(e):  # daidzein || new
   return tiso(e) and (E22 in e)

def tnon(e):  # 6,7,4'-Trihydroxyflavanone || new
   return dfv(e) and (E08 in e)

def tet2(e):  # 2,6,7,4'-Tetrahydroxyisoflavanone || new
   return tnon(e) and (E21 in e)

def hdai(e):  # 6-Hydroxydaidzein || new
   return tet2(e)

def ggt(e):  # Glycosaminoglycan galactosyltransferase || enzyme
   return E_GGT in e

def dec(e):  # Deleted entry || enzyme
   return E_DEC in e

def soa(e):  # Serine O-acetyltransferase || enzyme
   return E_SOA in e

def v1g(e):  # vanillate 1-glucosyltransferase || enzyme
   return E_V1G in e
