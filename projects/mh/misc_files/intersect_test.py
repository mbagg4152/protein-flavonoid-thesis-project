import sys
import os

kmp = ['Aegilops tauschii', 'Amborella trichopoda', 'Arabidopsis lyrata', 'Arabidopsis thaliana',
       'Arachis duranensis', 'Arachis ipaensis', 'Asparagus officinalis', 'Beta vulgaris',
       'Brachypodium distachyon', 'Brassica napus', 'Brassica oleracea', 'Brassica rapa',
       'Camelina sativa', 'Capsicum annuum', 'Carica papaya', 'Chenopodium quinoa',
       'Cicer arietinum', 'Citrus clementina', 'Citrus sinensis', 'Cucumis melo', 'Cucumis sativus',
       'Cucurbita maxima', 'Cucurbita moschata', 'Cucurbita pepo subsp. pepo',
       'Cynara cardunculus var. scolymus', 'Daucus carota', 'Durio zibethinus', 'Elaeis guineensis',
       'Eucalyptus grandis', 'Eutrema salsugineum', 'Fragaria vesca', 'Glycine max', 'Glycine soja',
       'Gossypium arboreum', 'Gossypium hirsutum', 'Gossypium raimondii', 'Helianthus annuus',
       'Hevea brasiliensis', 'Ipomoea nil', 'Jatropha curcas', 'Juglans regia', 'Lactuca sativa',
       'Lotus japonicus', 'Lupinus angustifolius', 'Malus domestica', 'Manihot esculenta',
       'Medicago truncatula', 'Momordica charantia', 'Musa acuminata', 'Nelumbo nucifera',
       'Nicotiana attenuata', 'Nicotiana sylvestris', 'Nicotiana tabacum',
       'Nicotiana tomentosiformis', 'Olea europaea v. sylvestris', 'Oryza brachyantha',
       'Oryza sativa japonica RAPDB', 'Oryza sativa japonica RefSeq', 'Papaver somniferum',
       'Phalaenopsis equestris', 'Phaseolus vulgaris', 'Phoenix dactylifera', 'Populus euphratica',
       'Populus trichocarpa', 'Prunus avium', 'Prunus mume', 'Prunus persica',
       'Pyrus x bretschneideri', 'Quercus suber', 'Raphanus sativus', 'Ricinus communis',
       'Rosa chinensis', 'Sesamum indicum', 'Setaria italica', 'Solanum lycopersicum',
       'Solanum pennellii', 'Solanum tuberosum', 'Spinacia oleracea', 'Tarenaya hassleriana',
       'Theobroma cacao', 'Vigna angularis', 'Vigna radiata', 'Vigna unguiculata', 'Vitis vinifera',
       'Zea mays', 'Ziziphus jujuba']
myc = ['Aegilops tauschii', 'Brachypodium distachyon', 'Capsicum annuum', 'Carica papaya',
       'Cicer arietinum', 'Citrus clementina', 'Citrus sinensis', 'Durio zibethinus',
       'Elaeis guineensis', 'Eucalyptus grandis', 'Glycine max', 'Glycine soja',
       'Gossypium arboreum', 'Gossypium hirsutum', 'Gossypium raimondii', 'Hevea brasiliensis',
       'Jatropha curcas', 'Juglans regia', 'Lotus japonicus', 'Lupinus angustifolius',
       'Manihot esculenta', 'Medicago truncatula', 'Musa acuminata', 'Nelumbo nucifera',
       'Nicotiana attenuata', 'Nicotiana sylvestris', 'Nicotiana tabacum',
       'Nicotiana tomentosiformis', 'Oryza brachyantha', 'Oryza sativa japonica RAPDB',
       'Oryza sativa japonica RefSeq', 'Phalaenopsis equestris', 'Phaseolus vulgaris',
       'Phoenix dactylifera', 'Populus euphratica', 'Populus trichocarpa', 'Quercus suber',
       'Setaria italica', 'Solanum lycopersicum', 'Solanum pennellii', 'Solanum tuberosum',
       'Theobroma cacao', 'Vigna angularis', 'Vigna radiata', 'Vigna unguiculata', 'Vitis vinifera',
       'Zea mays', 'Ziziphus jujuba']
que = ['Aegilops tauschii', 'Amborella trichopoda', 'Arabidopsis lyrata', 'Arabidopsis thaliana',
       'Arachis duranensis', 'Arachis ipaensis', 'Asparagus officinalis', 'Beta vulgaris',
       'Brachypodium distachyon', 'Brassica napus', 'Brassica oleracea', 'Brassica rapa',
       'Camelina sativa', 'Capsicum annuum', 'Carica papaya', 'Chenopodium quinoa',
       'Cicer arietinum', 'Citrus clementina', 'Citrus sinensis', 'Cucumis melo', 'Cucumis sativus',
       'Cucurbita maxima', 'Cucurbita moschata', 'Cucurbita pepo subsp. pepo',
       'Cynara cardunculus var. scolymus', 'Daucus carota', 'Durio zibethinus', 'Elaeis guineensis',
       'Eucalyptus grandis', 'Eutrema salsugineum', 'Fragaria vesca', 'Glycine max', 'Glycine soja',
       'Gossypium arboreum', 'Gossypium hirsutum', 'Gossypium raimondii', 'Helianthus annuus',
       'Hevea brasiliensis', 'Ipomoea nil', 'Jatropha curcas', 'Juglans regia', 'Lactuca sativa',
       'Lotus japonicus', 'Lupinus angustifolius', 'Malus domestica', 'Manihot esculenta',
       'Medicago truncatula', 'Musa acuminata', 'Nelumbo nucifera', 'Nicotiana attenuata',
       'Nicotiana sylvestris', 'Nicotiana tabacum', 'Nicotiana tomentosiformis',
       'Olea europaea v. sylvestris', 'Oryza brachyantha', 'Oryza sativa japonica RAPDB',
       'Oryza sativa japonica RefSeq', 'Papaver somniferum', 'Phalaenopsis equestris',
       'Phaseolus vulgaris', 'Phoenix dactylifera', 'Populus euphratica', 'Populus trichocarpa',
       'Prunus avium', 'Prunus mume', 'Prunus persica', 'Pyrus x bretschneideri', 'Quercus suber',
       'Raphanus sativus', 'Ricinus communis', 'Rosa chinensis', 'Sesamum indicum',
       'Setaria italica', 'Solanum lycopersicum', 'Solanum pennellii', 'Solanum tuberosum',
       'Spinacia oleracea', 'Tarenaya hassleriana', 'Theobroma cacao', 'Vigna angularis',
       'Vigna radiata', 'Vigna unguiculata', 'Vitis vinifera', 'Zea mays', 'Ziziphus jujuba']
kxn = ['Aegilops tauschii', 'Amborella trichopoda', 'Arachis duranensis', 'Asparagus officinalis',
       'Beta vulgaris', 'Brachypodium distachyon', 'Carica papaya', 'Chenopodium quinoa',
       'Cicer arietinum', 'Citrus clementina', 'Citrus sinensis', 'Durio zibethinus',
       'Elaeis guineensis', 'Eucalyptus grandis', 'Glycine max', 'Gossypium arboreum',
       'Gossypium hirsutum', 'Gossypium raimondii', 'Helianthus annuus', 'Hevea brasiliensis',
       'Jatropha curcas', 'Juglans regia', 'Lotus japonicus', 'Lupinus angustifolius',
       'Malus domestica', 'Manihot esculenta', 'Musa acuminata', 'Nelumbo nucifera',
       'Nicotiana attenuata', 'Nicotiana sylvestris', 'Nicotiana tabacum',
       'Nicotiana tomentosiformis', 'Olea europaea v. sylvestris', 'Oryza brachyantha',
       'Oryza sativa japonica RAPDB', 'Papaver somniferum', 'Phaseolus vulgaris',
       'Phoenix dactylifera', 'Populus trichocarpa', 'Populus euphratica', 'Prunus avium',
       'Prunus mume', 'Prunus persica', 'Pyrus x bretschneideri', 'Quercus suber',
       'Ricinus communis', 'Solanum pennellii', 'Solanum tuberosum', 'Spinacia oleracea',
       'Theobroma cacao', 'Vigna angularis', 'Vigna radiata', 'Vitis vinifera', 'Ziziphus jujuba',
       'Glycine soja', 'Medicago truncatula', 'Rosa chinensis', 'Vigna unguiculata']
egt = ['Aegilops tauschii', 'Carica papaya', 'Cicer arietinum', 'Citrus clementina',
       'Citrus sinensis', 'Durio zibethinus', 'Elaeis guineensis', 'Eucalyptus grandis',
       'Glycine max', 'Gossypium arboreum', 'Gossypium hirsutum', 'Gossypium raimondii',
       'Hevea brasiliensis', 'Jatropha curcas', 'Juglans regia', 'Lotus japonicus',
       'Lupinus angustifolius', 'Manihot esculenta', 'Musa acuminata', 'Nelumbo nucifera',
       'Nicotiana attenuata', 'Nicotiana sylvestris', 'Nicotiana tabacum',
       'Nicotiana tomentosiformis', 'Oryza brachyantha', 'Oryza sativa japonica RAPDB',
       'Phaseolus vulgaris', 'Phoenix dactylifera', 'Populus trichocarpa', 'Populus euphratica',
       'Quercus suber', 'Setaria italica', 'Solanum lycopersicum', 'Solanum pennellii',
       'Solanum tuberosum', 'Theobroma cacao', 'Vigna angularis', 'Vigna radiata', 'Zea mays',
       'Ziziphus jujuba', 'Glycine soja', 'Medicago truncatula', 'Vigna unguiculata']
ec = ['Aegilops tauschii', 'Amborella trichopoda', 'Arabidopsis lyrata', 'Arabidopsis thaliana',
      'Arachis duranensis', 'Asparagus officinalis', 'Beta vulgaris', 'Brassica napus',
      'Brassica oleracea', 'Brassica rapa', 'Camelina sativa', 'Carica papaya',
      'Chenopodium quinoa', 'Cicer arietinum', 'Citrus clementina', 'Citrus sinensis',
      'Durio zibethinus', 'Elaeis guineensis', 'Eucalyptus grandis', 'Eutrema salsugineum',
      'Glycine max', 'Gossypium arboreum', 'Gossypium hirsutum', 'Gossypium raimondii',
      'Hevea brasiliensis', 'Jatropha curcas', 'Juglans regia', 'Lotus japonicus',
      'Lupinus angustifolius', 'Malus domestica', 'Manihot esculenta', 'Musa acuminata',
      'Nelumbo nucifera', 'Nicotiana attenuata', 'Nicotiana sylvestris', 'Nicotiana tabacum',
      'Nicotiana tomentosiformis', 'Oryza brachyantha', 'Oryza sativa japonica RAPDB',
      'Papaver somniferum', 'Phaseolus vulgaris', 'Phoenix dactylifera', 'Populus trichocarpa',
      'Populus euphratica', 'Prunus avium', 'Prunus mume', 'Prunus persica',
      'Pyrus x bretschneideri', 'Quercus suber', 'Ricinus communis', 'Setaria italica',
      'Solanum lycopersicum', 'Solanum pennellii', 'Solanum tuberosum', 'Tarenaya hassleriana',
      'Theobroma cacao', 'Vigna angularis', 'Vigna radiata', 'Zea mays', 'Ziziphus jujuba',
      'Glycine soja', 'Medicago truncatula', 'Rosa chinensis', 'Vigna unguiculata',
      'Raphanus sativus']
gc = ['Aegilops tauschii', 'Brachypodium distachyon', 'Carica papaya', 'Cicer arietinum',
      'Citrus clementina', 'Citrus sinensis', 'Durio zibethinus', 'Elaeis guineensis',
      'Eucalyptus grandis', 'Glycine max', 'Gossypium arboreum', 'Gossypium hirsutum',
      'Gossypium raimondii', 'Hevea brasiliensis', 'Jatropha curcas', 'Juglans regia',
      'Lotus japonicus', 'Lupinus angustifolius', 'Manihot esculenta', 'Musa acuminata',
      'Nelumbo nucifera', 'Nicotiana attenuata', 'Nicotiana sylvestris', 'Nicotiana tabacum',
      'Nicotiana tomentosiformis', 'Oryza brachyantha', 'Oryza sativa japonica RAPDB',
      'Phaseolus vulgaris', 'Phoenix dactylifera', 'Populus trichocarpa', 'Populus euphratica',
      'Quercus suber', 'Solanum pennellii', 'Solanum tuberosum', 'Theobroma cacao',
      'Vigna angularis', 'Vigna radiata', 'Vitis vinifera', 'Ziziphus jujuba', 'Glycine soja',
      'Medicago truncatula', 'Vigna unguiculata']
agi = ['Aegilops tauschii', 'Arachis duranensis', 'Arachis ipaensis', 'Brachypodium distachyon',
       'Cajanus cajan', 'Cicer arietinum', 'Citrus clementina', 'Citrus sinensis',
       'Cynara cardunculus var. scolymus', 'Daucus carota', 'Durio zibethinus',
       'Eucalyptus grandis', 'Glycine max', 'Gossypium arboreum', 'Gossypium hirsutum',
       'Gossypium raimondii', 'Helianthus annuus', 'Hevea brasiliensis', 'Lactuca sativa',
       'Lupinus angustifolius', 'Manihot esculenta', 'Olea europaea v. sylvestris',
       'Oryza brachyantha', 'Oryza sativa japonica RAPDB', 'Oryza sativa japonica RefSeq',
       'Phaseolus vulgaris', 'Populus trichocarpa', 'Populus euphratica', 'Ricinus communis',
       'Sesamum indicum', 'Setaria italica', 'Sorghum bicolor', 'Theobroma cacao',
       'Vigna angularis', 'Zea mays', 'Glycine soja', 'Medicago truncatula', 'Vigna unguiculata']
lu2 = ['Aegilops tauschii', 'Arachis duranensis', 'Arachis ipaensis', 'Brachypodium distachyon',
       'Cajanus cajan', 'Cicer arietinum', 'Citrus clementina', 'Citrus sinensis',
       'Cynara cardunculus var. scolymus', 'Daucus carota', 'Durio zibethinus',
       'Eucalyptus grandis', 'Glycine max', 'Glycine soja', 'Gossypium arboreum',
       'Gossypium hirsutum', 'Gossypium raimondii', 'Helianthus annuus', 'Hevea brasiliensis',
       'Lactuca sativa', 'Lupinus angustifolius', 'Manihot esculenta', 'Medicago truncatula',
       'Olea europaea v. sylvestris', 'Oryza brachyantha', 'Oryza sativa japonica RAPDB',
       'Oryza sativa japonica RefSeq', 'Phaseolus vulgaris', 'Populus euphratica',
       'Populus trichocarpa', 'Ricinus communis', 'Sesamum indicum', 'Setaria italica',
       'Sorghum bicolor', 'Theobroma cacao', 'Vigna angularis', 'Vigna unguiculata', 'Zea mays']
sp = os.sep
ind = 'inter' + sp
esd = 'ec-stuff' + sp

def main():
    intersections()
    ec_ops()

def intersections():
    try: os.mkdir('inter')
    except (FileExistsError, OSError) as e: pass
    
    kmp_que_ec_chk = ''
    kmp_que = set(kmp).intersection(set(que))  # intersect of kmp & que
    kmp_que_only = list(set(kmp_que).difference(myc))  # intersect of kmp & que w/out vals from myc
    
    f3ols = set(ec).intersection(set(gc), set(egt), set(kxn))
    ec_vf3ols = set(ec).difference(f3ols)  # ec vs shared list of flav-3-ols
    just_ec_kxn = set(ec).intersection(kxn).difference(f3ols)  # ec, kxn v shared list of f-3-ols
    egt_vf3ols = set(egt).difference(f3ols)  # egt vs shared list of flav-3-ols
    gc_vf3ols = set(gc).difference(f3ols)  # gc vs shared list of flav-3-ols
    kxn_vf3ols = set(kxn).difference(f3ols)  # kxn vs shared list of flav-3-ols
    just_gc_egt = set(gc).intersection(egt).difference(f3ols)  # gc, egt v shared list of f-3-ols
    
    with open(ind + '_plant-ec-nums.tsv', 'r') as f:
        for item in f.readlines():
            name = item.split('\t')[0]
            if name in kmp_que_only:
                if '1.14.14.81' in item: kmp_que_ec_chk += name + ' YES 1.14.14.81'
                else: kmp_que_ec_chk += name + ' NO 1.14.14.81'
                kmp_que_ec_chk += '\n'
    
    write(ind + 'kmp-que-check.txt', kmp_que_ec_chk)
    write(ind + 'kmp-que-only.txt', '\n'.join(kmp_que_only))
    write(ind + 'flavan-3-ols.txt', '\n'.join(list(f3ols)))
    write(ind + 'ec-vs-f3ols.txt', '\n'.join(list(ec_vf3ols)))
    write(ind + 'egt-vs-f3ols.txt', '\n'.join(list(egt_vf3ols)))
    write(ind + 'gc-vs-f3ols.txt', '\n'.join(list(gc_vf3ols)))
    write(ind + 'kxn-vs-f3ols.txt', '\n'.join(list(kxn_vf3ols)))
    write(ind + 'just-gc-egt.txt', '\n'.join(list(just_gc_egt)))
    write(ind + 'just-ec-kxn.txt', '\n'.join(list(just_ec_kxn)))

def ec_ops():
    try: os.mkdir('ec-stuff')
    except (FileExistsError, OSError) as e: pass
    
    ec_nums = {}
    
    with open('inter' + sp + '_plant-ec-nums.tsv', 'r') as f:
        fdat = f.readlines()
        for item in fdat:
            name = item.split('\t')[0]
            for num in item.split('\t')[1:]:
                if num not in ec_nums.keys(): ec_nums[num] = []
        
        for key in ec_nums:
            for item in fdat:
                if key in item: ec_nums[key].append(item.split('\t')[0])
            
            fname = key.replace('EC:', '').replace('.', '-') + '.tsv'
            write(esd + fname, '\n'.join(sorted(ec_nums[key])))

def write(fname, content):
    with open(fname, 'w') as f: f.write(content)

def append(fname, content):
    with open(fname, 'a') as f: f.write(content)

if __name__ == '__main__':
    main()
