# -*- coding: utf-8 -*-
"""
Originally Created on Wed Apr 17 2019

@author: vmoorman and Jordan Wilson
Made using Python 2.7
KEGG.py - JW started the code to get information from KEGG about the species we were interested in - April 2019
v0p1 - VRM getting the iteration and parsing part of the code to work - April 2019
v0p2 - VRM making the count code work - April 2019
v0p3 - JW getting the file output to work, revising the species codes and adding a dictionary - April 2019
v0p4 - JW creating directory change function for Gene_data and adding additional plant species - multiple things still broken - May 2019
v0p5 - VRM cleaned up some duplicate issues, folder locations, and ensured that species are written out in full  - May 2019
v0p6 - JW creating function for master FASTA file (error occurs when running every entry) - May 2019
v0p7 - JW editing lists for master fasta function and updating fasta function - May 2019
v0p8 - VRM fixing and creating gene list for master fasta file - June 2019
v0p8p2 - JW made a fasta file for each EC#- June 2019
v1p0 - VRM cleaned up code
v1p1p1 - JW added ReadMe
v1p2 - JW added function to get plants that make Epicatechin - incorrectly
v1p3 - VRM coded in epicatechin, catechin, eriodictyol, luteolin, naringenin listed in Chemical_Data --- leutolin actually wrong
v1p3p1 - VRM tried to fix leutolin


"""

'''This are required for the script to work - it has all of the required dependencies'''
#bioservices must be installed on your computer 
from bioservices.kegg import KEGG
import os
import sys
import datetime
k = KEGG()
now=datetime.datetime.now()
'''Make the necessary folders'''
foldername = os.getcwd()+ "\\" + raw_input("Input a save folder name: ")#gets current working directory and asks for user input for a new foldername
if os.path.exists(foldername): #determines if the new foldername already exists
    decision = raw_input("Warning, this folder already exists. Press the return key to continue anyway, or type anything to try again: ")
    if decision == "": #if return key hit will continue running th code, will overwrite anything in the prexisting folder
        pass
    else:
        sys.exit("Try an unused folder name next time") #stops code completely and the code will need to be restarted with a different foldername
genefoldername = foldername+"\Gene_Data" #variable used to create the folder for the gene data
fastafoldername = foldername+"\FASTA_Data" #variable used to create the folder for the fasta data
chemicalfoldername=foldername+"\Chemical_Data"
try: os.mkdir(foldername) #creates the main working folder for the code
except WindowsError: pass #stops code from making folder if this error occurs
try: os.mkdir(genefoldername) #creates folder were gene data files will be inserted
except WindowsError: pass
try: os.mkdir(fastafoldername) #creates folder where fasta files will be inserted
except WindowsError: pass
try: os.mkdir(chemicalfoldername) #creates folder where fasta files will be inserted
except WindowsError: pass


'''Python code to remove duplicate elements --- needs to be up here because of testing code '''
def Remove(listwithduplicates): 
    listwithoutduplicates = [] #creates an empty list
    for item in listwithduplicates: 
        if item not in listwithoutduplicates: 
            listwithoutduplicates.append(item) #adds item to empty list if it's not already in the list
    return listwithoutduplicates 

'''These are the list of the codes that you need to iterate over. 
 If you are needing just one code, just have one item in the list, but keep it as a list.
 You must have both species codes and pathway codes for this script to work.'''
speciescode_list = ["ats", "atr", "aly", "ath", "adu", "aip", "aof", "apro", "bpg", "bvg", "bdi", "bna", 
                    "boe", "brp", "ccaj", "csat", "crb", "cann", "cpap", "cqi", "cre", "cvr", "ccp", "cam", "cic", "cit", 
                    "csl", "cmo", "csv", "cmax", "cmos", "cpep", "cme", "ccav", "dcr", "dct", "dzi", "egu", "egr", 
                    "eus", "fve", "gsl", "gmx", "gab", "ghi", "gra", "han", "hbr", "ini", "jcu", "jre", 
                    "lsv", "lja", "lang", "mdm", "mesc", "mis", "mpp", "mcha", "mng", "mus", "nnu", "nau", "nsy", "nta", 
                    "nto", "oeu", "obr", "dosa", "osa", "olu", "ota", "psom", "peq", "pvu", "pda", "ppp", "pop", "peu", "pavi", "pmum", 
                    "pper", "pxb", "qsu", "rcu", "smo", "sind", "sita", "sly", "spen", "sot", "sbi", 
                    "soe", "thj", "tcc", "var", "vra", "vvi", "vcn", "zma", "zju"]#species codes for plants of interest,  "mtr" cut due to code errors
speciescode_dictionary = {"ats": "Aegilops tauschii", "atr": "Amborella trichopoda", "aly": "Arabidopsis lyrata",
                          "ath": "Arabidopsis thaliana", "adu": "Arachis duranensis", "aip": "Arachis ipaensis",
                          "aof": "Asparagus officinalis", "apro": "Auxenochlorella protothecoides",
                          "bpg": "Bathycoccus prasinos", "bvg": "Beta vulgaris", "bdi": "Brachypodium distachyon",
                          "bna": "Brassica napus", "boe": "Brassica oleracea", "brp": "Brassica rapa",
                          "ccaj": "Cajanus cajan", "csat": "Camelina sativa", "crb": "Capsella rubella",
                          "cann": "Capsicum annuum", "cpap": "Carica papaya", "cqi": "Chenopodium quinoa",
                          "cre": "Chlamydomonas reinhardtii", "cvr": "Chlorella variabilis", "ccp": "Chondrus crispus", "cam": "Cicer arietinum", 
                          "cic": "Citrus clementina", "cit": "Citrus sinensis", "csl": "Coccomyxa subellipsoidea", "cmo": "Cucumis melo", 
                          "csv": "Cucumis sativus", "cmax": "Cucurbita maxima", "cmos": "Cucurbita moschata", 
                          "cpep": "Cucurbita pepo subsp. pepo", "cme": "Cyanidioschyzon merolae",  "ccav": "Cynara cardunculus var. scolymus",
                          "dcr": "Daucus carota", "dct": "Dendrobium catenatum", "dzi": "Durio zibethinus", 
                          "egu": "Elaeis guineensis", "egr": "Eucalyptus grandis", "eus": "Eutrema salsugineum", 
                          "fve": "Fragaria vesca", "gsl": "Galdieria sulphuraria", "gmx": "Glycine max", "gab": "Gossypium arboreum", "ghi": "Gossypium hirsutum", 
                          "gra": "Gossypium raimondii", "han": "Helianthus annuus", "hbr": "Hevea brasiliensis", 
                          "ini": "Ipomoea nil", "jcu": "Jatropha curcas", "jre": "Juglans regia", 
                          "lsv": "Lactuca sativa", "lja": "Lotus japonicus", "lang": "Lupinus angustifolius", "mdm": "Malus domestica", 
                          "mesc": "Manihot esculenta", "mis": "Micromonas commoda",
                          "mpp": "Micromonas pusilla", "mcha": "Momordica charantia", "mng": "Monoraphidium neglectum", "mus": "Musa acuminata", 
                          "nnu": "Nelumbo nucifera", "nau": "Nicotiana attenuata", "nsy": "Nicotiana sylvestris", "nta": "Nicotiana tabacum", 
                          "nto": "Nicotiana tomentosiformis", "oeu": "Olea europaea v. sylvestris", "obr": "Oryza brachyantha", 
                          "dosa": "Oryza sativa japonica (RAPDB)", "osa": "Oryza sativa japonica (RefSeq)",
                          "olu": "Ostreococcus lucimarinus", "ota": "Ostreococcus tauri", "psom": "Papaver somniferum",
                          "peq": "Phalaenopsis equestris", "pvu": "Phaseolus vulgaris", "pda": "Phoenix dactylifera",
                          "ppp": "Physcomitrella patens subsp. Patens", "pop": "Populus trichocarpa", "peu": "Populus euphratica",
                          "pavi": "Prunus avium", "pmum": "Prunus mume", "pper": "Prunus persica", 
                          "pxb": "Pyrus x bretschneideri", "qsu": "Quercus suber", "rcu": "Ricinus communis", 
                          "smo": "Selaginella moellendorffii", "sind": "Sesamum indicum", "sita": "Setaria italica", 
                          "sly": "Solanum lycopersicum", "spen": "Solanum pennellii", "sot": "Solanum tuberosum", 
                          "sbi": "Sorghum bicolor", "soe": "Spinacia oleracea", "thj": "Tarenaya hassleriana", 
                          "tcc": "Theobroma cacao", "var": "Vigna angularis", "vra": "Vigna radiata", 
                          "vvi": "Vitis vinifera", "vcn": "Volvox carteri f. nagariensis", "zma": "Zea mays", "zju": "Ziziphus jujuba"}#Dictionary that defines the appropriate genus species for each code
pathwaycode_list = ["00940", "00941", "00942", "00943", "00944"]#list of pathways of interest
pathwaycode_dictionary = {"00940": "phenylpropanoids", "00941": "flavonoids", "00942": "anthocyanins",  
                          "00943": "isoflavonoids", "00944": "flavones/flavonols", "00945": "stilbenoids"}#Dictionary defining each pathway by each chemical they're responisble for, "mtr": "Medicago truncatula" cut due to code error

pathwayID_list = [i+j for i in speciescode_list for j in pathwaycode_list] #this is the full list of every pathway and species from the above lists
'''pathwayID_list = ["ats00940", "ats00941", "brp00940", "brp00941", "aip00940", "aip00941"] #use for testing only; comment out when actually running the script
speciescode_list=[] #use for testing only; comment out when actually running the script
for i in pathwayID_list: speciescode_list.append(filter(str.isalpha, i)) #use for testing only; comment out when actually running the script
speciescode_list = Remove(speciescode_list) #use for testing only; comment out when actually running the script'''



'''defining the function that will actually get all of the data that is required'''      
def gene_pathway_data(pathwayID):
    #print "Looking up data from: "+pathwayID #good for testing but takes up a lot of time in actuality
    entrylines_list = k.get(pathwayID).split("\n") #gets all of the data and splits it by line
    #print genes
    linecount=0
    for line in entrylines_list: #finds the places that have the genes listed
        entrylines_list[linecount]=line.strip()#removes extra unicode/spaces and replaces the original entry
        #print str(linecount) + ": " + line
        if line.startswith('GENE'): #finds where GENE is at in each entry
            #print line 
            entrylines_list[linecount]=line.replace("GENE","").strip()#replaces GENE with a blank and removes extra spaces at the beginning and end
            GENElocator=linecount #gene locator is now the item in the list that has GENE
        if line.startswith('COMPOUND'): #finds where COMPOUND is in the entry
            #print line
            COMPOUNDlocator=linecount#compoundlocator is now the item in the list that begins with compound
        linecount+=1
    geneline_list= entrylines_list[GENElocator:COMPOUNDlocator] #makes a list that is just the gene entry lines
    linecount=0
    for i in geneline_list: #this section makes a list of lists that are appreciately separated
        i= i.replace("  ","^*^").replace(";","^*^").replace("[","^*^").replace("]","")#makes ^*^ the signifier for splitting
        i=i.split("^*^")
        i.insert(0,speciescode_dictionary[filter(str.isalpha, pathwayID)])#pathwayID[0:3]])#replaces the species code with the genus species from the dictionary value for each key 
        jcount=0
        for j in i: #cleans up list of lists of extra spaces at the beginning and end of each list
            i[jcount]=j.strip()
            jcount+=1
        geneline_list[linecount]= i #replaces the list with the new cleaned list of lists
        linecount+=1#iterates through each list in the entry
    return geneline_list
        
'''function that saves each file with a name that includes pathwayID  using the data from genes_listoflists'''
def get_lists(whatlist, outname, outfolder=os.getcwd()): #whatlist should be a list of lists, outname should include an appropriate extension
    os.chdir(outfolder)#changes directory to current working directory
    writedoc = file(outname,"w") #gives write privilages for the function to write to the file
    for line in whatlist:
        for item in line:
            item=str(item).replace("\n","") #removes the new lines in each list of list
            if item == "":
                writedoc.write("-") #if the list in the list of list is empty writes a dash
            else:
                writedoc.write(item) #writes the entry in the list of lists to the file
            writedoc.write(",")#tab deliniated; use "," for csv files
        writedoc.write("\n") #creates a new line after the above portion of teh function
    writedoc.write("\n")
    writedoc.close() #closes file, required
    print  outname, "saved in: ", outfolder

'''iterate over pathwayID_list and use the fuction defined above'''
masterlist=[]
for pathwayID in pathwayID_list:
    notpresent=0
    try: currentlist=gene_pathway_data(pathwayID)  #need to ignore everything if there is no pathway for that species
    except AttributeError: 
        #print " No data in "+ pathwayID #good for testing but takes up a lot of time in actuality
        notpresent=1
    if notpresent==0: 
        masterlist.extend(currentlist)
    
    #run the function that saves each file    
    notpresent2=0
    try: get_gene_lists = get_lists(gene_pathway_data(pathwayID), "Gene_data_"+pathwayID+".csv", genefoldername) #try actually does it if it works
    except AttributeError:
        #print " No data in "+ pathwayID #good for testing but takes up a lot of time in actuality
        notpresent2=1

'''remove duplicates'''    
#print masterlist
masterlist_nodup=Remove(masterlist)
masterlist_nodup = list(filter(None,masterlist_nodup)) #removes false values from masterlist_nodup and turns it into a list
count=0
for i in masterlist_nodup: #removes false values and iterates through the list of lists
    masterlist_nodup[count]=list(filter(None,masterlist_nodup[count]))
    count+=1
    

'''find unique EC numbers so have a generic function and run it'''
def UniqueElementList(listname, locationnumber): #be careful as there are cases of one less item - use "last" to fix that problem here
    originallocationnumber=locationnumber
    ElementList=[]
    for i in listname:
        #print 
        if originallocationnumber == "last": locationnumber=len(i)-1 #assigns the string "last" to the very last list in the list of lists
        #print locationnumber
        if i[int(locationnumber)] not in ElementList: #finds unique EC number not in the list
            ElementList.append(i[int(locationnumber)]) #adds it to the list
    return ElementList
EC_list=UniqueElementList(masterlist_nodup, "last")
#print EC_list
#print speciescode_list

'''creating the matrix and adding up the counts'''
mastercount_list = [["Species"]]
mastercount_list[0].extend(EC_list) #adds the Unique EC numbers to the end of the mastercount list of lists
for i in speciescode_list: #@#
    mastercount_list.append([i]) #first item in each row (but first) is the speciescode
#print EC_list
#print mastercount_list

icount=0
for i in mastercount_list[0]: #for each EC number
    #print "i: "+ i
    if icount != 0: #first column isn't actually an EC#,
        EC=i
        #print EC
        jcount=0
        for j in mastercount_list: #for each species (using the speciecode)
            #print j
            if jcount != 0: #first row isn't actually a species
                #print j
                species=speciescode_dictionary[j[0]] 
                #print species
                lcount=0
                for l in masterlist_nodup: #iterate over the culled masterlist to check for matching sets
                    if l[0]==species and l[len(l)-1]==EC:
                        lcount+=1 
                mastercount_list[jcount].append(lcount)
            #else: print "help"
            jcount+=1 ####  
    icount+=1
#print mastercount_list#[0:3] #print this to test, but hide when actually running

#change mastercount_list to be actual species:
icount=0
for i in mastercount_list:
    if icount != 0: mastercount_list[icount][0]=speciescode_dictionary[i[0]] #replaces species code with genus specie names
    icount+=1

    

'''Make Master Files'''
get_masterlist = get_lists(masterlist_nodup, "Master_List.csv", foldername)
get_mastercount = get_lists(mastercount_list, "Master_Count.csv", foldername)



'''make a master fasta file saved in foldername'''
rev_dict={v:k for k,v in speciescode_dictionary.iteritems()}#reverses dictionary keys and values
genelist_frommaster=[]
for i in masterlist_nodup:
    #print rev_dict[i[0]]+":"+i[1]
    genelist_frommaster.append(rev_dict[i[0]]+":"+i[1]) #combines species codes and gene numbers in a list to be used for the master fasta function
 
def get_master_fasta(gene):
    DNA_info_list=[]
    for gene in genelist_frommaster:
        #print gene
        gene_fasta_data = k.get(gene).split("\n") #calls the entry from KEGG and splits it into new lines
        linecount = 0
        for line in gene_fasta_data:
            gene_fasta_data[linecount]=line.strip() #removes blank spaces at the beginning and end of each line
            if line.startswith('ORGANISM'): #finds where the entry that begins with organism is
                gene_fasta_data[linecount]=line.replace("ORGANISM","").strip() #remoes organism and removes blank spaces ta beginning and end
                find_blankspace=gene_fasta_data[linecount].find(" ") #finds where thedouble blank space is in the organism line
                organism_name= gene_fasta_data[linecount][find_blankspace:] #the genus specie name is left from the orginal entry
            linecount+=1
        linecount=0
        for line in gene_fasta_data:
            gene_fasta_data[linecount]=line.strip()
            if line.startswith('ORTHOLOGY'):
                gene_fasta_data[linecount]=line.strip
                find_EC=line.find("EC:") #finds within the line where the EC number is
                gene_fasta_data[linecount]=line[find_EC:-1].replace("[", "").replace("EC:", "") #removes the beginning bracket and EC: leaving just the number
                EC_number="EC "+ gene_fasta_data[linecount] #adds EC back
                joined_organism_EC=[">"+str(organism_name).strip()+"%"+EC_number+"%"+gene.split(":")[1]]#adds > to beginning to find the beginning of each entry more easily, adds % between EC number and the rest of the entry to help separate the EC number for later and removes semicolon from the gene entry and adds the gene number
                DNA_info_list.append(joined_organism_EC) #adds the entry to the blank list
            linecount+=1
        linecount=0
        for line in gene_fasta_data:
            gene_fasta_data[linecount]=line.strip()
            if line.startswith('NTSEQ'):
                gene_fasta_data[linecount]=line.replace("NTSEQ","").strip()
                NTSEQlocator=linecount
            linecount+=1
        DNA_data_list=gene_fasta_data[NTSEQlocator:]
        DNA_Seq=DNA_data_list[1:len(DNA_data_list)-2] #Takes just the DNA sequence
        separator=""
        joined_DNA_seq=[separator.join(DNA_Seq)] #combines the separate DNA sequence lines into one string and turns that into a single entry list
        DNA_info_list.append(joined_DNA_seq) #adds single entry list to the list of lists
    return DNA_info_list
              
get_lists(get_master_fasta(genelist_frommaster), "Master_FASTA.csv", fastafoldername)

masterfasta=get_master_fasta(genelist_frommaster)#should actually go above the first time it gets called

'''make a fasta file for each EC number saved in fastafoldername'''
ECorderlist=[]
fasta_byEC=[]
icount=0
for i in masterfasta:
    #print i
    if i[0].startswith(">"):
        isplit=i[0].split("%") #finds the EC numbers using the % added previously
        #print isplit
        EC=isplit[1]
        EC_tf="false"
        jcount=0
        for j in ECorderlist: #if the EC number is already present sets it to true and continues
            if EC==j:
                ECcount=jcount
                EC_tf="true"
            jcount+=1
        if EC_tf =="false": #if the EC number is not there, it will be added to the list
            ECorderlist.append(EC)
            fasta_byEC.append(masterfasta[icount:icount+2])
        else:
            fasta_byEC[ECcount].extend(masterfasta[icount:icount+2])
        
    icount+=1

#print  ECorderlist
#print fasta_byEC
    
'''creates the FASTA files by EC numbers'''
icount=0
for i in ECorderlist:
    name=i.replace(".","p").replace("EC ","").replace(" " ,"")
    FASTAbyEC= get_lists(fasta_byEC[icount], name+".csv", fastafoldername)
    print name+".csv" + " saved in: "+fastafoldername
    icount+=1

'''Creates ReadMe file'''     
with open(foldername+"\ReadMe.txt", "w") as ReadMe:
    ReadMe.write("KEGG_v1p1.py\n")
    ReadMe.write(now.strftime("%m-%d-%Y")+"\n")
    ReadMe.write(foldername+"\n")
    ReadMe.write("This script creates a series of files related to the genes associated with plant flavonoids from various species of plants. This script first creates the MasterCount and MasterList files; the MasterCount counts the number genes each plant species have that correspond with each EC number; while the MasterList lists every gene with number for each plant specie. These are located in "+os.getcwd())
    ReadMe.write(". The script also creates files that only contains the genes of a single plant species biochemical pathway which are located in "+genefoldername+". The script also creates a Master FASTA files which contains the DNA sequence of each gene and FASTA files organized by EC number, these are located in "+fastafoldername)
    ReadMe.close
    
'''Gets plants that make epicatechin'''

def ECandor(listname):
    icount=0
    codestring=''
    for i in listname:
        if icount>0:
            codestring+='"'+i+'"'
            try: 
                listname[icount+1]
            except IndexError: 
                break
            codestring+=' '+listname[0]+' '
        icount+=1
    #print codestring
    return codestring

masterEC_list = [["species","EC#s"]]
icount=0
for i in mastercount_list: #for each species
    #print "species :", i[0]
    speciesEClist=[]
    if icount==0:
        pass
        
    else:
        jcount=0
        for j in i: #for EC in species
            if jcount==0:
                speciesEClist.append(j)
                print j
            else:
                if str(j)=="0": 
                    print " "+ str(j)+ " = no enzyme"
                    pass
                else: 
                    print mastercount_list[0][jcount] + " = " + str(j)
                    speciesEClist.append(mastercount_list[0][jcount])
            jcount+=1    
    masterEC_list.append(speciesEClist)
    icount+=1
    
print masterEC_list
masterEC_list=masterEC_list[1:]

phenylalanineTOcinnamicacid=["or","EC:4.3.1.24","EC:4.3.1.25"]
cinnamicacidTOpcoumaroyllcoa=["and","EC:6.2.1.12","EC:1.14.14.91"]
pcoumaroyllcoaTOcaffeoylcoa1=["and","EC:1.14.13.-"]
pcoumaroyllcoaTOcaffeoylcoa2=["and","EC:2.3.1.133","EC:1.14.14.96"]
pcoumaroyllcoaTOnaringenin=["and","EC:2.3.1.74","EC:5.5.1.6"]
naringeninTOeriodictyol=["or","EC:1.14.14.81","EC:1.14.14.82"]
caffeoylcoaTOeriodictyol=["and","EC:2.3.1.74"]
eriodictyolTOleucocyanidin=["and","EC:1.14.11.9","EC:1.1.12.19"]
leucocyanidinTOcatechin=["and","EC:1.17.1.3"]
leucocyanidinTOcyanidin=["and","EC:1.14.20.4"]
cyanidinTOepicatechin=["and","EC:1.3.1.77"]
pcoumaroyllcoaTOnaringenin=["and","EC:2.3.1.74","EC:5.5.1.6"]
eriodictyolTOluteolin=["or","EC:1.14.20.5", "EC:1.14.19.76"]
naringeninTOapigenin=["or","EC:1.14.20.5","EC:1.14.19.76"]
apigeninTOluteolin=["or","EC:1.14.14.81","EC:1.14.14.82"]


epicatechinlist=[]
catechinlist=[]
eriodictyollist=[]
naringeninlist=[]
luteolinlist=[]
#Be careful in making these of parenthases
for i in masterEC_list:
    #print i
    
    epicatechin= "if (("+ECandor(phenylalanineTOcinnamicacid) + ") and " + ECandor(cinnamicacidTOpcoumaroyllcoa)+ " and (" + ECandor(pcoumaroyllcoaTOcaffeoylcoa1)+" or ("+ECandor(
        pcoumaroyllcoaTOcaffeoylcoa2)+"))"+" and "+ECandor(caffeoylcoaTOeriodictyol)+" and "+ECandor(eriodictyolTOleucocyanidin)+" and "+ECandor(
        leucocyanidinTOcyanidin)+" and "+ECandor(cyanidinTOepicatechin)+") in i: epicatechinlist.append([i[0]])"
    exec epicatechin   
    
    catechin= "if (("+ECandor(phenylalanineTOcinnamicacid) + ") and " + ECandor(cinnamicacidTOpcoumaroyllcoa)+ " and (" + ECandor(pcoumaroyllcoaTOcaffeoylcoa1)+" or ("+ECandor(
        pcoumaroyllcoaTOcaffeoylcoa2)+"))"+" and "+ECandor(caffeoylcoaTOeriodictyol)+" and "+ECandor(eriodictyolTOleucocyanidin)+" and "+ECandor(
        leucocyanidinTOcatechin)+") in i: catechinlist.append([i[0]])"
    exec catechin

    eriodictyol= "if (("+ECandor(phenylalanineTOcinnamicacid) + ") and " + ECandor(cinnamicacidTOpcoumaroyllcoa)+ " and (" + ECandor(pcoumaroyllcoaTOcaffeoylcoa1)+" or ("+ECandor(
        pcoumaroyllcoaTOcaffeoylcoa2)+"))"+" and "+ECandor(caffeoylcoaTOeriodictyol)+") in i: eriodictyollist.append([i[0]])"
    exec eriodictyol
    
    luteolin= "if (("+ECandor(phenylalanineTOcinnamicacid) + ") and " + ECandor(cinnamicacidTOpcoumaroyllcoa)+ " and (" + ECandor(pcoumaroyllcoaTOcaffeoylcoa1)+" or ("+ECandor(
        pcoumaroyllcoaTOcaffeoylcoa2)+"))"+" and "+ECandor(caffeoylcoaTOeriodictyol)+" and ("+ECandor(eriodictyolTOluteolin)+")) in i: luteolinlist.append([i[0]])"
    exec luteolin
    '''luteolin= "if (("+ECandor(phenylalanineTOcinnamicacid) + ") and " + ECandor(cinnamicacidTOpcoumaroyllcoa)+ " and (" + ECandor(pcoumaroyllcoaTOcaffeoylcoa1)+" or ("+ECandor(
        pcoumaroyllcoaTOcaffeoylcoa2)+"))"+" and "+ECandor(caffeoylcoaTOeriodictyol)+" and "+ECandor(eriodictyolTOluteolin)+") in i: luteolinlist.append([i[0]])"
    exec luteolin'''#didn't work

    naringenin= "if (("+ECandor(phenylalanineTOcinnamicacid) + ") and " + ECandor(cinnamicacidTOpcoumaroyllcoa) + " and "+ ECandor(pcoumaroyllcoaTOnaringenin) +") in i: naringeninlist.append([i[0]])"
    exec naringenin

get_lists(epicatechinlist, "epicatechinspecies.txt", outfolder=foldername+"\Chemical_Data")
get_lists(catechinlist, "catechinspecies.txt", outfolder=foldername+"\Chemical_Data")
get_lists(eriodictyollist, "eriodictyolspecies.txt", outfolder=foldername+"\Chemical_Data")
get_lists(luteolinlist, "luteolinspecies.txt", outfolder=foldername+"\Chemical_Data")
get_lists(naringeninlist, "naringeninspecies.txt", outfolder=foldername+"\Chemical_Data")
    