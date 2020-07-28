# -*- coding: utf-8 -*-

# bioservices must be installed on your computer.
# use command `pip3 install bioservices`

try:
    from bioservices.kegg import KEGG
except ImportError as e:
    print("you need to have bioservices installed to run this program.\n" +
          "install using command `pip3 install bioservices`")
    exit(0)

import os
import sys
import datetime
import json
from pathlib import Path

kegg = KEGG()
now = datetime.datetime.now()

# get active directory and user specified new directory name
userInput = input("Input a save folder name: ")
# prefixPath = os.getcwd()
# prefixPath = str(Path.home()) + "/Desktop"
newDirName = os.getcwd() + "/" + userInput

# determines if the new folder name already exists
if os.path.exists(newDirName):
    decision = input("Warning, this folder already exists." +
                     "Hit enter to continue anyway or try another:")
    # hitting return key removes previous directory contents
    if decision == "":
        pass
    # stops code completely and the code will need to be restarted with a different folder name
    else:
        sys.exit("Try an unused folder name next time")

# make path vars for gene, fasta & chemical data
geneDirPath = newDirName + '/Gene_Data'
fastaDirPath = newDirName + '/FASTA_Data'
chemDirPath = newDirName + '/Chemical_Data'
geneListFromMaster = []


def makeNewDirs():
    try:
        os.mkdir(newDirName)  # creates the main working folder for the code
    except OSError:
        print("could not make dir" + newDirName)  # stops code from making folder if this error occurs

    try:
        os.mkdir(geneDirPath)
    except OSError:
        print("could not make dir " + geneDirPath)

    try:
        os.mkdir(fastaDirPath)
    except OSError:
        print("could not make dir " + fastaDirPath)

    try:
        os.mkdir(chemDirPath)
    except OSError:
        print("could not make dir " + chemDirPath)


def getJsonData(fname, key):
    with open(fname) as jsonFile:
        data = json.load(jsonFile)
        # print(data[key])
        return data[key]


makeNewDirs()


# Python code to remove duplicate elements --- needs to be up here because of testing code
def removeDupes(hasDupes):
    uniqueList = []  # creates an empty list
    for item in hasDupes:
        if item not in uniqueList:
            uniqueList.append(item)  # adds item to empty list if it's not already in the list
    return uniqueList


# These are the list of the codes that you need to iterate over.
# If you are needing just one code, just have one item in the list, but keep it as a list.
# You must have both species codes and pathway codes for this script to work.
# species codes for plants of interest,  "mtr" cut due to code errors
speciesCodeList = getJsonData("species-list.json", 'list')
speciesCodeDict = getJsonData("species-pairs.json", 'pairs')
pathwayCodeList = getJsonData("path-codes.json", 'codes')  # list of pathways of interest
# Dictionary defining each pathway by each chemical they're responsible for
# "mtr": "Medicago truncatula" cut due to code error
pathwayCodeDict = getJsonData("path-code-pairs.json", 'pairs')

# this is the full list of every pathway and species from the above lists
pathwaySpeciesList = [i + j
                      for i in speciesCodeList
                      for j in pathwayCodeList]

# use for testing only; comment out when actually running the script
pathwayID_list = ["ats00940", "ats00941", "brp00940", "brp00941", "aip00940", "aip00941"]
speciesCodeList = []  # use for testing only; comment out when actually running the script
for i in pathwayID_list: speciesCodeList.append(filter(str.isalpha, i))
speciesCodeList = removeDupes(speciesCodeList)  # use for testing only; comment out when actually running the script


# defining the function that will actually get all of the data that is required
def genePathwayData(pathId):
    # initialize gene & compound locator since they're modified inside of loop
    geneLoc = None
    compoundLoc = None
    entryLines = kegg.get(pathId).split("\n")  # gets all of the data and splits it by line

    lineCount = 0

    for line in entryLines:  # finds the places that have the genes listed
        entryLines[lineCount] = line.strip()  # remove extra whitespace

        if line.startswith('GENE'):  # finds where GENE is at in each entry
            # remove 'GENE' & extra whitespace
            entryLines[lineCount] = line.replace("GENE", "").strip()
            geneLoc = lineCount  # gene locator is now the item in the list that has GENE

        if line.startswith('COMPOUND'):  # finds where COMPOUND is in the entry
            compoundLoc = lineCount  # compound locator is now the item in the list that begins with compound
        lineCount += 1

    geneLineList = entryLines[geneLoc:compoundLoc]  # makes a list that is just the gene entry lines
    lineCount = 0
    for i in geneLineList:  # this section makes a list of lists that are appreciately separated
        # makes ^*^ the signifier for splitting
        i = i.replace("  ", "^*^").replace(";", "^*^").replace("[", "^*^").replace("]", "")
        i = i.split("^*^")
        # pathwayID[0:3]])
        # replaces the species code with the genus species from the dictionary value for each key
        # print("path id: " + pathId)
        try:
            i.insert(0, speciesCodeDict[filter(str.isalpha, pathId)])
            # i.insert(0, pathwayCodeDict[filter(str.isalpha, pathId)])
        except KeyError:
            print("issue filtering. " + pathId + " not found in species code dictionary")
        jc = 0
        for j in i:  # cleans up list of lists of extra spaces at the beginning and end of each list
            i[jc] = j.strip()
            jc += 1
        geneLineList[lineCount] = i  # replaces the list with the new cleaned list of lists
        lineCount += 1  # iterates through each list in the entry
    return geneLineList


# function that saves each file with a name that includes pathwayID  using the data from genes_listoflists


# whatList should be a list of lists, outName should include an appropriate extension
def saveListsToFile(whatList, outName, outDir=os.getcwd()):
    os.chdir(outDir)  # changes directory to current working directory
    # os.chmod(outName, 222)  # give write permissions
    writedoc = open(outName, "w")

    for line in whatList:
        for item in line:
            item = str(item).replace("\n", "")  # removes the new lines in each list of list
            if item == "":
                writedoc.write("-")  # if the list in the list of list is empty writes a dash
            else:
                writedoc.write(item)  # writes the entry in the list of lists to the file
            writedoc.write(",")  # tab delineated; use "," for csv files
        writedoc.write("\n")  # creates a new line after the above portion of teh function
    writedoc.write("\n")
    writedoc.close()  # closes file, required
    print(outName, "saved in: ", outDir)


# iterate over pathwayID_list and use the function defined above
masterList = []
for pathwayID in pathwaySpeciesList:
    notPresent = 0
    currentList = None
    try:
        currentList = genePathwayData(pathwayID)  # need to ignore everything if there is no pathway for that species
    except AttributeError:
        # print " No data in "+ pathwayID #good for testing but takes up a lot of time in actuality
        notPresent = 1
    if notPresent == 0:
        masterList.extend(currentList)

    # run the function that saves each file
    notPresent2 = 0
    try:
        # try actually does it if it works
        saveListsToFile(genePathwayData(pathwayID), "Gene_data_" + pathwayID + ".csv", geneDirPath)
    except AttributeError:
        # print " No data in "+ pathwayID #good for testing but takes up a lot of time in actuality
        notPresent2 = 1

'''remove duplicates'''
masterListNoDupes = removeDupes(masterList)

# removes false values from masterListNoDupes and turns it into a list
masterListNoDupes = list(filter(None, masterListNoDupes))
count = 0
for i in masterListNoDupes:  # removes false values and iterates through the list of lists
    masterListNoDupes[count] = list(filter(None, masterListNoDupes[count]))
    count += 1

'''find unique EC numbers so have a generic function and run it'''


# be careful as there are cases of one less item - use "last" to fix that problem here
def uniqueElementList(listName, locationNum):
    originalLocNum = locationNum
    ElementList = []
    for i in listName:
        # assigns the string "last" to the very last list in the list of lists
        if originalLocNum == "last":
            locationNum = len(i) - 1

        # finds unique EC number not in the list
        if i[int(locationNum)] not in ElementList:
            ElementList.append(i[int(locationNum)])  # adds it to the list
    return ElementList


EC_list = uniqueElementList(masterListNoDupes, "last")

# creating the matrix and adding up the counts
masterCountList = [["Species"]]
masterCountList[0].extend(EC_list)  # adds the Unique EC numbers to the end of the master count list of lists
for i in speciesCodeList:  # @#
    masterCountList.append([i])  # first item in each row (but first) is the species code

iCount = 0
for i in masterCountList[0]:  # for each EC number
    if iCount != 0:  # first column isn't actually an EC#,
        EC = i
        jCount = 0
        for j in masterCountList:  # for each species (using the specie code)
            if jCount != 0:  # first row isn't actually a species
                species = speciesCodeDict[j[0]]
                kCount = 0
                for k in masterListNoDupes:  # iterate over the culled master list to check for matching sets
                    if k[0] == species and k[len(k) - 1] == EC:
                        kCount += 1
                masterCountList[jCount].append(str(kCount))
            jCount += 1
    iCount += 1

# change master count list to be actual species:
iCount = 0
for i in masterCountList:
    if iCount != 0: masterCountList[iCount][0] = speciesCodeDict[
        i[0]]  # replaces species code with genus specie names
    iCount += 1


# Make Master Files
def makeMasterFiles():
    saveListsToFile(masterListNoDupes, "Master_List.csv", newDirName)
    saveListsToFile(masterCountList, "Master_Count.csv", newDirName)


makeMasterFiles()


# make a master fasta file saved in folder name
def makeMasterFasta():
    rev_dict = {v: k for k, v in speciesCodeDict.items()}  # reverses dictionary keys and values
    for i in masterListNoDupes:
        # combines species codes and gene numbers in a list to be used for the master fasta function
        geneListFromMaster.append(rev_dict[i[0]] + ":" + i[1])


makeMasterFasta()


def getMasterFasta(gene):
    DNA_info_list = []
    for gene in geneListFromMaster:
        # print gene
        gene_fasta_data = kegg.get(gene).split("\n")  # calls the entry from KEGG and splits it into new lines
        linecount = 0
        for line in gene_fasta_data:
            gene_fasta_data[linecount] = line.strip()  # removes blank spaces at the beginning and end of each line
            if line.startswith('ORGANISM'):  # finds where the entry that begins with organism is
                gene_fasta_data[linecount] = line.replace("ORGANISM",
                                                          "").strip()  # remoes organism and removes blank spaces ta beginning and end
                find_blankspace = gene_fasta_data[linecount].find(
                    " ")  # finds where thedouble blank space is in the organism line
                organism_name = gene_fasta_data[linecount][
                                find_blankspace:]  # the genus specie name is left from the orginal entry
            linecount += 1
        linecount = 0
        for line in gene_fasta_data:
            gene_fasta_data[linecount] = line.strip()
            if line.startswith('ORTHOLOGY'):
                gene_fasta_data[linecount] = line.strip
                find_EC = line.find("EC:")  # finds within the line where the EC number is
                gene_fasta_data[linecount] = line[find_EC:-1].replace("[", "").replace("EC:",
                                                                                       "")  # removes the beginning bracket and EC: leaving just the number
                EC_number = "EC " + gene_fasta_data[linecount]  # adds EC back
                joined_organism_EC = [">" + str(organism_name).strip() + "%" + EC_number + "%" + gene.split(":")[
                    1]]  # adds > to beginning to find the beginning of each entry more easily, adds % between EC number and the rest of the entry to help separate the EC number for later and removes semicolon from the gene entry and adds the gene number
                DNA_info_list.append(joined_organism_EC)  # adds the entry to the blank list
            linecount += 1
        linecount = 0
        for line in gene_fasta_data:
            gene_fasta_data[linecount] = line.strip()
            if line.startswith('NTSEQ'):
                gene_fasta_data[linecount] = line.replace("NTSEQ", "").strip()
                NTSEQlocator = linecount
            linecount += 1
        DNA_data_list = gene_fasta_data[NTSEQlocator:]
        DNA_Seq = DNA_data_list[1:len(DNA_data_list) - 2]  # Takes just the DNA sequence
        separator = ""
        joined_DNA_seq = [separator.join(
            DNA_Seq)]  # combines the separate DNA sequence lines into one string and turns that into a single entry list
        DNA_info_list.append(joined_DNA_seq)  # adds single entry list to the list of lists
    return DNA_info_list


saveListsToFile(getMasterFasta(geneListFromMaster), "Master_FASTA.csv", fastaDirPath)

masterfasta = getMasterFasta(geneListFromMaster)  # should actually go above the first time it gets called

'''make a fasta file for each EC number saved in fastafoldername'''
ECorderlist = []
fasta_byEC = []
iCount = 0
for i in masterfasta:
    # print i
    if i[0].startswith(">"):
        isplit = i[0].split("%")  # finds the EC numbers using the % added previously
        # print isplit
        EC = isplit[1]
        EC_tf = "false"
        jCount = 0
        for j in ECorderlist:  # if the EC number is already present sets it to true and continues
            if EC == j:
                ECcount = jCount
                EC_tf = "true"
            jCount += 1
        if EC_tf == "false":  # if the EC number is not there, it will be added to the list
            ECorderlist.append(EC)
            fasta_byEC.append(masterfasta[iCount:iCount + 2])
        else:
            fasta_byEC[ECcount].extend(masterfasta[iCount:iCount + 2])

    iCount += 1

# print  ECorderlist
# print fasta_byEC

'''creates the FASTA files by EC numbers'''
iCount = 0
for i in ECorderlist:
    name = i.replace(".", "p").replace("EC ", "").replace(" ", "")
    saveListsToFile(fasta_byEC[iCount], name + ".csv", fastaDirPath)
    print(name + ".csv" + " saved in: " + fastaDirPath)
    iCount += 1

'''Creates ReadMe file'''
with open(newDirName + "/ReadMe.txt", "w") as ReadMe:
    ReadMe.write("KEGG_v1p1.py\n")
    ReadMe.write(now.strftime("%m-%d-%Y") + "\n")
    ReadMe.write(newDirName + "\n")
    ReadMe.write("This script creates a series of files related to the genes associated with " +
                 "plant flavonoids from various species of plants. This script first creates " +
                 "the MasterCount and MasterList files; the MasterCount counts the number genes " +
                 "each plant species have that correspond with each EC number; while the " +
                 "MasterList lists every gene with number for each plant specie. " +
                 "These are located in " + os.getcwd() + ". The script also creates files " +
                 "that only contains the genes of a single plant species biochemical pathway which " +
                 "are located in " + geneDirPath + ". The script also creates a Master FASTA files " +
                 "which contains the DNA sequence of each gene and FASTA files organized by EC " +
                 "number, these are located in " + fastaDirPath)
    ReadMe.close()

'''Gets plants that make epicatechin'''


def ECandor(listname):
    icount = 0
    codestring = ''
    for i in listname:
        if icount > 0:
            codestring += '"' + i + '"'
            try:
                listname[icount + 1]
            except IndexError:
                break
            codestring += ' ' + listname[0] + ' '
        icount += 1
    # print codestring
    return codestring


masterEC_list = [["species", "EC#s"]]
iCount = 0
for i in masterCountList:  # for each species
    # print "species :", i[0]
    speciesEClist = []
    if iCount == 0:
        pass

    else:
        jCount = 0
        for j in i:  # for EC in species
            if jCount == 0:
                speciesEClist.append(j)
                print(j)
            else:
                if str(j) == "0":
                    print(" " + str(j) + " = no enzyme")
                    pass
                else:
                    print(masterCountList[0][jCount] + " = " + str(j))
                    speciesEClist.append(masterCountList[0][jCount])
            jCount += 1
    masterEC_list.append(speciesEClist)
    iCount += 1

print(masterEC_list)
masterEC_list = masterEC_list[1:]

apigenToLute = ["or", "EC:1.14.14.81", "EC:1.14.14.82"]  # apigenin -> luteolin
cCoaToEriod = ["and", "EC:2.3.1.74"]  # caffeoyl coa -> eriodoctyol
cinnToPCoa = ["and", "EC:6.2.1.12", "EC:1.14.14.91"]  # cinnamic acid -> p coumaroyll coa
cyanToEpicat = ["and", "EC:1.3.1.77"]  # cyanidin -> epicatechin
eriodToLeuco = ["and", "EC:1.14.11.9", "EC:1.1.12.19"]  # eriodictyol -> leucocyanidin
eriodToLute = ["or", "EC:1.14.20.5", "EC:1.14.19.76"]  # eriodictyol -> luteolin
leucoToCate = ["and", "EC:1.17.1.3"]  # luteolin -> catechin
leucoToCyan = ["and", "EC:1.14.20.4"]  # luteolin -> cyanidin
narinToApigen = ["or", "EC:1.14.20.5", "EC:1.14.19.76"]  # naringenin -> apigenin
narinToEriod = ["or", "EC:1.14.14.81", "EC:1.14.14.82"]  # naringenin -> eriodictyol
pcoaToCcoa1 = ["and", "EC:1.14.13.-"]  # p coumaroyll coa -> caffeoyl coa
pcoaToCcoa2 = ["and", "EC:2.3.1.133", "EC:1.14.14.96"]  # p coumaroyll coa -> caffeoyl coa
pcoaToNarin = ["and", "EC:2.3.1.74", "EC:5.5.1.6"]  # p coumaroyll coa -> naringenin
phenylToCinn = ["or", "EC:4.3.1.24", "EC:4.3.1.25"]  # phenylalanine -> cinnamic acid

epicatList = []
cateList = []
eriodList = []
narinList = []
luteList = []
# Be careful in making these of parentheses
for i in masterEC_list:
    # print i

    epicatechin = "if ((" + ECandor(phenylToCinn) + ") and " + ECandor(cinnToPCoa) + " and (" + ECandor(
        pcoaToCcoa1) + " or (" + ECandor(
        pcoaToCcoa2) + "))" + " and " + ECandor(cCoaToEriod) + " and " + ECandor(eriodToLeuco) + " and " + ECandor(
        leucoToCyan) + " and " + ECandor(cyanToEpicat) + ") in i: epicatechinlist.append([i[0]])"
    exec(epicatechin)

    catechin = "if ((" + ECandor(phenylToCinn) + ") and " + ECandor(cinnToPCoa) + " and (" + ECandor(
        pcoaToCcoa1) + " or (" + ECandor(
        pcoaToCcoa2) + "))" + " and " + ECandor(cCoaToEriod) + " and " + ECandor(eriodToLeuco) + " and " + ECandor(
        leucoToCate) + ") in i: catechinlist.append([i[0]])"
    exec(catechin)

    eriodictyol = "if ((" + ECandor(phenylToCinn) + ") and " + ECandor(cinnToPCoa) + " and (" + ECandor(
        pcoaToCcoa1) + " or (" + ECandor(
        pcoaToCcoa2) + "))" + " and " + ECandor(cCoaToEriod) + ") in i: eriodictyollist.append([i[0]])"
    exec(eriodictyol)

    luteolin = "if ((" + ECandor(phenylToCinn) + ") and " + ECandor(cinnToPCoa) + " and (" + ECandor(
        pcoaToCcoa1) + " or (" + ECandor(
        pcoaToCcoa2) + "))" + " and " + ECandor(cCoaToEriod) + " and (" + ECandor(
        eriodToLute) + ")) in i: luteolinlist.append([i[0]])"
    exec(luteolin)
    '''luteolin= "if (("+ECandor(phenylToCinn) + ") and " + ECandor(cinnToPCoa)+ " and (" + ECandor(pcoaToCcoa1)+" or ("+ECandor(
        pcoaToCcoa2)+"))"+" and "+ECandor(cCoaToEriod)+" and "+ECandor(eriodToLute)+") in i: luteolinlist.append([i[0]])"
    exec luteolin'''  # didn't work

    naringenin = "if ((" + ECandor(phenylToCinn) + ") and " + ECandor(cinnToPCoa) + " and " + ECandor(
        pcoaToNarin) + ") in i: naringeninlist.append([i[0]])"
    exec(naringenin)

saveListsToFile(epicatList, "epicatechinspecies.txt", outDir=newDirName + '/Chemical_Data')
saveListsToFile(cateList, "catechinspecies.txt", outDir=newDirName + '/Chemical_Data')
saveListsToFile(eriodList, "eriodictyolspecies.txt", outDir=newDirName + '/Chemical_Data')
saveListsToFile(luteList, "luteolinspecies.txt", outDir=newDirName + '/Chemical_Data')
saveListsToFile(narinList, "naringeninspecies.txt", outDir=newDirName + '/Chemical_Data')
