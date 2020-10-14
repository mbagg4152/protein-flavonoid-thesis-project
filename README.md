# Chem KEGG Py-Project
Former repository for the research methods class. The goal was to modernize, clean and adapt a former student's thesis project.  
Currently this is repository is holding code used in research.

---
## Original Code Information
Originally Created on Wed Apr 17 2019. @author: vmoorman and Jordan Wilson  
KEGG.py, Made using Python 2.7  
JW started the code to get information from KEGG about the species we were interested in - April 2019  

---
## Project Structure
> chem-py-project - main folder  
chem-py-project/original/ - contains the code from when Jordan Wilson was developing  
chem-py-project/project/ - current code  
chem-py-project/project/data/ - output data folder  
chem-py-project/project/data/Chemical_Data - contains the files for the list of plants predicted per flavonoid  
chem-py-project/project/data/FASTA_Data - contains FASTA EC number data  
chem-py-project/project/data/Gene_Data - contains the data pulled from KEGG for each plant pathway  
chem-py-project/project/json-data -  holds the lists of plants and pathways used in the KEGG program (in JSON format).  
chem-py-project/project/lib - contains the library/helper files  
chem-py-project/project/test-output - contains output from testing different KEGG functions in ktester.py  
  
---
# Important Files and Functions
__Note:__ Not all files or functions are covered below.  


## ```chem-py-project/project/kegg-prog.py```
This is the main program of the code project.  

#### ```gene_pathway_data```
This function is called and passed a pathway for a specific organism (e.g. adu00941). For each pathway that is passed in, the code uses the kegg plugin in order to get the gene data for the specific pathway, which is done using the following lines:  
```python
    raw = kegg.get(pathway_id)
    gd = kegg.parse(raw)
    ...
    fetched_genes = gd.get('GENE')
```
Once the gene data is collected,the data is split up into different list items. For reference, a single gene entry for cam00941 is shown below (the other entries in 'GENE' take similar form):  
```
'GENE': {'101489106': 'chalcone synthase 1 [KO:K00660] [EC:2.3.1.74]', ...}
```
Then the data is added to a list such that this gene data is associated with the appropriate plant.


---
## ```chem-py-project/project/lib/compoundinfo.py```
This contains the labeled EC numbers as well as the logic used in order to make the predictions.  
The EC numbers are labeled as such:
```python
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
E17_2 = 'EC:1.1.1.219 1.1.1.234' # unique number for bifunctional dihydroflavonol 4-reductase/flavanone 4-reductase
E18 = 'EC:1.14.20.4'
E19 = 'EC:1.3.1.77'
E20 = 'EC:1.17.1.3'
E21 = 'EC:1.14.14.87'
E22 = 'EC:4.2.1.105'
E23 = 'EC:2.4.1.74'
E24 = 'EC:2.3.1.70'
E25 = 'EC:2.3.1.30'
E26 = 'EC:2.4.1.136'
```
### ```or_in```
Takes in a list and elements and returns true if at least ONE element is present in the list.

### ```and_in```
Takes in a list and elements and returns true only if all of the passed elements are in the list.

### The Logical Functions
Different logical functions have been written not only for the flavonoids of interest, but also for the prerequisite compounds which are found on the map. The prerequisite functions are used to get the total result for the specific compound. If the prerequisite returns False, the compound logic function will also return false.  
Each function is abbreviated in order to keep code lines at a decent length, but each function does have its compound's full name commented at the end of its respective line. For example, the function for Cinnamic acid, which requires EC:4.3.1.24 OR EC:4.3.1.25:

```
def cia(e): return or_in(e, E1, E2)  # cinnamic acid
```
The logical functions are used in the function ```finish_up``` in kegg-prog.py. For each plant's total EC list, the program will loop through each of the flavonoids' logical requirements, which are held held in the list of ```ChemData``` objects called ```data_lists```, where each item ```chem_data``` has a property ```chem_data.label``` that is passed to a function in ```compoundinfo.py``` called ```flav_check```, which then determines the logical function to be called. If the result of ```flav_check``` returns as true, then the current plant's name will be appended to the list of plants, which is held in ```chem_data.species```.   

---
## ```chem-py-project/project/lib/knapsackinfo.py```
This program uses the 'wget' command in order to pull the information for each of the species from KNApSAcK. For each plant, it compiles a list of entries from the database if the entry line contains one of the flavonoids of interest. This program was written with the purpose of making data collection for both the compounds and their relatives easier.  
More searching will be done in order to determine how easy or difficult it would be to update this code to work with other databases.

__Note:__ It will not work correctly if the wget system command is not installed (currently unsure if wget is available for Windows). 



---
## ```chem-py-project/project/lib/datatypes.py```
This contains the custom data types that are or have been used in the program.  
  

---
## ``` chem-py-project/project/lib/jsondata.py```
This file calls the ```get_json_data(filename,key)``` function from util.py, which reads in the list of plant and pathway codes as well as the file containing the scientific name for each plant and the full name of each pathway map.




---
## ```chem-py-project/project/lib/pathstrings.py```
This file just contains the strings which hold the dedicated output folder and file names.



---
## ```chem-py-project/project/lib/util.py```
This file contains various utility functions used throughout the program.



#### ```get_json_data```  
Reads in a JSON (JavaScript Object Notation) file and converts the data into usable python variables.




#### ```remove_dupes```   
Removes duplicate elements from a list.



#### ```write_readme```  
Writes the program's original README file.


--
## Original Changelog 
Version | Change |
:-------|:-------|
<sub>v0p1 | <sub>VRM getting the iteration and parsing part of the code to work - April 2019 |
<sub>v0p2 | <sub>VRM making the count code work - April 2019 |
<sub>v0p3 | <sub>JW getting the file output to work, revising the species codes and adding a dictionary - April 2019 |
<sub>v0p4 | <sub>JW creating directory change function for Gene_data and adding additional plant species - multiple things still broken - May 2019 |
<sub>v0p5 | <sub>VRM cleaned up some duplicate issues, folder locations, and ensured that species are written out in full  - May 2019 |
<sub>v0p6 | <sub>JW creating function for master FASTA file (error occurs when running every entry) - May 2019 |
<sub>v0p7 | <sub>JW editing lists for master fasta function and updating fasta function - May 2019 |
<sub>v0p8 | <sub>VRM fixing and creating gene list for master fasta file - June 2019 |
<sub>v0p8p2 | <sub>JW made a fasta file for each EC#- June 2019 |
<sub>v1p0 | <sub>VRM cleaned up code |
<sub>v1p1p1 | <sub>JW added ReadMe |
<sub>v1p2 | <sub>JW added function to get plants that make Epicatechin - incorrectly |
<sub>v1p3 | <sub>VRM coded in epicatechin, catechin, eriodictyol, luteolin, naringenin listed in Chemical_Data --- leutolin actually wrong |
<sub>v1p3p1 | <sub>VRM tried to fix leutolin |




