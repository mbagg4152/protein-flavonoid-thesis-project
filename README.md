# Chem KEGG Py-Project
Repo for the research methods class. The goal is to modernize, clean and adapt a former student's thesis project.

## Original Code Information
Originally Created on Wed Apr 17 2019. @author: vmoorman and Jordan Wilson  
KEGG.py, Made using Python 2.7  
JW started the code to get information from KEGG about the species we were interested in - April 2019  

### Project Structure
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




