################################################################################
#                                                                              #
#  README - Program: 1_Obtain_Balanced_and_Suitable_Reactions.py               #
#                                                                              #
#  - Paper: [tbd]                                                              #
#                                                                              #
#  - Github: [tbd]                                                             #
#                                                                              #
#  - Contributor(s):                                                           #
#  * @MarcosLaffitte - Marcos E. Gonzalez Laffitte                             #
#                                                                              # 
#  - Description: This is a version of the script                              #
#                                                                              #
#                3_Maps_over_suitable_from_Golden_Set.py,                      #
#                                                                              #
#    but made specifically to filter-out unbalanced and unsuitable reactions   #
#    from a list of arbitrary reaction-SMILES. Here "suitable" refers to those #
#    (balanced) reactions that at the same time CAN be completely mapped by    #
#    the 3 atom-mapping tools: RXNMapper, RDT mapper and GraphormerMapper,     #
#    thus allowing us to obtain a (bijective) aam to carry on with our study.  #
#    For information on these mappers see the  3_Maps_over_...py  script.      #
#                                                                              #
#  - Input: plain-text file with *.smiles extension whose lines are each a     #
#    single reaction-SMILES (mapped or unmapped).                              #
#                                                                              #
#  - Output: (1) plain-text *_suitable.smiles file containing balanced and     #
#    suitable reactions, without any mapping. i.e., not-annotated SMILLES.     #
#                                                                              #
#  - Run with:                                                                 #
#      python  1_Obtain_Balanced_and_Suitable_Reactions.py  [myFile.smiles]    #
#                                                                              #
#  - Expected output:                                                          #
#    (1)   myFile_suitable.smiles                                              #
#                                                                              #
#  - Date: 12 January 2023                                                     #
#                                                                              #
#  ***** NOTE: found BAD SMILES in original Golden Data Set *****              #
#        see file in line:                                                     #
#        * 367 - empty products                                                #
#        * 1121, 1139, 1849, 1850 - contain "~" symbol                         #
#  Supporting information: https://pubs.acs.org/doi/10.1021/acs.jcim.2c00344   #
#                                                                              #
################################################################################


# Requires / Used Versions #####################################################
"""
> Language: python 3.9.13
> Anaconda: conda 22.9.0
> Required packages could be manually installed preferably in the following order
> Packages installed with anaconda's pip:
***** pysmiles 1.0.2
***** rxnmapper 0.2.4
***** rdkit 2022.9.3 (required by rxnmapper)
***** chytorch-rxnmap 1.3 
      (https://github.com/chython/chytorch-rxnmap)
***** chython 1.50 
      (with chython library DEV version as in https://chython.readthedocs.io/en/latest/)
***** chytorch 1.27 
      (required by chytorch-rxnmap; https://github.com/chython/chytorch)
***** transformers 4.24.0 
      (mind the version by running: pip install transformers==4.24.0)
> Java (OpenJDK) and Ubuntu [java --version]
***** openjdk 11.0.17 2022-10-18
***** OpenJDK Runtime Environment (build 11.0.17+8-post-Ubuntu-1ubuntu222.04)
***** OpenJDK 64-Bit Server VM (build 11.0.17+8-post-Ubuntu-1ubuntu222.04, mixed mode)
"""


# Dependencies #################################################################


# installed with conda's pip ---------------------------------------------------
import pysmiles as ps
from chython import smiles
from rxnmapper import RXNMapper


# already in python ------------------------------------------------------------
import os
from sys import argv
from copy import deepcopy
from math import factorial, modf


# only report erros and not deprecation warnings -------------------------------
import warnings
from transformers import logging
warnings.filterwarnings("ignore")
logging.set_verbosity_error()


# variables ####################################################################


# input ------------------------------------------------------------------------
inputFileName = argv[1]
inputSMILES = []
cleanSMILES = []
originalSMILES = dict()
inputFile = None


# output -----------------------------------------------------------------------
outputFileName = inputFileName.replace(".smiles", "_suitable.smiles")
outputFile = None
resultsByReaction = []
resultsFormat = ""
fileContent = ""


# data -------------------------------------------------------------------------
reactantSide = ""
productSide = ""
nonH = 0
badSMILES = 0
totBrackets = 0
RXNaam = ""
RDTaam = ""
CHYaam = ""
mapsRXN = []
mapsRDT = []
mapsCHY = []
consistencyRXN = False
consistencyRDT = False
consistencyCHY = False
myrxn = None
myrxnUnmapped = ""
rxn_mapper = None
unsuitableReactions = []
badInput = ""
badRXN = ""
badRDT = ""
badCHY = ""


# Functions ####################################################################


# function: print custom progress bar ------------------------------------------
def printProgress(casePercentage, caseNum, totCases):
    # local variables
    tail = "".join(10*[" "])
    base = "-"
    done = "="
    bar = ""        
    pile = []
    finished = ""
    percentageInt = 0
    # generate bar
    percentageInt = int(modf(casePercentage/10)[1])
    for i in range(1, 11):
        if(i <= percentageInt):
            pile.append(done)
        else:
            pile.append(base)
    finished = "".join(pile)
    bar = "- progress:   0%  [" + finished + "]  100%" + " ;  done: " + str(caseNum) + " / " + str(totCases)
    # message    
    print(bar + tail, end = "\r")


# Main #########################################################################


# initial message
print("\n")
print(">>> 1_Obtain_Balanced_and_Suitable_Reactions")


# get reactions from input
inputFile = open(inputFileName, "r")
inputSMILES = inputFile.read().splitlines()
inputFile.close()


# print BAD SMILES Golden Set
# print(inputSMILES[366])
# print(inputSMILES[1120])
# print(inputSMILES[1138])
# print(inputSMILES[1848])
# print(inputSMILES[1849])


# task message
print("\n")
print("* received " + str(len(inputSMILES)) + " reaction-SMILES ...")


# task message
print("\n")
print("* removing any previous maps in the given SMILES ...")


# remove any previous maps using the __str__ provided by GraphormerMapper objects
badSMILES = 0
cleanSMILES = []
for i in range(len(inputSMILES)):
    # check non-empty strings of reactants and products and no anomalous symbols
    reactantSide = inputSMILES[i].split(">>")[0]
    productSide = inputSMILES[i].split(">>")[1]
    # save data and print line in file of bad smiles
    if((reactantSide == "") or (productSide == "") or ("~" in inputSMILES[i])):
        # print(i+1, inputSMILES[i])
        badSMILES = badSMILES + 1
        continue
    else:
        myrxn = smiles(inputSMILES[i])    
        myrxnUnmapped = myrxn.__str__()
        cleanSMILES.append(myrxnUnmapped)
        originalSMILES[myrxnUnmapped] = inputSMILES[i]
    # print progress
    printProgress(round((i+1)*100/len(inputSMILES), 2), i+1, len(inputSMILES))
inputSMILES = deepcopy(cleanSMILES)


# task message
print("\n")
print("* found " + str(badSMILES) + " BAD reaction-SMILES ...")


# task message
print("\n")
print("* applying RXN mapper ...")


# call RXNMapper (for each individual reaction so as to reduce memory consumption)
for i in range(len(inputSMILES)):
    # reinitialize RXNMapper
    rxn_mapper = RXNMapper()
    # map reaction
    RXNaam = rxn_mapper.get_attention_guided_atom_maps([inputSMILES[i]])[0]["mapped_rxn"]
    mapsRXN.append(RXNaam)
    # print progress
    printProgress(round((i+1)*100/len(inputSMILES), 2), i+1, len(inputSMILES))    

    
# task message
print("\n")
print("* applying RDT mapper ...")
    

# call RDT (in silet mode)
for i in range(len(inputSMILES)):
    os.system("nohup java -jar RDT_2.4.1.jar -Q SMI -q " + "\"" + inputSMILES[i] + "\"" + " -c -j AAM -f TEXT >/dev/null 2>&1")
    inputFile = open("ECBLAST_smiles_AAM.txt", "r")
    RDTaam = (inputFile.read().splitlines())[3]
    inputFile.close()
    os.system("rm *.txt")
    os.system("rm *.rxn")
    mapsRDT.append(RDTaam)
    printProgress(round((i+1)*100/len(inputSMILES), 2), i+1, len(inputSMILES))


# task message
print("\n")
print("* applying Graphormer mapper ...")


# call Graphormer
for i in range(len(inputSMILES)):
    myrxn = smiles(inputSMILES[i])
    myrxn.reset_mapping()
    CHYaam = format(myrxn, "m")
    # get result
    mapsCHY.append(CHYaam)
    printProgress(round((i+1)*100/len(inputSMILES), 2), i+1, len(inputSMILES))
   
 
# task message
print("\n")
print("* checking all atoms in each reaction were mapped by all mappers and")
print("  that the given reactions where in fact initially balanced ...")


# checking all atoms were mapped by all mappers
for i in range(len(inputSMILES)):
    # reinitizlize consistency
    consistencyRXN = False
    consistencyRDT = False
    consistencyCHY = False
    mol = None
    nonH = 0
    reactants = []
    # get number of non-hydrogen atoms involved in reaction
    reactants = list((inputSMILES[i].split(">>")[0]).split("."))
    for eachReactant in reactants:
        mol = ps.read_smiles(eachReactant)
        nonH = nonH + len(list(mol.nodes()))
    # reactions should be initially balanced
    nonH = 2*nonH
    # get consistency of RXN map based on number of bracket-atoms
    totBrackets = len([char for char in mapsRXN[i] if char == "["])
    if(totBrackets == nonH):
        consistencyRXN = True
    # get consistency of RDT map based on number of bracket-atoms
    totBrackets = len([char for char in mapsRDT[i] if char == "["])
    if(totBrackets == nonH):
        consistencyRDT = True
    # get consistency of CHY map based on number of bracket-atoms
    totBrackets = len([char for char in mapsCHY[i] if char == "["])
    if(totBrackets == nonH):
        consistencyCHY = True
    # evaluate consistency of all mappers
    if(not(consistencyRXN and consistencyRDT and consistencyCHY)):
        unsuitableReactions.append((inputSMILES[i], mapsRXN[i], mapsRDT[i], mapsCHY[i]))

        
# remove unsuitable reactions
for (badInput, badRXN, badRDT, badCHY) in unsuitableReactions:
    inputSMILES.remove(badInput)
    mapsRXN.remove(badRXN)
    mapsRDT.remove(badRDT)
    mapsCHY.remove(badCHY)

        
# task message
print("\n")
print("* printing SMILES of " + str(len(inputSMILES)) + " balanced and suitable (completely mapped) reactions ...")
    

# print final results
resultsByReaction = []
for i in range(len(inputSMILES)):
    resultsByReaction.append(originalSMILES[inputSMILES[i]] + "\n")
fileContent = "".join(resultsByReaction)    
outputFile = open(outputFileName, "w")
outputFile.writelines(fileContent)
outputFile.close()


# final message
print("\n")
print(">>> Finished")
print("\n")


# End ##########################################################################
################################################################################
