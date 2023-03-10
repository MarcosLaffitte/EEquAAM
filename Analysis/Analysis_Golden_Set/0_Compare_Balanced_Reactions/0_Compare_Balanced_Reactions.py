################################################################################
#                                                                              #
#  README - Program: 0_Compare_Balanced_Reactions.py                           #
#                                                                              #
#  - Paper: https://match.pmf.kg.ac.rs/issues/m90n1/m90n1_75-102.html          #
#                                                                              #
#  - Github repository: https://github.com/MarcosLaffitte/EEquAAM              #
#                                                                              #
#  - Contributor(s):                                                           #
#  * @MarcosLaffitte - Marcos E. Gonzalez Laffitte                             #
#                                                                              # 
#  - Description: obtains balanced reactions from a list of reaction smiles    #
#                                                                              #
#  - Input: plain-text file with *.smiles extension whose lines are each a     #
#    single reaction-SMILES (mapped or unmapped).                              #
#                                                                              #
#  - Output: (1) plain-text *_balanced.smiles file containing balanced         #
#    reactions without any mapping. i.e., not-annotated SMILLES.               #
#                                                                              #
#  - Run with:                                                                 #
#          python  0_Compare_Balanced_Reactions.py  [myFile.smiles]            #
#                                                                              #
#  - Expected output:                                                          #
#    (1)   myFile_balanced.smiles                                              #
#                                                                              #
#  - Date: 12 January 2023                                                     #
#                                                                              #
################################################################################


# Requires / Used Versions #####################################################
"""
> Language: python 3.9.13
> Anaconda: conda 22.9.0
> Required packages could be manually installed preferably in the following order
> Packages installed with anaconda's pip:
***** networkx 2.8.4
***** pysmiles 1.0.2
***** chytorch-rxnmap 1.3 
      (https://github.com/chython/chytorch-rxnmap)
***** chython 1.50 
      (with chython library DEV version as in https://chython.readthedocs.io/en/latest/)
***** chytorch 1.27 
      (required by chytorch-rxnmap; https://github.com/chython/chytorch)
***** transformers 4.24.0 
      (mind the version by running: pip install transformers==4.24.0)
"""


# Dependencies #################################################################


# installed with conda's pip ---------------------------------------------------
import networkx as nx
import pysmiles as ps
from chython import smiles


# already in python ------------------------------------------------------------
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
outputFileName = inputFileName.replace(".smiles", "_balanced.smiles")
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


# function: determin if a given reaction SMILES is balanced --------------------
def isBalanced(someSMILES):
    # local variables
    mol = None
    reactants = []
    products = []
    typesR = dict()
    typesP = dict()
    chemElement = dict()
    balancedAns = True        
    # only carbon and heteroatoms are considered for balance-check
    # aam of e.g. protons [H+] depends on input SMILES and mappers
    ignoreElements = ["H"]
    # get number of non-hydrogen atoms in reactants    
    reactants = list((someSMILES.split(">>")[0]).split("."))
    for eachReactant in reactants:
        # read smiles
        mol = ps.read_smiles(eachReactant)
        # obtain atom type and number per type
        chemElement = nx.get_node_attributes(mol, "element")
        # the only atoms without element are "*" unknowns
        for atom in list(chemElement.keys()):
            if(not(chemElement[atom] in ignoreElements)):
                if(chemElement[atom] in list(typesR.keys())):
                    typesR[chemElement[atom]] = typesR[chemElement[atom]] + 1
                else:
                    typesR[chemElement[atom]] = 1
    # get number of non-hydrogen atoms in products
    products = list((someSMILES.split(">>")[1]).split("."))
    for eachProduct in products:
        # read smiles
        mol = ps.read_smiles(eachProduct)
        # obtain atom type and number per type
        chemElement = nx.get_node_attributes(mol, "element")
        # the only atoms without element are "*" unknowns
        for atom in list(chemElement.keys()):
            if(not(chemElement[atom] in ignoreElements)):
                if(chemElement[atom] in list(typesP.keys())):
                    typesP[chemElement[atom]] = typesP[chemElement[atom]] + 1
                else:
                    typesP[chemElement[atom]] = 1
    # evaluate consistency of nonH atoms
    intersecTypes = list(set(typesR.keys()).intersection(set(typesP.keys())))
    # case 1: new element types appear, disappear or change, e.g. CN>>CNO or CNS>>CNP
    if(not (len(intersecTypes) == len(list(typesR.keys())) and len(intersecTypes) == len(list(typesP.keys())))):
        balancedAns = False        
    else:
        # case 2: number of elements of each type changes        
        for eachType in intersecTypes:
            if(not typesR[eachType] == typesP[eachType]):
                balancedAns = False
                break    
    # end of function
    return(balancedAns)
    

# Main #########################################################################


# initial message
print("\n")
print(">>> 0_Compare_Balanced_Reactions")


# get reactions from input
inputFile = open(inputFileName, "r")
inputSMILES = inputFile.read().splitlines()
inputFile.close()


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
print("* filtering out unbalanced reactions ...")


# filtering out unbalanced reactions
for i in range(len(inputSMILES)):
    balanced = isBalanced(inputSMILES[i])
    if(not balanced):
        unsuitableReactions.append((inputSMILES[i]))

        
# remove unsuitable reactions
for (badInput) in unsuitableReactions:
    inputSMILES.remove(badInput)

        
# task message
print("\n")
print("* printing SMILES of " + str(len(inputSMILES)) + " balanced reactions ...")
    

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
