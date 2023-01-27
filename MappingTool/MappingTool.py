################################################################################
#                                                                              #
#  README - Program: MappingTool.py                                            #
#                                                                              #
#  - Paper: [tba]                                                              #
#                                                                              #
#  - Github repository: https://github.com/MarcosLaffitte/EEquAAM              #
#                                                                              #
#  - Date: 12 January 2023                                                     #
#                                                                              #
#  - Contributor(s):                                                           #
#    * @MarcosLaffitte - Marcos E. Gonzalez Laffitte                           #
#                                                                              #
#  - Description: this script receives a set of reaction-SMILES (mapped or     #
#    unmapped), removes any map from them through the string representation of #
#    the objects build with GraphormerMapper (see below), and evaluates the    #
#    following 3 atom-mapping tools over them (last checked 12 January 2023):  #
#                                                                              #
#    * RXNMapper 0.2.4 [RXN]                                                   #
#      https://pypi.org/project/rxnmapper/                                     #
#    * Reaction Decoder Tool version 2.4.1  [RDT]                              #
#      requires java; is the *.jar file in same folder                         #
#      https://github.com/asad/ReactionDecoder/releases                        #
#    * GraphormerMapper [CHY]                                                  #
#      also pip package python chytorch-rxnmap 1.3                             #
#      https://github.com/chython/chytorch-rxnmap                              #
#      https://pypi.org/project/chytorch-rxnmap/                               #
#                                                                              #
#    Finally it returns a plain text file containing the mapped-SMILES of      #
#    only those reactions that where originally balanced and that where        #
#    completely mapped by the three mappers, together with id's (RXNmap,       #
#    RDTmap, CHYmap) indicating which atom-mapper tool obtained which          #
#    mapped-SMILES string. This script calls the RDT application from the      #
#    command line with the os package, since that mapper is written in java    #
#    and thus RDT_2.4.1.jar should always be in same dir as this script.       #
#                                                                              #
#  - Input: plain-text file with *.smiles extension whose lines are each a     #
#    single reaction-SMILES. For general specifications of both anotated and   #
#    unmapped SMILES strings see: http://opensmiles.org/opensmiles.html        #
#                                                                              #
#  - Output: (1) plain-text (csv) *_aam.smiles file containing only those      #
#    reaction-SMILES that where balanced and which where mapped by the three   #
#    tools, followed by the 3 resulting maps including the ID of the           #
#    corresponding mapper. This information is ordered in columns separated    #
#    by commas as show in the following toy example:                           #
#                                                                              #
#    # , reaction_1:CCCCCCCCCC>>CCCCCCCCC.C                                    #
#    RXNmap , [CH3:1][...][CH:10]>>[CH3:1][...][CH3:9].[CH4:10]                #
#    RDTmap , [CH3:2][...][CH:4]>>[CH3:1][...][CH3:9].[CH4:2]                  #
#    CHYmap , [CH3:6][...][CH:7]>>[CH3:7][...][CH3:1].[CH4:6]                  #
#    # , reaction_2:CCCCCCCCCC>>CCCC.CCCCCC                                    #
#    RXNmap , [CH3:1][...][CH:10]>>[CH3:1][...][CH3:4].[CH3:5][...][CH4:10]    #
#    RDTmap , [CH3:2][...][CH:4]>>[CH3:1][...][CH3:9].[CH4:2][...][CH4:8]      #
#    CHYmap , [CH3:6][...][CH:7]>>[CH3:7][...][CH3:1].[CH3:8][...][CH3:4]      #
#    ...                                                                       #
#                                                                              #
#    * the line with the original SMILES will always start with "#"            #
#    * there is no limit to the number of reactions that can be analized       #
#    * "reaction_i:" is an arbitrary identifier and can be changed in-code     #
#    * blank spaces between commas are only for presentation purposes          #
#      and should not be present in the output file                            #
#                                                                              #
#  - Run with (after activating eequaam conda environment):                    #
#                 python  MappingTool.py  [myFile.smiles]                      #
#                                                                              #
#  - Expected output:                                                          #
#    (1)   myFile_aam.smiles                                                   #
#                                                                              #
#  - Notes:                                                                    #
#    * if the program found a "bad" SMILES it means this is a SMILES string    #
#      with either empty reactants, empty products, or containing the symbol   #
#      "~" which is not currently supported by some of the implmented mappers. #
#                                                                              #
#    * please note that the output file is a plain-text file and so you have   #
#      to zoom it out to PROPERLY see the actual output lines :)               #
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
inputFile = None


# output -----------------------------------------------------------------------
outputFileName = inputFileName.replace(".smiles", "_aam.smiles")
outputFile = None
resultsByReaction = []
resultsFormat = ""
fileContent = ""


# data -------------------------------------------------------------------------
reactantSide = ""
productSide = ""
badSMILES = 0
nonH = 0
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
print(">>> MappingTool - EEquAAM Github Repository")


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
    # print progress
    printProgress(round((i+1)*100/len(inputSMILES), 2), i+1, len(inputSMILES))
inputSMILES = deepcopy(cleanSMILES)


# task message
print("\n")
print("* found " + str(badSMILES) + " bad reaction-SMILES ...")


# task message
print("\n")
print("* applying RXN mapper ...")


# call RXNMapper (for each individual reaction so as to reduce memory consumption)
for i in range(len(inputSMILES)):
    # reinitialize RNXMapper
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
print("* determining balanced reactions with complete maps ...")


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
    # reactions should be balanced
    nonH = 2*nonH
    # get consistency of RXN map based on bracket-atoms
    totBrackets = len([char for char in mapsRXN[i] if char == "["])
    if(totBrackets == nonH):
        consistencyRXN = True
    # get consistency of RDT map based on bracket-atoms
    totBrackets = len([char for char in mapsRDT[i] if char == "["])
    if(totBrackets == nonH):
        consistencyRDT = True
    # get consistency of CHY map based on bracket-atoms
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
print("* printing SMILES of " + str(len(inputSMILES)) + " balanced and completely mapped reactions ...")
    

# print final results
resultsByReaction = []
for i in range(len(inputSMILES)):
    originalReactionID = "reaction_" + str(i+1) + ":"
    resultsFormat = "#," + originalReactionID  + inputSMILES[i] + "\n"
    resultsFormat = resultsFormat + "RXNmap" + "," + mapsRXN[i] + "\n"
    resultsFormat = resultsFormat + "RDTmap" + "," + mapsRDT[i] + "\n"
    resultsFormat = resultsFormat + "CHYmap" + "," + mapsCHY[i] + "\n"
    resultsByReaction.append(resultsFormat)
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
