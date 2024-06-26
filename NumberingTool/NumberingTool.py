################################################################################
#                                                                              #
#  README - Program: NumberingTool.py                                          #
#                                                                              #
#  - Paper: https://match.pmf.kg.ac.rs/issues/m90n1/m90n1_75-102.html          #
#                                                                              #
#  - Github repository: https://github.com/MarcosLaffitte/EEquAAM              #
#                                                                              #
#  - Date: 12 January 2023                                                     #
#                                                                              #
#  - Contributor(s) to this script:                                            #
#    * @MarcosLaffitte - Marcos E. Gonzalez Laffitte                           #
#                                                                              # 
#  - Description: this script receives a plain-text file with "*.smiles"       #
#    extension containing a list of reaction SMILES like in the following      #
#    toy-example:                                                              #
#                                                                              #
#    CCCCCCCCCC>>CCCCCCCCC.C                                                   #
#    CCCCCCCCCC>>CCCC.CCCCCC                                                   #
#    ...                                                                       #
#                                                                              #
#    * there is no limit to the number of reactions that can provided          #
#                                                                              #
#    Then the program asigns identifiers to each reaction as follows:          #
#                                                                              #
#    #,reaction_1                                                              #
#    CCCCCCCCCC>>CCCCCCCCC.C                                                   #
#    #,reaction_2                                                              #
#    CCCCCCCCCC>>CCCC.CCCCCC                                                   #
#    ...                                                                       #
#                                                                              #
#    * lines containing the identifier always start with "#,"                  #
#    * the id "reaction_i" is an arbitrary identifier and can be changed to    #
#      be e.g. the common name of a reaction, including more commas if         #
#      needed but no newline characters "\n", i.e, this is a single line       #
#    * this format is the input for the MappingTool.py in the repo             #
#                                                                              #
#    For general specifications of both anotated and unmapped SMILES strings   #
#    see: http://opensmiles.org/opensmiles.html                                #
#                                                                              #
#  - Input: plain-text file with *.smiles extension whose lines are            #            
#    reaction SMILES strings.                                                  #
#                                                                              #
#  - Output: (1) plain-text *_id.smiles file whose lines are alternating       #
#    identifiers and reaction SMILES.                                          #
#                                                                              #
#  - Run with (after activating eequaam conda environment):                    #
#                python  NumberingTool.py  [myFile.smiles]                     #
#                                                                              #
#  - Expected output:                                                          #
#    (1)   myFile_id.smiles                                                    #
#                                                                              #
#  - Notes:                                                                    #
#                                                                              #
#    * the output of this script is the input for MappingTool.py in the repo   #
#                                                                              #
#  --------------------------------------------------------------------------  #
#                                                                              #
# - LICENSE:                                                                   #
#                                                                              #
#   This file is part of the work published in                                 #
#            https://match.pmf.kg.ac.rs/issues/m90n1/m90n1_75-102.html         #
#   and it is released under                                                   #
#            MIT License Copyright (c) 2023 Marcos E. González Laffitte        #
#   See LICENSE file in                                                        #
#            https://github.com/MarcosLaffitte/EEquAAM                         #
#   for full license details.                                                  #
#                                                                              #
################################################################################


# Dependencies #################################################################


# already in python ------------------------------------------------------------
from sys import argv, exit
from copy import deepcopy
from math import factorial, modf


# only report erros and not deprecation warnings -------------------------------
import warnings
warnings.filterwarnings("ignore")


# variables ####################################################################


# input ------------------------------------------------------------------------
inputFileName = ""
inputLines = []
inputSMILES = dict()
inputFile = None


# check input ------------------------------------------------------------------
if(".smiles" in argv[1]):
    remainder = (argv[1].split(".smiles"))[-1]
    if(not remainder == ""):
        errorStr = "\n >> NumberingTool: Wrong input extension.\n"
        errorStr = errorStr + "- Expected: *.smiles\n"
        errorStr = errorStr + "- Received: *.smiles" + remainder + "\n"
        exit(errorStr)
    else:
        inputFileName = argv[1]
else:
    exit("\n >> NumberingTool: missing *.smiles extension.\n")


# output -----------------------------------------------------------------------
outputFileName = inputFileName.replace(".smiles", "_id.smiles")
outputLines = []
fileContent = []
outputFile = None


# Main #########################################################################


# initial message
print("\n")
print(">>> NumberingTool - EEquAAM Github Repository")


# get reactions from input
inputFile = open(inputFileName, "r")
inputLines = inputFile.read().splitlines()
inputFile.close()


# here one can define the identifiers (should be the same length of given list)
myIdentifiers = []
# comment-out the following loop if defining the myIdentifiers list manually
for i in range(len(inputLines)):
    eachIdentifier = "reaction_" + str(i+1)
    myIdentifiers.append(eachIdentifier)


# task message
print("\n")
print("* assigning identifiers ...")

    
# read smiles and add identifiers
for i in range(len(inputLines)):
    inputSMILES[myIdentifiers[i]] = inputLines[i]


# task message
print("\n")
print("* saving data ...")

    
# print results
if(len(list(inputSMILES.keys())) > 0):
    fileContent = []
    outputStr = ""
    for eachReaction in list(inputSMILES.keys()):
        fileContent.append("#," + eachReaction + "\n")
        fileContent.append(inputSMILES[eachReaction] + "\n")
    outputStr = "".join(fileContent)
    outputFile = open(outputFileName, "w")
    outputFile.writelines(outputStr)
    outputFile.close()
    

# final message
print("\n")
print(">>> Finished")
print("\n")


# End ##########################################################################
################################################################################
