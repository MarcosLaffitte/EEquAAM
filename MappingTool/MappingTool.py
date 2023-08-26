################################################################################
#                                                                              #
#  README - Program: MappingTool.py                                            #
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
#  - Description: this script receives a set of UNMAPPED reaction-SMILES,      #
#    checks that they are suitable and balanced, and evaluates the following   #
#    three atom-mapping tools over them (last checked 12 January 2023):        #
#                                                                              #
#    * [RXN] RXNMapper 0.2.4                                                   #
#      https://pypi.org/project/rxnmapper/                                     #
#    * [RDT] Reaction Decoder Tool version 2.4.1                               #
#      requires java; is the *.jar file in same folder                         #
#      https://github.com/asad/ReactionDecoder/releases                        #
#    * [CHY] GraphormerMapper                                                  #
#      also pip package python chytorch-rxnmap 1.3                             #
#      https://github.com/chython/chytorch-rxnmap                              #
#      https://pypi.org/project/chytorch-rxnmap/                               #
#                                                                              #
#    Finally it returns a plain text file containing the mapped-SMILES of      #
#    only those suitable and balanced reactions that where completely mapped   #
#    by the three mappers, together with id's (RXNmap, RDTmap, CHYmap)         #
#    indicating which mapping tool obtained which (bijective) atom-to-atom     #
#    map. This program calls the RDT application from the command line with    #
#    the os package, since that mapper is written in java and thus the file    #
#    RDT_2.4.1.jar should always be in the same dir as this script.            #
#                                                                              #
#  - Input: plain-text file with *.smiles extension whose lines alternate      #
#    between identifiers and their corresponding UNMAPPED reaction SMILES.     #
#    For general specifications on both anotated and unmapped SMILES           #
#    strings see: http://opensmiles.org/opensmiles.html                        #
#    The input follows the format in the next toy-example:                     #
#                                                                              #
#    #,reaction_1                                                              #
#    CCCCCCCCCC>>CCCCCCCCC.C                                                   #
#    #,reaction_2                                                              #
#    CCCCCCCCCC>>CCCC.CCCCCC                                                   # 
#    ...                                                                       #
#                                                                              #
#  - Output: (1) plain-text (csv) *_aam.smiles file containing only those      #
#    reaction-SMILES that where suitable, balanced and which where mapped by   #
#    the three tools, followed by the 3 resulting maps including the ID of the #
#    corresponding mapper. This information is ordered in columns separated    #
#    by commas as show in the following toy-example:                           #
#                                                                              #
#    #,reaction_1                                                              #
#    #,CCCCCCCCCC>>CCCCCCCCC.C                                                 #
#    RXNmap,[CH3:1][...][CH:10]>>[CH3:1][...][CH3:9].[CH4:10]                  #
#    RDTmap,[CH3:2][...][CH:4]>>[CH3:1][...][CH3:9].[CH4:2]                    #
#    CHYmap,[CH3:6][...][CH:7]>>[CH3:7][...][CH3:1].[CH4:6]                    #
#    #,reaction_2                                                              #
#    #,CCCCCCCCCC>>CCCC.CCCCCC                                                 #
#    RXNmap,[CH3:1][...][CH:10]>>[CH3:1][...][CH3:4].[CH3:5][...][CH4:10]      #
#    RDTmap,[CH3:2][...][CH:4]>>[CH3:1][...][CH3:9].[CH4:2][...][CH4:8]        #
#    CHYmap,[CH3:6][...][CH:7]>>[CH3:7][...][CH3:1].[CH3:8][...][CH3:4]        #
#    ...                                                                       #
#                                                                              #
#    * lines with the original SMILES or identifiers always start with "#,"    #
#    * there is no limit to the number of reactions that can be analized       #
#    * the id "reaction_i" is an arbitrary identifier and can be changed to be #
#      e.g. the common name of a reaction, including more commas if needed     #
#      but not newline characters "\n", i.e, this is a single line             #
#    * WARNING! reactions with same identifiers will be silently overwritten   #
#                                                                              #
#    Also three other files are returned after running the script: (2) the     #
#    SMILES refered to as "unsuitable" which: have no reactans or no products, #
#    are missing the ">>" symbol, or were previously mapped. Then (3) includes #
#    the unbalanced SMILES, and finally (4) contains those reactions that      #
#    where not completeley mapped by the three atom-mapping tools. This last   #
#    file has the same format as (1) but displaying the uneven maps.           #
#                                                                              #
#  - Run with (after activating eequaam conda environment):                    #
#                 python  MappingTool.py  [myFile.smiles]                      #
#                                                                              #
#  - Expected output:                                                          #
#    (1)   myFile_aam.smiles                                                   #
#    (2)   myFile_unsuitable.smiles                                            #
#    (3)   myFile_unbalanced.smiles                                            #
#    (4)   myFile_uneven.smiles                                                #
#                                                                              #
#  - Notes:                                                                    #
#                                                                              #
#    * the output file *_aam.similes is the input of EEquAAM.py                #
#                                                                              #
#    * please note that the output file is a plain-text file and so you have   #
#      to zoom it out to PROPERLY see each of the lines in the output          #
#                                                                              #
#    * WARNING! ONLY the balanced reactions that where mapped by ALL the 3     #
#      mappers are considered to be suitable by this script. This is done in   #
#      this way since incomplete maps may depend on the respective mapper      #
#      methodology (RDT, RXN, Graphormer) over such reactions, or more         #
#      IMPORTANTLY there may be a problem with the underlying chemistry        #
#                                                                              #
#    * WARNING! you may use uneven maps under your own risk and consideration  #
#                                                                              #
#    * WARNING! reactions with same identifiers will be silently overwritten   #
#                                                                              #
#    * WARNING! identifiers without reactions wont be considered for analysis  #
#                                                                              #
#    * WARNING! symbol "*" cannot be processed by Graphormer, and symbol "~"   #
#      cannot be processed by RDT Mapper, on the implemented versions          #
#                                                                              #
#    * WARNING! annotations of the type "|1^...|" will be silently removed     #
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


# Requires / Used Versions #####################################################
"""
> Language: python 3.9.13
> Anaconda: conda 22.9.0
> Required packages could be manually installed preferably in the following order
> Packages installed with anaconda's pip:
***** networkx 2.8.4
***** pysmiles 1.0.2
***** rxnmapper 0.2.4
***** rdkit 2022.9.3 (required by rxnmapper)
***** chytorch-rxnmap 1.3 
      (https://github.com/chython/chytorch-rxnmap)
***** chython 1.54
      (https://chython.readthedocs.io/en/latest/)
***** chytorch 1.27 
      (required by chytorch-rxnmap; https://github.com/chython/chytorch)
***** transformers 4.24.0 
      (mind the version by running: pip install transformers==4.24.0)
***** torch 1.11.0
> Java (OpenJDK) and Ubuntu [java --version]
***** openjdk 11.0.17 2022-10-18
***** OpenJDK Runtime Environment (build 11.0.17+8-post-Ubuntu-1ubuntu222.04)
***** OpenJDK 64-Bit Server VM (build 11.0.17+8-post-Ubuntu-1ubuntu222.04, mixed mode)
"""


# Dependencies #################################################################


# installed with conda's pip ---------------------------------------------------
import networkx as nx
import pysmiles as ps
from chython import smiles
from rxnmapper import RXNMapper


# already in python ------------------------------------------------------------
import os
from sys import argv, exit
from copy import deepcopy
from math import factorial, modf


# only report erros and not deprecation warnings -------------------------------
import warnings
import logging
import transformers
warnings.filterwarnings("ignore")
transformers.logging.set_verbosity_error()
logging.getLogger("pysmiles").setLevel(logging.CRITICAL)


# variables ####################################################################


# input ------------------------------------------------------------------------
inputFileName = ""
inputFile = None


# check input ------------------------------------------------------------------
if(".smiles" in argv[1]):
    remainder = (argv[1].split(".smiles"))[-1]
    if(not remainder == ""):
        errorStr = "\n >> MappingTool: Wrong input extension.\n"
        errorStr = errorStr + "- Expected: *.smiles\n"
        errorStr = errorStr + "- Received: *.smiles" + remainder + "\n"
        exit(errorStr)
    else:
        inputFileName = argv[1]
else:
    exit("\n >> MappingTool: missing *.smiles extension.\n")

    
# output -----------------------------------------------------------------------
outputFileName = inputFileName.replace(".smiles", "_aam.smiles")
outputFileNameBad = inputFileName.replace(".smiles", "_unsuitable.smiles")
outputFileNameUnbalanced = inputFileName.replace(".smiles", "_unbalanced.smiles")
outputFileNameUneven = inputFileName.replace(".smiles", "_uneven.smiles")
outputStr = ""
fileContent = []
inputSMILES = dict()
outputFile = None


# data -------------------------------------------------------------------------
RXNaam = ""
RDTaam = ""
CHYaam = ""
mapsRXN = dict()
mapsRDT = dict()
mapsCHY = dict()
eachLineTuple = ()
consistencyRXN = False
consistencyRDT = False
consistencyCHY = False
myrxn = None
myrxnUnmapped = ""
rxn_mapper = None
badInput = ""
badRXN = ""
badRDT = ""
badCHY = ""
currentReaction = ""
originalSMILES = dict()
badSMILES = dict()
unbalancedSMILES = dict()
incompleteMaps = []
incompleteRXN = dict()
incompleteRDT = dict()
incompleteCHY = dict()


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


# function: determin if a given reaction SMILES already has a map --------------
def isAnnotated(someSMILES):
    # local variables
    mol = None
    reactants = []
    products = []
    previousMap = dict()
    annotatedAns = False        
    # check if reactants are annotated
    reactants = list((someSMILES.split(">>")[0]).split("."))
    for eachReactant in reactants:
        # read smiles
        mol = ps.read_smiles(eachReactant)
        # determine if annotated
        previousMap = nx.get_node_attributes(mol, "class")
        if(len(list(previousMap.keys())) > 0):
            annotatedAns = True
            break                        
    # check if products are annotated
    if(not annotatedAns):
        products = list((someSMILES.split(">>")[1]).split("."))
        for eachProduct in products:
            # read smiles
            mol = ps.read_smiles(eachProduct)
            # determine if annotated
            previousMap = nx.get_node_attributes(mol, "class")
            if(len(list(previousMap.keys())) > 0):
                annotatedAns = True
                break                                
    # end of function
    return(annotatedAns)

    
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
        # the only atoms without element are "*" unknowns and thus they are not in the keys
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
        # the only atoms without element are "*" unknowns and thus they are not in the keys        
        for atom in list(chemElement.keys()):        
            if(not(chemElement[atom] in ignoreElements)):
                if(chemElement[atom] in list(typesP.keys())):
                    typesP[chemElement[atom]] = typesP[chemElement[atom]] + 1
                else:
                    typesP[chemElement[atom]] = 1
    # evaluate consistency of nonH atoms
    intersecTypes = list(set(typesR.keys()).intersection(set(typesP.keys())))
    # case 1: element types appear, disappear or change, e.g. CN>>CNO or CNS>>CNP
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


# function: determin if a given atom map completely maps all atoms -------------
def isEven(someMap):
    # local variables
    reactants = []
    products = []
    mapAttribute = dict()
    mappingR = []
    mappingP = []
    atomsR = 0
    atomsP = 0
    completeAns = True
    # get mapping on reactants side
    reactants = list((someMap.split(">>")[0]).split("."))
    for eachReactant in reactants:
        mol = ps.read_smiles(eachReactant)
        mapAttribute = nx.get_node_attributes(mol, "class")
        for atom in list(mol.nodes()):
            atomsR = atomsR + 1
            if(atom in list(mapAttribute.keys())):
                mappingR.append(mapAttribute[atom])
    mappingR.sort()
    # get mapping on products side
    products = list((someMap.split(">>")[1]).split("."))
    for eachProduct in products:
        mol = ps.read_smiles(eachProduct)
        mapAttribute = nx.get_node_attributes(mol, "class")
        for atom in list(mol.nodes()):
            atomsP = atomsP + 1            
            if(atom in list(mapAttribute.keys())):
                mappingP.append(mapAttribute[atom])
    mappingP.sort()
    # determine consistency
    if(not mappingR == mappingP):
        completeAns = False
    else:
        if(not(len(mappingR) == atomsR and len(mappingR) == atomsP)):
            completeAns = False            
    # end of function
    return(completeAns)


# Main #########################################################################


# initial message
print("\n")
print(">>> MappingTool - EEquAAM Github Repository")


# task message
print("\n")
print("* retrieving input file ...")


# get reactions from input
inputFile = open(inputFileName, "r")
inputLines = inputFile.read().splitlines()
inputFile.close()


# task message
print("\n")
print("* collecting the SMILES and their identifiers ...")


# read smiles and their identifiers
currentReaction = ""
for eachLine in inputLines:
    # split line
    eachLineTuple = tuple(eachLine.split(","))
    # get identifier
    if(eachLineTuple[0] == "#"):
        currentReaction = ",".join(eachLineTuple[1:])
        inputSMILES[currentReaction] = ""
        originalSMILES[currentReaction] = ""
    # save reaction
    else:
        # remove annotations of the type "|1^...|"
        if(" " in eachLine):
            inputSMILES[currentReaction] = eachLine.split(" ")[0]
        else:
            inputSMILES[currentReaction] = eachLine        
        originalSMILES[currentReaction] = eachLine        


# remove identifirs without reactions (if any)
theList = list(inputSMILES.keys())
for eachReaction in theList:
    if(inputSMILES[eachReaction] == ""):
        print("--- Warning - found identifier without reaction - removing it:")
        print("--> " + eachReaction)
        inputSMILES.pop(eachReaction)
        
        
# task message
print("\n")
print("* received " + str(len(list(inputSMILES.keys()))) + " reactions ...")


# task message
print("\n")
print("* determining balanced and suitable reactions ...")


# looking for balanced and suitable reactions
i = 0
theList = list(inputSMILES.keys())
for eachReaction in theList:
    if(not ">>" in inputSMILES[eachReaction]):
        badSMILES[eachReaction] = deepcopy(originalSMILES[eachReaction])
        inputSMILES.pop(eachReaction)        
    else:
        reactantSide = inputSMILES[eachReaction].split(">>")[0]
        productSide = inputSMILES[eachReaction].split(">>")[1]
        if((reactantSide == "") or (productSide == "")):        
            badSMILES[eachReaction] = deepcopy(originalSMILES[eachReaction])
            inputSMILES.pop(eachReaction)        
        else:
            annotated = isAnnotated(inputSMILES[eachReaction])
            if(annotated):
                badSMILES[eachReaction] = deepcopy(originalSMILES[eachReaction])
                inputSMILES.pop(eachReaction)        
            else:
                balanced = isBalanced(inputSMILES[eachReaction])        
                if(not balanced):
                    unbalancedSMILES[eachReaction] = deepcopy(originalSMILES[eachReaction])        
                    inputSMILES.pop(eachReaction)
    # print progress
    i = i + 1
    printProgress(round(i*100/len(theList), 2), i, len(theList))


# task message
if(len(list(badSMILES.keys())) > 0):
    print("\n")
    print("* found " + str(len(list(badSMILES.keys()))) + " UNSUITABLE reaction SMILES.")
    print("  - these may be missing the symbol (>>), they may have no reactants or no products,")
    print("  - or they may be already annotated, so they won't be processed")
    print("  - they will be included in the file *_unsuitable.similes")
    print("  - for more information see accompanying README")


# task message
if(len(list(unbalancedSMILES.keys())) > 0):
    print("\n")
    print("* found " + str(len(list(unbalancedSMILES.keys()))) + " UNBALANCED reaction SMILES.")    
    print("  - these cannot provide a (bijectve) atom-to-atom map")
    print("  - they will be included in the file *_unbalanced.similes")
    print("  - for more information see accompanying README")


# task message
if(len(list(inputSMILES.keys())) > 0):
    print("\n")
    print("* analyzing " + str(len(list(inputSMILES.keys()))) + " suitable and balanced reaction SMILES...")


# task message
if(len(list(inputSMILES.keys())) > 0):
    print("\n")
    print("* applying RXN mapper ...")


# call RXNMapper (for each individual reaction so as to reduce memory consumption)
i = 0
theList = list(inputSMILES.keys())
for eachReaction in theList:
    # reinitialize RNXMapper
    rxn_mapper = RXNMapper()
    # map reaction if possible
    try:
        RXNaam = rxn_mapper.get_attention_guided_atom_maps([inputSMILES[eachReaction]], canonicalize_rxns = False)[0]["mapped_rxn"]
        if(" " in RXNaam):
            mapsRXN[eachReaction] = RXNaam.split(" ")[0]
        else:
            mapsRXN[eachReaction] = RXNaam        
    except:
        mapsRXN[eachReaction] = deepcopy(inputSMILES[eachReaction])
    # print progress    
    i = i + 1
    printProgress(round(i*100/len(theList), 2), i, len(theList))    


# task message
if(len(list(inputSMILES.keys())) > 0):
    print("\n")
    print("* applying RDT mapper ...")
    print("* (takes the longest time out of the three mappers)")


# call RDT (in silet mode)
i = 0
theList = list(inputSMILES.keys())
for eachReaction in theList:
    # map reaction if possible
    try:
        os.system("nohup java -jar RDT_2.4.1.jar -Q SMI -q " + "\"" + inputSMILES[eachReaction] + "\"" + " -c -j AAM -f TEXT >/dev/null 2>&1")
        inputFile = open("ECBLAST_smiles_AAM.txt", "r")
        RDTaam = (inputFile.read().splitlines())[3]
        inputFile.close()
        os.system("rm ECBLAST_smiles_AAM.txt")
        os.system("rm ECBLAST_smiles_AAM.rxn")
        if(" " in RDTaam):
            mapsRDT[eachReaction] = RDTaam.split(" ")[0]
        else:
            mapsRDT[eachReaction] = RDTaam
    except:
        mapsRDT[eachReaction] = deepcopy(inputSMILES[eachReaction])
    # print progress
    i = i +1
    printProgress(round(i*100/len(theList), 2), i, len(theList))    


# task message
if(len(list(inputSMILES.keys())) > 0):
    print("\n")
    print("* applying Graphormer mapper ...")


# call Graphormer
i = 0
theList = list(inputSMILES.keys())
for eachReaction in theList:
    # map reaction if possible
    try:
        myrxn = smiles(inputSMILES[eachReaction])
        myrxn.reset_mapping()
        CHYaam = format(myrxn, "m")
        if(" " in CHYaam):
            mapsCHY[eachReaction] = CHYaam.split(" ")[0]
        else:       
            mapsCHY[eachReaction] = CHYaam
    except:
        mapsCHY[eachReaction] = deepcopy(inputSMILES[eachReaction])
    # print progress
    i = i + 1
    printProgress(round(i*100/len(theList), 2), i, len(theList))    
   

# task message
if(len(list(inputSMILES.keys())) > 0):
    print("\n")
    print("* determining reactions with the 3 complete maps (RXN, RDT and CHY) ...")


# checking all atoms were mapped by all mappers
theList = list(inputSMILES.keys())
for eachReaction in theList:
    # get consistency of maps
    consistencyRXN = isEven(mapsRXN[eachReaction])
    consistencyRDT = isEven(mapsRDT[eachReaction])
    consistencyCHY = isEven(mapsCHY[eachReaction])
    # evaluate consistency of all mappers
    if(not(consistencyRXN and consistencyRDT and consistencyCHY)):
        incompleteMaps.append(eachReaction)
        incompleteRXN[eachReaction] = deepcopy(mapsRXN[eachReaction])
        incompleteRDT[eachReaction] = deepcopy(mapsRDT[eachReaction])
        incompleteCHY[eachReaction] = deepcopy(mapsCHY[eachReaction])


# remove unsuitable reactions
for eachReaction in incompleteMaps:
    inputSMILES.pop(eachReaction)
    mapsRXN.pop(eachReaction)
    mapsRDT.pop(eachReaction)
    mapsCHY.pop(eachReaction)
    

# task message
print("\n")
if(len(list(inputSMILES.keys())) > 0):
    print("* saving results obtained with RXN, RDT and Graphormer ...")    
    if(len(incompleteMaps) > 0):    
        print("- reactions with the 3 complete maps (*_aam.smiles): " + str(len(list(inputSMILES.keys()))))
        print("- reactions with some incomplete maps (*_uneven.simles): " + str(len(incompleteMaps)))
    else:
        print("- all suitable and balanced reactions were completely mapped (*_aam.smiles): " + str(len(list(inputSMILES.keys()))))
else:
    print("- no suitable reactions with 3 complete maps were found; *_aam.similes file won't be generated")
    

# print results suitable
if(len(list(inputSMILES.keys())) > 0):
    fileContent = []
    outputStr = ""    
    for eachReaction in list(inputSMILES.keys()):        
        fileContent.append("#," + eachReaction + "\n")
        fileContent.append("RXNmap," + mapsRXN[eachReaction] + "\n")
        fileContent.append("RDTmap," + mapsRDT[eachReaction] + "\n")
        fileContent.append("CHYmap," + mapsCHY[eachReaction] + "\n")
    outputStr = "".join(fileContent)
    outputFile = open(outputFileName, "w")
    outputFile.write(outputStr)
    outputFile.close()
        

# print results unsuitable
if(len(list(badSMILES.keys())) > 0):
    fileContent = []
    outputStr = ""    
    for eachReaction in list(badSMILES.keys()):        
        fileContent.append("#," + eachReaction + "\n")
        fileContent.append(originalSMILES[eachReaction] + "\n")
    outputStr = "".join(fileContent)
    outputFile = open(outputFileNameBad, "w")
    outputFile.write(outputStr)
    outputFile.close()


# print results unbalanced
if(len(list(unbalancedSMILES.keys())) > 0):
    fileContent = []
    outputStr = ""    
    for eachReaction in list(unbalancedSMILES.keys()):        
        fileContent.append("#," + eachReaction + "\n")
        fileContent.append(originalSMILES[eachReaction] + "\n")
    outputStr = "".join(fileContent)
    outputFile = open(outputFileNameUnbalanced, "w")
    outputFile.write(outputStr)
    outputFile.close()

    
# print results uneven
if(len(incompleteMaps) > 0):
    fileContent = []
    outputStr = ""    
    for eachReaction in incompleteMaps:
        fileContent.append("#," + eachReaction + "\n")
        fileContent.append("RXNmap," + incompleteRXN[eachReaction] + "\n")
        fileContent.append("RDTmap," + incompleteRDT[eachReaction] + "\n")
        fileContent.append("CHYmap," + incompleteCHY[eachReaction] + "\n")
    outputStr = "".join(fileContent)
    outputFile = open(outputFileNameUneven, "w")
    outputFile.write(outputStr)
    outputFile.close()


# final message
print("\n")
print(">>> Finished")
print("\n")


# End ##########################################################################
################################################################################
