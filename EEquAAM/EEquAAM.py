################################################################################
#                                                                              #
#  README - Program: EEquAAM.py                                                #
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
#  - Description: this script receives a plain-text file with *_aam.smiles     #
#    format containing COMPLETELY mapped reaction SMILES and identifiers for   #
#    them, as shown in the following toy-example:                              #
#                                                                              #
#    #,reaction_1                                                              #
#    RXNmap,[CH3:1][...][CH:8]>>[CH3:1][...][CH3:9].[CH4:8]                    #
#    RDTmap,[CH3:2][...][CH:4]>>[CH3:1][...][CH3:9].[CH4:2]                    #
#    CHYmap,[CH3:6][...][CH:7]>>[CH3:7][...][CH3:1].[CH4:6]                    #
#    #,reaction_2                                                              #
#    RXNmap,[CH3:1][...][CH:8]>>[CH3:1][...][CH3:4].[CH3:5][...][CH4:8]        #
#    RDTmap,[CH3:2][...][CH:4]>>[CH3:1][...][CH3:9].[CH4:2][...][CH4:8]        #
#    CHYmap,[CH3:6][...][CH:7]>>[CH3:7][...][CH3:1].[CH3:8][...][CH3:4]        #
#    ...                                                                       #
#                                                                              #
#    * WARNING! if the program founds an INCOMPLETE map it will terminate      #
#      before analyzing any reaction. An incomplete atom-map is that in which  #
#      there are atoms on one side of the reaction not being mapped to atoms   #
#      on the  other side, in other words, it is not a bijective atom-map.     #
#    * the line with the identifier should always start with "#,"              #
#    * there is no limit to the number of reactions that can be analized       #
#    * the number of maps per rxn doesn't need to be the same for all rxns     #
#    * identifiers of the mappers [RXNmap...] can be changend if needed        #
#    * the id "reaction_i" is an arbitrary identifier and can be changed to be #
#      e.g. the common name of a reaction, including more commas if needed     #
#      but not newline characters "\n", i.e, this is a single line             #
#    * this format is the output of the MappingTool found in the EEquAAM repo  #
#      already implementing RXNmapper, RDT and GraphormerMapper                #
#    * general specifications on both anotated and unmapped SMILES strings     #
#      can be found at: http://opensmiles.org/opensmiles.html                  #
#                                                                              #
#    The program can run a "--sanity-check" involving the three mathematically #
#    equivalent methods for comparing atom maps: (AUX) comparing the           #
#    auxiliary graphs, (CGR or ITS) comparing Fujita's (1980) imaginary        #
#    transition state and (ISO) enumerating the isomorphisms associated to     #
#    each reaction. The default mode only runs the ITS (the fastest method).   #
#    The output includes a summary, three pkls, and a pdf with the boxplots    #
#    of the time taken by each method to analyze each reaction once.           #
#                                                                              #
#  - Input: plain-text csv *_aam.smiles file as specified above.               #
#                                                                              #
#  - Output: (1) pdf with the mentioned time (seconds) boxplots. (2) pkl with  #
#    the auxiliary graphs and equivalence classes found by AUX. (3) pkl with   #
#    the condensed graphs of the reaction and the equivalence classes found by #
#    the CGR method. (4) pkl with the reactants graph G, the products graph H  #
#    and the equivalence classes found by ISO. (5) a txt file with a summary   #
#    indicating if the atom maps where equivalent or not under each method, as #
#    well as the number of equivalence classes found and the time (seconds)    #
#    taken by each method for each reaction. Note that the equivalence classes #
#    of the atom maps can be retrieved from the PKL's. The time boxplots show  #
#    Log_10 of time; linear time is always included in the summary.            #
#                                                                              #
#  - Run with (after activating eequaam conda environment):                    #
#      * default ITS  python EEquAAM.py [myFile_aam.smiles]                    #
#      * or           python EEquAAM.py --sanity-check [myFile_aam.smiles]     #
#                                                                              #
#  - Expected output:                                                          #
#    (1) myFile_aam_times.pdf                                                  #
#    (2) myFile_aam_aux.pkl      (only with --sanity-check option)             #
#    (3) myFile_aam_its.pkl                                                    #
#    (4) myFile_aam_iso.pkl      (only with --sanity-check option)             #
#    (5) myFile_aam_summary.txt  (tab separated file, may be changed in-code)  #
#                                                                              #
#  - Notes:                                                                    #
#                                                                              #
#    * WARNING! if running --sanity-check take into account that ISO's running #
#      time grows exponentially w.r.t the number of symmetries in the          #
#      molecules, leading to minutes or hours (or more time) required for      #
#      completion. On the other hand AUX and ITS (or CGR) have more tolerance  #
#      for these cases. For a full discussion on this see the cited paper.     #
#                                                                              #
#    * WARNING! if the program founds an INCOMPLETE map it will terminate      #
#      before analyzing any reaction                                           #
#                                                                              #
#    * WARNING! identifiers without maps wont be considered for analysis       #
#                                                                              #
#                                                                              #
#  - DISCLAIMER: This code is provided "AS IS". You may use it uder your own   #
#    risk and consideration. This script was built for research purposes and   #
#    its developers claim no responsabilty on further third-party use.         #
#                                                                              #
################################################################################


# Requires / Used Versions #####################################################
"""
> Lenguaje: python 3.9.13
> Anaconda: conda 22.9.0
> Packages installed with anaconda:
***** pysmiles 1.0.2
***** networkx 2.8.4
***** matplotlib 3.5.2
"""


# Dependencies #################################################################


# installed with conda ---------------------------------------------------------
import pysmiles as ps
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from networkx.algorithms import isomorphism
from networkx.algorithms.isomorphism import generic_node_match, generic_edge_match


# already in python ------------------------------------------------------------
import time
import pickle
from sys import argv, exit
from operator import eq
from copy import deepcopy
from math import factorial, modf, log10


# turn off warnings and matplotlib to run in the server ------------------------
import warnings
import logging
matplotlib.use("Agg")
warnings.filterwarnings("ignore")
logging.getLogger("pysmiles").setLevel(logging.CRITICAL)


# variables ####################################################################


# input ------------------------------------------------------------------------
inputFileName = ""
inputSMILES = []
allMaps = []
inputFile = None
sanityCheck = False


# check user input options -----------------------------------------------------
if(len(argv) in [2, 3]):
    if(len(argv) == 2):
        if(".smiles" in argv[1]):
            inputFileName = argv[1]
            sanityCheck = False
        else:
            exit("\n >> EEquAAM: Wrong input format.\n")
    if(len(argv) == 3):
        if((argv[1] == "--sanity-check") and (".smiles" in argv[2])):
            inputFileName = argv[2]
            sanityCheck = True
        else:
            exit("\n >> EEquAAM: Wrong input format.\n")                    
else:
    exit("\n >> EEquAAM: Wrong input format.\n")


# output -----------------------------------------------------------------------
outputFileNamePklAUX = inputFileName.replace(".smiles", "_aux.pkl")
outputFileNamePklCGR = inputFileName.replace(".smiles", "_its.pkl")
outputFileNamePklISO = inputFileName.replace(".smiles", "_iso.pkl")
outputFileNameSummary = inputFileName.replace(".smiles", "_summary.txt")
outputFileNameBoxPlot = inputFileName.replace(".smiles", "_times.pdf")
outputFile = None
summaryOutput = ""
subsetsOutput = ""


# data -------------------------------------------------------------------------
productsSide = ""
reactantsSide = ""
currentReaction = ""
tempResults = []
productsList = []
reactantsList = []
eachLineTuple = ()
resultsAUX = dict()
resultsCGR = dict()
resultsISO = dict()
graphsByMap = dict()
mappedSMILES = dict()
originalSMILES = dict()
# graphs
G = None
H = None
molecule = None
# times
finalTime = 0
initialTime = 0
timeAUX = dict()
timeCGR = dict()
timeISO = dict()


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
    

# function: obtain and compare auxiliary graphs --------------------------------
def analyzeAuxiliaryGraphs(someResults):
    # local variables
    theClass = 0
    lastClass = 0
    someClass = 0    
    someID = ""
    someMap = ""
    classified = []        
    someTempAUX = []
    nodeLabelNames = []
    nodeLabelDefault = []
    nodeLabelOperator = []
    foundEquivalent = False
    nodeMatch = None
    edgeMatch = None        
    auxG = None
    auxH = None    
    someG = None
    someH = None
    someAUX = None    
    namesG = dict()
    namesH = dict()
    coordG = dict()
    coordH = dict()
    # obtain auxiliary graphs
    for (someID, someMap, someG, someH) in someResults:
        # copy original graphs
        auxG = deepcopy(someG)
        auxH = deepcopy(someH)
        # add property to preserve during iso check
        coordG = {v:1 for v in list(auxG.nodes())}
        coordH = {v:2 for v in list(auxH.nodes())}
        nx.set_node_attributes(auxG, coordG, "GorH")
        nx.set_node_attributes(auxH, coordH, "GorH")
        # rename nodes to make (disjoint) union
        namesG = {v:(v, 1) for v in list(auxG.nodes())}
        namesH = {v:(v, 2) for v in list(auxH.nodes())}
        auxG = deepcopy(nx.relabel_nodes(auxG, namesG))
        auxH = deepcopy(nx.relabel_nodes(auxH, namesH))
        # make (disjoint) union
        someAUX = deepcopy(nx.union(auxG, auxH))
        # add the edges representing the aam
        for v in list(someG.nodes()):
            someAUX.add_edge((v, 1), (v, 2), order = "*")
        # save auxiliary graph
        someTempAUX.append((someID, someMap, someG, someH, someAUX))
    # compare auxiliary graphs
    nodeLabelNames = ["GorH", "element", "aromatic", "hcount", "charge"]
    nodeLabelDefault = [0, "*", False, 0, 0]
    nodeLabelOperator = [eq, eq, eq, eq, eq]
    nodeMatch = generic_node_match(nodeLabelNames, nodeLabelDefault, nodeLabelOperator)
    edgeMatch = generic_edge_match("order", 1, eq)
    for (someID, someMap, someG, someH, someAUX) in someTempAUX:
        # restart control variables
        foundEquivalent = False
        # compare with classified auxiliary graphs (if any yet)
        for (someClass, someIDc, someMapc, someGc, someHc, someAUXc) in classified:
            if(nx.is_isomorphic(someAUX, someAUXc, node_match = nodeMatch, edge_match = edgeMatch)):
                foundEquivalent = True
                theClass = someClass
                break
        # if not equivalent make new mapping class
        if(not foundEquivalent):
            theClass = lastClass + 1
            lastClass = theClass
        # save results
        classified.append((theClass, someID, someMap, someG, someH, someAUX))
    # end of function    
    return(classified)


# function: obtain and compare condensed graph of the reaction -----------------
def analyzeCGRs(someResults):
    # local variables
    theClass = 0
    lastClass = 0
    someClass = 0    
    someID = ""
    someMap = ""
    edgesG = []
    edgesH = []
    edgesCGR = []    
    classified = []        
    someTempCGR = []
    nodeLabelNames = []
    nodeLabelDefault = []
    nodeLabelOperator = []
    typesG = ()
    typesH = ()
    typesDict = dict()
    foundEquivalent = False
    nodeMatch = None
    edgeMatch = None        
    someG = None
    someH = None
    someCGR = None    
    # obtain CGRs
    for (someID, someMap, someG, someH) in someResults:
        # creat null from copy of G in order to preserve attributes
        someCGR = deepcopy(someG)
        someCGR.remove_edges_from(list(someCGR.edges()))
        # add typeG and typeH attributes, or default attributes for "*" unknown elements
        for v in list(someCGR.nodes()):
            # get attributes in G or default if missing
            try:
                elementAttr = someG.nodes[v]["element"]
            except:
                elementAttr = "*"
            try:
                aromaticAttr = someG.nodes[v]["aromatic"]
            except:
                aromaticAttr = False
            try:
                hcountAttr = someG.nodes[v]["hcount"]
            except:
                hcountAttr = 0
            try:
                chargeAttr = someG.nodes[v]["charge"]
            except:
                chargeAttr = 0                
            typesG = (elementAttr, aromaticAttr, hcountAttr, chargeAttr)
            # get attributes in H or default if missing
            try:
                elementAttr = someH.nodes[v]["element"]
            except:
                elementAttr = "*"
            try:
                aromaticAttr = someH.nodes[v]["aromatic"]
            except:
                aromaticAttr = False
            try:
                hcountAttr = someH.nodes[v]["hcount"]
            except:
                hcountAttr = 0
            try:
                chargeAttr = someH.nodes[v]["charge"]
            except:
                chargeAttr = 0                
            typesH = (elementAttr, aromaticAttr, hcountAttr, chargeAttr)
            # make new label to be preserved
            typesDict[v] = (typesG, typesH)
        nx.set_node_attributes(someCGR, typesDict, "typesGH")        
        # copy original graphs
        edgesG = list(someG.edges())
        edgesH = list(someH.edges())
        # add edges from G
        for (u, v) in edgesG:
            # check if not already an edge of CGR
            edgesCGR = list(someCGR.edges())
            if((not (u, v) in edgesCGR) and (not (v, u) in edgesCGR)):
                # check if also an edge in H
                if(((u, v) in edgesH) or ((v, u) in edgesH)):
                    if((u, v) in edgesH):
                        someLabel = (someG[u][v]["order"], someH[u][v]["order"])
                        someCGR.add_edge(u, v, order = someLabel)                        
                    if((v, u) in edgesH):
                        someLabel = (someG[u][v]["order"], someH[v][u]["order"])
                        someCGR.add_edge(u, v, order = someLabel)
                # if not then add "*" to the H-part of the label 
                else:
                    someLabel = (someG[u][v]["order"], "*")
                    someCGR.add_edge(u, v, order = someLabel)                        
        # add edges from H
        for (u, v) in edgesH:
            # check if not already an edge of CGR
            edgesCGR = list(someCGR.edges())
            if((not (u, v) in edgesCGR) and (not (v, u) in edgesCGR)):
                # check if also an edge in G
                if(((u, v) in edgesG) or ((v, u) in edgesG)):
                    if((u, v) in edgesG):
                        someLabel = (someG[u][v]["order"], someH[u][v]["order"])
                        someCGR.add_edge(u, v, order = someLabel)                        
                    if((v, u) in edgesG):
                        someLabel = (someG[v][u]["order"], someH[u][v]["order"])
                        someCGR.add_edge(u, v, order = someLabel)
                # if not then add "*" to the G-part of the label 
                else:
                    someLabel = ("*", someH[u][v]["order"])
                    someCGR.add_edge(u, v, order = someLabel)                        
        # save auxiliary graph
        someTempCGR.append((someID, someMap, someG, someH, someCGR))        
    # compare CGRs
    nodeLabelNames = ["element", "aromatic", "hcount", "charge", "typesGH"]
    nodeLabelDefault = ["*", False, 0, 0, ()]
    nodeLabelOperator = [eq, eq, eq, eq, eq]
    nodeMatch = generic_node_match(nodeLabelNames, nodeLabelDefault, nodeLabelOperator)
    edgeMatch = generic_edge_match("order", 1, eq)   
    for (someID, someMap, someG, someH, someCGR) in someTempCGR:
        # restart control variables
        foundEquivalent = False
        # compare with classified CGRs (if any yet)
        for (someClass, someIDc, someMapc, someGc, someHc, someCGRc) in classified:
            if(nx.is_isomorphic(someCGR, someCGRc, node_match = nodeMatch, edge_match = edgeMatch)):
                foundEquivalent = True
                theClass = someClass
                break
        # if not equivalent make new mapping class
        if(not foundEquivalent):
            theClass = lastClass + 1
            lastClass = theClass
        # save results
        classified.append((theClass, someID, someMap, someG, someH, someCGR))
    # end of function    
    return(classified)


# function: obtain and compare isomorphisms ------------------------------------
def analyzeISOs(someResults):
    # local variables
    theClass = 0
    lastClass = 0
    someClass = 0    
    someID = ""
    someMap = ""
    classified = []
    nodeLabelNames = []
    nodeLabelDefault = []
    nodeLabelOperator = []
    foundEquivalent = False
    nodeMatch = None
    edgeMatch = None        
    someG = None
    someH = None
    isoMatcherG = None
    isoMatcherH = None    
    listIsosG = []
    listIsosH = []    
    # obtain and compare ISOs
    nodeLabelNames = ["element", "aromatic", "hcount", "charge"]
    nodeLabelDefault = ["*", False, 0, 0]
    nodeLabelOperator = [eq, eq, eq, eq]
    nodeMatch = generic_node_match(nodeLabelNames, nodeLabelDefault, nodeLabelOperator)
    edgeMatch = generic_edge_match("order", 1, eq)   
    for (someID, someMap, someG, someH) in someResults:
        # restart control variables
        foundEquivalent = False
        # compare with classified graph's isos (if any yet)
        for (someClass, someIDc, someMapc, someGc, someHc) in classified:
            isoMatcherG = isomorphism.GraphMatcher(someG, someGc, node_match = nodeMatch, edge_match = edgeMatch)
            isoMatcherH = isomorphism.GraphMatcher(someH, someHc, node_match = nodeMatch, edge_match = edgeMatch)
            listIsosG = list(isoMatcherG.isomorphisms_iter())
            listIsosH = list(isoMatcherH.isomorphisms_iter())                        
            for isoG in listIsosG:
                if(isoG in listIsosH):
                    foundEquivalent = True
                    theClass = someClass
                    break
            if(foundEquivalent):
                break
        # if not equivalent make new mapping class
        if(not foundEquivalent):
            theClass = lastClass + 1
            lastClass = theClass
        # save results
        classified.append((theClass, someID, someMap, someG, someH))
    # end of function    
    return(classified)


# Main #########################################################################


# initial message
print("\n")
print(">>> EEquAAM - EEquAAM Github Repository")


# sanity check message
if(sanityCheck):
    print("\n+++++ received --sanity-check argument")
    print("+++++ running the three methods: AUX, ITS and ISO")
    print("+++++ take into account that the ISO method may be more time consuming")
else:
    print("\n+++++ running only ITS method (default)")
    print("+++++ for other options see accompanying README")
    

# task message
print("\n")
print("* retreiving input file ...")


# load mapped smiles file (*_aam.smiles)
inputFile = open(inputFileName, "r")
inputSMILES = inputFile.read().splitlines()
inputFile.close()


# task message
print("\n")
print("* collecting SMILES and their identifiers ...")


# read smiles and their identifiers
currentReaction = ""
for eachLine in inputSMILES:
    # split line
    eachLineTuple = tuple(eachLine.split(","))
    # get identifier
    if(eachLineTuple[0] == "#"):
        currentReaction = ",".join(eachLineTuple[1:])
        mappedSMILES[currentReaction] = []
        originalSMILES[currentReaction] = []
    # save reaction
    else:
        allMaps.append([currentReaction, eachLineTuple])
        mappedSMILES[currentReaction].append(eachLineTuple)
        originalSMILES[currentReaction].append(eachLineTuple)


# remove identifirs without maps (if any)
for eachReaction in list(mappedSMILES.keys()):
    if(len(mappedSMILES[eachReaction]) == 0):
        print("--- Warning - found identifier without maps - removing it:")
        print("--> " + eachReaction)
        mappedSMILES.pop(eachReaction)

        
# task message
print("\n")
print("* received " + str(len(list(mappedSMILES.keys()))) + " reaction SMILES ...")
print("* and a total of " + str(len(allMaps)) + " maps over these reactions ...")


# task message
print("\n")
print("* evaluating that all the given atom maps are complete (bijective) ...")


# looking for incomplete maps
i = 0
for [eachReaction, eachMap] in allMaps:
    even = isEven(eachMap[1])
    if(not even):        
        answer = "***** Found an incomplete map, i.e., not bijective. Terminating before completing the analysis.\n"
        answer = answer + "***** This cannot be properly compared to other maps, please see accompanying literature.\n"
        answer = answer + "***** The incomplete map is the following...\n"
        answer = answer + "***** ReactionID: " + eachReaction + "\n"
        answer = answer + "***** MapID, Map: " + eachMap[0] + ", " + eachMap[1] + "\n"
        print("\n")
        exit(answer)
    # print progress
    i = i + 1
    printProgress(round(i*100/len(allMaps), 2), i, len(allMaps))

    
# task message
print("\n")
print("* all given maps are complete; building reactants-graph and products-graph ...")

        
# turn mapped smiles into networkx graphs G of reactants and H of products
i = 0
for eachReaction in list(mappedSMILES.keys()):
    tempResults = []
    for eachMap in mappedSMILES[eachReaction]:
        # split reaction
        reactantsSide = (eachMap[1].split(">>"))[0]
        productsSide = (eachMap[1].split(">>"))[1]
        # get reactants graph G
        reactantsList = reactantsSide.split(".")        
        G = nx.Graph()
        for eachReactant in reactantsList:
            # read smiles into networkx graph
            molecule = ps.read_smiles(eachReactant)
            # rename vertices using atom map            
            atomMapping = nx.get_node_attributes(molecule, "class")
            molecule = nx.relabel_nodes(molecule, atomMapping)
            # make (disjoint) union over G            
            G = deepcopy(nx.union(G, molecule))
        # get products graph H
        productsList = productsSide.split(".")        
        H = nx.Graph()
        for eachProduct in productsList:
            # read smiles into networkx graph            
            molecule = ps.read_smiles(eachProduct)
            # rename vertices using atom map            
            atomMapping = nx.get_node_attributes(molecule, "class")
            molecule = nx.relabel_nodes(molecule, atomMapping)
            # make (disjoint) union over H            
            H = deepcopy(nx.union(H, molecule))
        # save map data
        tempResults.append((eachMap[0], eachMap[1], deepcopy(G), deepcopy(H)))
    # save reaction data
    graphsByMap[eachReaction] = deepcopy(tempResults)
    # print progress
    i = i + 1
    printProgress(round(i*100/len(list(mappedSMILES.keys())), 2), i, len(list(mappedSMILES.keys())))


# AUX analysis
if(sanityCheck):
    # task message
    print("\n")
    print("* running AUX: comparing auxiliary graphs of each reaction ...")
    # analyze auxiliary graphs for the maps of each reaction
    count = 0
    for eachReaction in list(graphsByMap.keys()):
        initialTime = time.time()
        resultsAUX[eachReaction] = analyzeAuxiliaryGraphs(graphsByMap[eachReaction])
        finalTime = time.time()    
        timeAUX[eachReaction] = finalTime-initialTime
        count = count + 1
        printProgress(round(count*100/len(list(graphsByMap.keys())), 2), count, len(list(graphsByMap.keys())))


# task message
print("\n")
print("* running ITS: comparing imaginary transition state graphs of each reaction ...")    
# analyze CGRs for the maps of each reaction
count = 0
for eachReaction in list(graphsByMap.keys()):
    initialTime = time.time()
    resultsCGR[eachReaction] = analyzeCGRs(graphsByMap[eachReaction])
    finalTime = time.time()    
    timeCGR[eachReaction] = finalTime-initialTime
    count = count + 1
    printProgress(round(count*100/len(list(graphsByMap.keys())), 2), count, len(list(graphsByMap.keys())))


# ISO analysis
if(sanityCheck):    
    # task message
    print("\n")
    print("* running ISO: comparing isomorphisms associated to each reaction ...")
    print("- WARNING! This method can take more time due to its mathematical properties.")
    # analyze ISOs for the maps of each reaction
    count = 0
    for eachReaction in list(graphsByMap.keys()):
        initialTime = time.time()
        resultsISO[eachReaction] = analyzeISOs(graphsByMap[eachReaction])
        finalTime = time.time()    
        timeISO[eachReaction] = finalTime-initialTime
        count = count + 1
        printProgress(round(count*100/len(list(graphsByMap.keys())), 2), count, len(list(graphsByMap.keys())))
    

# task message
print("\n")
print("* making boxplots ...")

    
# define data to plot
if(sanityCheck):
    logTimeCGR = [log10(t) for t in list(timeCGR.values())]    
    logTimeISO = [log10(t) for t in list(timeISO.values())]
    logTimeAUX = [log10(t) for t in list(timeAUX.values())]
    timeData = [logTimeISO, logTimeAUX, logTimeCGR]
    timeLabels = [r"ISO-$\equiv$", r"AUX-$\Gamma$", r"ITS-$\Upsilon$"]
    whiskers = 1.25
    widthsBoxes = 0.4
    hatchB = ["....", "////", "oo"]    
else:
    logTimeCGR = [log10(t) for t in list(timeCGR.values())]    
    timeData = [logTimeCGR]
    timeLabels = [r"ITS-$\Upsilon$"]
    whiskers = 1.25
    widthsBoxes = 0.4
    hatchB = ["////"]
# make plot
fig, ax = plt.subplots()
bps = ax.boxplot(timeData, whis = whiskers, widths = widthsBoxes, labels = timeLabels, sym = ".", patch_artist = True)
# set colors
for i in range(len(bps["fliers"])):
    bps["fliers"][i].set(markeredgecolor = "grey", markerfacecolor = "grey")
for i in range(len(bps["caps"])):
    bps["caps"][i].set(color = "dimgrey")
for i in range(len(bps["boxes"])):
    bps["boxes"][i].set(facecolor = "w", edgecolor = "grey", hatch = hatchB[i])
for i in range(len(bps["medians"])):    
    bps["medians"][i].set(color = "dimgrey")
for i in range(len(bps["whiskers"])):
    bps["whiskers"][i].set(color = "dimgrey")
# set some attributes and save figures
plt.ylabel(r"$Log_{10}$" + " of running time [s]\n", fontsize = 15)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 12)
plt.title("Analyzed " + str(len(list(graphsByMap.keys()))) + " reactions", fontsize = 12)
plt.grid(color = "lightgray", axis = "y", linestyle = "--", linewidth = 0.7)
plt.tight_layout()
plt.savefig(outputFileNameBoxPlot)
plt.close()


# task message
print("* making summary ...")


# save summary in txt file
dataStr = ""
summaryOutput = ""
sep = "\t"
totEquivalent = 0
totNotEquivalent = 0
# add lines per reaction
for eachReaction in list(graphsByMap.keys()):
    # add reaction smiles
    summaryOutput = summaryOutput + "#," + eachReaction + "\n"
    # add number of input maps for reaction
    dataStr = str(len(graphsByMap[eachReaction]))
    summaryOutput = summaryOutput + "maps given" + sep + dataStr + "\n"
    # add number of classes obtained when comparing auxiliary graphs
    if(sanityCheck):
        dataStr = str(max([classR for (classR, idR, mapR, GR, HR, auxR) in resultsAUX[eachReaction]]))
        cA = int(dataStr)
        summaryOutput = summaryOutput + "classes AUX" + sep + dataStr + "\n"    
    # add number of classes obtained when comparing CGRs
    dataStr = str(max([classR for (classR, idR, mapR, GR, HR, cgrR) in resultsCGR[eachReaction]]))
    cB = int(dataStr)
    summaryOutput = summaryOutput + "classes ITS" + sep + dataStr + "\n"
    # add number of classes obtained when comparing ISOs
    if(sanityCheck):    
        dataStr = str(max([classR for (classR, idR, mapR, GR, HR) in resultsISO[eachReaction]]))
        cC = int(dataStr)
        summaryOutput = summaryOutput + "classes ISO" + sep + dataStr + "\n"    
    # add result AUX
    if(sanityCheck):    
        if(max([classR for (classR, idR, mapR, GR, HR, auxR) in resultsAUX[eachReaction]]) == 1):
            dataStr = "equivalent maps"
            rA = "equivalent maps"
        else:
            dataStr = "not equivalent maps"
            rA = "not equivalent maps"
        summaryOutput = summaryOutput + "result AUX" + sep + dataStr + "\n"
    # add result CGR
    if(max([classR for (classR, idR, mapR, GR, HR, cgrR) in resultsCGR[eachReaction]]) == 1):
        dataStr = "equivalent maps"
        rB = "equivalent maps"
        totEquivalent = totEquivalent + 1        
    else:
        dataStr = "not equivalent maps"
        rB = "not equivalent maps"
        totNotEquivalent = totNotEquivalent + 1
    summaryOutput = summaryOutput + "result ITS" + sep + dataStr + "\n"
    # add result ISO
    if(sanityCheck):    
        if(max([classR for (classR, idR, mapR, GR, HR) in resultsISO[eachReaction]]) == 1):
            dataStr = "equivalent maps"
            rC = "equivalent maps"
        else:
            dataStr = "not equivalent maps"
            rC = "not equivalent maps"
        summaryOutput = summaryOutput + "result ISO" + sep + dataStr + "\n"    
    # add time AUX
    if(sanityCheck):    
        dataStr = str(timeAUX[eachReaction])
        summaryOutput = summaryOutput + "time AUX" + sep + dataStr + "\n"
    # add time CGR        
    dataStr = str(timeCGR[eachReaction])
    summaryOutput = summaryOutput + "time ITS" + sep + dataStr + "\n"
    # add time ISO
    if(sanityCheck):    
        dataStr = str(timeISO[eachReaction])
        summaryOutput = summaryOutput + "time ISO" + sep + dataStr + "\n"
    # consistency info
    if(sanityCheck):
        if(len(list(set([cA, cB, cC]))) > 1):
            print("*** Wrong, please check:", eachReaction, list(set([rA, rB, rC])), list(set([cA, cB, cC])))
# add super-summary at the beginning of summary
superSummary = "+++ Total of reactions with equivalent maps: " + str(totEquivalent) + "\n"
superSummary = superSummary + "+++ Total of reactions with non-equivalent maps: " + str(totNotEquivalent) + "\n"
summaryOutput = superSummary + summaryOutput
# save summary
outputFile = open(outputFileNameSummary, "w")
outputFile.writelines(summaryOutput)
outputFile.close()


# task message
print("* saving data into pkl files ...")


# save pkl data of AUX
if(sanityCheck):
    outputFile = open(outputFileNamePklAUX, "wb")
    pickle.dump(resultsAUX, outputFile)
    outputFile.close()


# save pkl data of CGR
outputFile = open(outputFileNamePklCGR, "wb")
pickle.dump(resultsCGR, outputFile)
outputFile.close()


# save pkl data of ISO
if(sanityCheck):
    outputFile = open(outputFileNamePklISO, "wb")
    pickle.dump(resultsISO, outputFile)
    outputFile.close()


# final message
print("\n")
print(">>> Finished")
print("\n")


# End ##########################################################################
################################################################################
