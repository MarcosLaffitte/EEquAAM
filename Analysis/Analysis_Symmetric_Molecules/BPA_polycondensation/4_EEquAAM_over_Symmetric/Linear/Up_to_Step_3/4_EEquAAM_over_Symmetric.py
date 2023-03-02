################################################################################
#                                                                              #
#  README - Program: 4_EEquAAM_over_Symmetric.py                               #
#                                                                              #
#  - Paper: https://match.pmf.kg.ac.rs/issues/m90n1/m90n1_75-102.html          #
#                                                                              #
#  - Github repository: https://github.com/MarcosLaffitte/EEquAAM              #
#                                                                              #
#  - Contributor(s):                                                           #
#  * @MarcosLaffitte - Marcos E. Gonzalez Laffitte                             #
#                                                                              #
#  - Description: this script receives a tab-separated txt (\t) containing     #
#    both the unmapped and mapped SMILES of a set of reactions, and ids for    #
#    them as in the following toy example:                                     #
#                                                                              #
#    #         CCCCCCCCCC>>CCCCCCCCC.C                                         #
#    RXNmap    [CH3:1][...][CH:10]>>[CH3:1][...][CH3:9].[CH4:10]               #
#    RDTmap    [CH3:2][...][CH:4]>>[CH3:1][...][CH3:9].[CH4:2]                 #
#    CHYmap    [CH3:6][...][CH:7]>>[CH3:7][...][CH3:1].[CH4:6]                 #
#    #         CCCCCCCCCC>>CCCC.CCCCCC                                         #
#    RXNmap    [CH3:1][...][CH:10]>>[CH3:1][...][CH3:4].[CH3:5][...][CH4:10]   #
#    RDTmap    [CH3:2][...][CH:4]>>[CH3:1][...][CH3:9].[CH4:2][...][CH4:8]     #
#    CHYmap    [CH3:6][...][CH:7]>>[CH3:7][...][CH3:1].[CH3:8][...][CH3:4]     #
#    ...                                                                       #
#                                                                              #
#    * the line with the unmapped SMILES should always start with "#"          #
#    * there is no limit to the number of reactions that can be analized       #
#    * the number of maps per rxn doesn't need to be the same for all rxns     #
#                                                                              #
#    Then the program runs a "sanity check" over the three alternative         #
#    methods for evaluating the equivalence of aams: (AUX) comparing the       #
#    auxiliary graphs, (CGR or ITS) comparing the condensed graphs of the      #
#    reactions and (ISO) comparing the isomorphisms associated each reaction.  #
#    It returns a pdf with boxplots of the time (seconds) taken by each method #
#    to analyze each reaction once, another pdf with lineplots of the averag   #
#    time (seconds) taken by the methods to analyze each reaction 10 times.    #
#    Also three pickles each containig the eqivalence relation of the maps per #
#    reaction, and a txt summary indicating whether the provided maps for each #
#    reaction were equivalent or not. Note that the atom map equivalence       #
#    classes of the maps can be retrieved from the information in the PKL's.   #
#                                                                              #
#  - Input: plain-text *_aam.smiles file containing unmapped reactions-SMILES  #
#    as well as the corresponding four atom-atom maps: RXN, RDT, CHY.          #
#                                                                              #
#  - Output: (1) pdf with the mentioned time (seconds) boxplots. (2) pdf with  #
#    the mentioned time (seconds) lineplots. (3) pkl with the auxiliary graphs #
#    and equivalence classes found by AUX. (4) pkl with the condensed grpahs   #
#    of each reaction and the equivalnce classes found by the CGR method.      #
#    (5) pkl with the reactants graph G, the products graph H and the          #
#    equivalence classes found by ISO. (6) a txt file with a summary           #
#    indicating if the atom maps where equivalent or not under each method, as #
#    well as the number of equivalence classes found and the time (in seconds) #
#    taken by each method for each reaction. Since the 3 methods are           #
#    mathematically equivalent they should always return the same answer       #
#    regarding the equivalence of the 3 provided atom maps.                    #
#                                                                              #
#  - Run with:                                                                 #
#             python  4_EEquAAM_over_Symmetric.py  [myFile_aam.smiles]         #
#                                                                              #
#  - Expected output:                                                          #
#    (1)   myFile_aam_bxp.pdf                                                  #
#    (2)   myFile_aam_lnp.pdf                                                  #
#    (3)   myFile_aam_aux.pkl                                                  #
#    (4)   myFile_aam_its.pkl                                                  #
#    (5)   myFile_aam_iso.pkl                                                  #
#    (6)   myFile_aam_summary.txt                                              #
#                                                                              #
#  - Date: 12 January 2023                                                     #
#                                                                              #
################################################################################


# Requires #####################################################################
"""
> Lenguaje: python 3.9.13
> Anaconda: conda 22.9.0
> Packages installed with anaconda:
***** numpy 1.21.5
***** pysmiles 1.0.2
***** networkx 2.8.4
***** matplotlib 3.5.2
"""


# Dependencies #################################################################


# installed with conda ---------------------------------------------------------
import numpy as np
import pysmiles as ps
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from networkx.algorithms import isomorphism
from networkx.algorithms.isomorphism import generic_node_match, generic_edge_match


# already in python ------------------------------------------------------------
import time
import pickle
from sys import argv
from operator import eq
from copy import deepcopy
from math import factorial, modf, log10


# turn off warnings and matplotlib to run in the server ------------------------
import warnings
matplotlib.use("Agg")
warnings.filterwarnings("ignore")


# variables ####################################################################


# input ------------------------------------------------------------------------
inputFileName = argv[1]
inputSMILES = []
inputFile = None


# output -----------------------------------------------------------------------
outputFileNameBoxPlot = inputFileName.replace(".smiles", "_bxp.pdf")
outputFileNameLinePlot = inputFileName.replace(".smiles", "_lnp.pdf")
outputFileNamePklAUX = inputFileName.replace(".smiles", "_aux.pkl")
outputFileNamePklCGR = inputFileName.replace(".smiles", "_its.pkl")
outputFileNamePklISO = inputFileName.replace(".smiles", "_iso.pkl")
outputFileNameSummary = inputFileName.replace(".smiles", "_summary.txt")
outputFile = None
summaryOutput = ""


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
# line plots
repetitions = 10
timesReactionAUX = dict()
avgReactionAUX = []
stdReactionAUX = []
timesReactionISO = dict()
avgReactionISO = []
stdReactionISO = []
timesReactionCGR = dict()
avgReactionCGR = []
stdReactionCGR = []


# Functions ####################################################################


# function: print custom progress bar ------------------------------------------
def printProgress(casePercentage, caseNum, totCases, subCase):
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
    bar = "- progress:   0%  [" + finished + "]  100%" + " ;  case: " + str(caseNum) + " / " + str(totCases) + " ;  subcase: " + str(subCase)
    # message    
    print(bar + tail, end = "\r")

    
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
    foundEquivalent = False
    nodeMatch = None
    edgeMatch = None        
    someG = None
    someH = None
    someCGR = None    
    # obtain auxiliary graphs
    for (someID, someMap, someG, someH) in someResults:
        # creat null from copy of G in order to preserve attributes
        someCGR = deepcopy(someG)
        someCGR.remove_edges_from(list(someCGR.edges()))        
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
    nodeLabelNames = ["element", "aromatic", "hcount", "charge"]
    nodeLabelDefault = ["*", False, 0, 0]
    nodeLabelOperator = [eq, eq, eq, eq]
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
print(">>> 4_EEquAAM_over_Symmetric")


# load predefined number of mapped smiles file (*_aam.smiles file)
inputFile = open(inputFileName, "r")
inputSMILES = inputFile.read().splitlines()
inputFile.close()


# get atom-atom mappings per reaction together with their ids
currentReaction = ""
for eachLine in inputSMILES:
    eachLineTuple = tuple(eachLine.split("\t"))
    if(eachLineTuple[0] == "#"):
        currentReaction = eachLineTuple[1]
        mappedSMILES[currentReaction] = []
    else:
        mappedSMILES[currentReaction].append(eachLineTuple)


# task message
print("\n")
print("* Received " + str(len(list(mappedSMILES.keys()))) + " reaction SMILES ...")
        
        
# turn mapped smiles into networkx graphs G of reactants and H of products
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
for someL in range(len(list(graphsByMap.keys()))):
    timesReactionAUX[1 + someL] = []    
    timesReactionISO[1 + someL] = []    
    timesReactionCGR[1 + someL] = []    

    
# task message
print("\n")
print("* running AUX: comparing auxiliary graphs of each reaction ...")


# analyze auxiliary graphs for the maps of each reaction
count = 0
lengthC = 1
for eachReaction in list(graphsByMap.keys()):
    initialTime = time.time()
    resultsAUX[eachReaction] = analyzeAuxiliaryGraphs(graphsByMap[eachReaction])
    finalTime = time.time()    
    timeAUX[eachReaction] = finalTime-initialTime
    count = count + 1    
    for i in range(repetitions):
        initialTime = time.time()
        resultsAUX[eachReaction] = analyzeAuxiliaryGraphs(graphsByMap[eachReaction])
        finalTime = time.time()    
        timeAUX[eachReaction] = finalTime-initialTime
        timesReactionAUX[lengthC].append(finalTime-initialTime)
        printProgress(count*100/len(list(graphsByMap.keys())), count, len(list(graphsByMap.keys())), i+1)        
    lengthC = lengthC + 1

    
# task message
print("\n")
print("* running ITS: comparing imaginary transition state graphs of each reaction ...")

    
# analyze CGRs for the maps of each reaction
count = 0
lengthC = 1
for eachReaction in list(graphsByMap.keys()):
    initialTime = time.time()
    resultsCGR[eachReaction] = analyzeCGRs(graphsByMap[eachReaction])
    finalTime = time.time()    
    timeCGR[eachReaction] = finalTime-initialTime
    count = count + 1
    for i in range(repetitions):
        initialTime = time.time()
        resultsCGR[eachReaction] = analyzeCGRs(graphsByMap[eachReaction])
        finalTime = time.time()    
        timeCGR[eachReaction] = finalTime-initialTime
        timesReactionCGR[lengthC].append(finalTime-initialTime)
        printProgress(count*100/len(list(graphsByMap.keys())), count, len(list(graphsByMap.keys())), i+1) 
    lengthC = lengthC + 1


# task message
print("\n")
print("* running ISO: comparing isomorphisms associated to each reaction ...")
    

# analyze ISOs for the maps of each reaction
count = 0
lengthC = 1
for eachReaction in list(graphsByMap.keys()):
    initialTime = time.time()
    resultsISO[eachReaction] = analyzeISOs(graphsByMap[eachReaction])
    finalTime = time.time()    
    timeISO[eachReaction] = finalTime-initialTime
    count = count + 1    
    for i in range(repetitions):
        initialTime = time.time()
        resultsISO[eachReaction] = analyzeISOs(graphsByMap[eachReaction])
        finalTime = time.time()    
        timeISO[eachReaction] = finalTime-initialTime
        timesReactionISO[lengthC].append(finalTime-initialTime)
        printProgress(count*100/len(list(graphsByMap.keys())), count, len(list(graphsByMap.keys())), i+1)                
    lengthC = lengthC + 1
    

# obtain average and standard deviations of time data
for someL in range(len(list(graphsByMap.keys()))):
    avgReactionAUX.append(np.mean(timesReactionAUX[1 + someL]))
    stdReactionAUX.append(np.std(timesReactionAUX[1 + someL]))
    avgReactionISO.append(np.mean(timesReactionISO[1 + someL]))
    stdReactionISO.append(np.std(timesReactionISO[1 + someL]))
    avgReactionCGR.append(np.mean(timesReactionCGR[1 + someL]))
    stdReactionCGR.append(np.std(timesReactionCGR[1 + someL]))


# plot average and standard deviation
fig, ax = plt.subplots()
X = list(timesReactionAUX.keys())
ax.errorbar(X, avgReactionISO, yerr = stdReactionISO, color = "k", fmt = "-.o", label = "ISO",
             capsize = 3, elinewidth = 0.1, capthick = 0.3, linewidth = 1, markersize = 2)
ax.errorbar(X, avgReactionAUX, yerr = stdReactionAUX, color = "dimgrey", fmt = "--o", label = r"AUX $\Gamma$",
             capsize = 3, elinewidth = 0.1, capthick = 0.3, linewidth = 1, markersize = 2)
ax.errorbar(X, avgReactionCGR, yerr = stdReactionCGR, color = "grey", fmt = ":o", label = r"ITS $\Upsilon$",
             capsize = 3, elinewidth = 0.1, capthick = 0.3, linewidth = 1, markersize = 2)
ax.legend()
# logAvTimeISO = [log10(t) for t in avgReactionISO]
# logAvTimeAUX = [log10(t) for t in avgReactionAUX]
# logAvTimeCGR = [log10(t) for t in avgReactionCGR]
# plt.plot(X, logAvTimeISO, color = "k", linestyle = "-.", marker = "o", label = "ISO", linewidth = 1, markersize = 2)
# plt.plot(X, logAvTimeAUX, color = "dimgrey", linestyle = "--", marker = "o", label = r"AUX $\Gamma$", linewidth = 1, markersize = 2)
# plt.plot(X, logAvTimeCGR, color = "grey", linestyle = ":", marker = "o", label = r"ITS $\Upsilon$", linewidth = 1, markersize = 2)
# plt.ylabel("$Log_{10}$ of average running time (measured in seconds) \n over " + str(repetitions) + " iterations", fontsize = 9)
# plt.legend()
plt.yticks(fontsize = 8)
plt.xticks(ticks = X, labels = ["Step 1", "Step 2", "Step 3", "Step 4"], fontsize = 9, rotation = 45)
plt.ylabel("Average running time (measured in seconds) \n over " + str(repetitions) + " iterations", fontsize = 9)
plt.grid(color = "lightgray", axis = "y", linestyle = "--", linewidth = 0.7)
plt.tight_layout()
plt.savefig(outputFileNameLinePlot)
plt.close()


# task message
print("* making plots ...")

    
# define data to plot
logTimeCGR = [log10(t) for t in list(timeCGR.values())]
logTimeISO = [log10(t) for t in list(timeISO.values())]
logTimeAUX = [log10(t) for t in list(timeAUX.values())]
timeData = [logTimeISO, logTimeAUX, logTimeCGR]
timeLabels = ["ISO", r"AUX $\Gamma$", r"ITS $\Upsilon$"]
whiskers = 1.25
widthsBoxes = 0.4
hatchB = ["....", "////", "oo"]
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
plt.ylabel(r"$Log_{10}$ of running time (measured in seconds)", fontsize = 10)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 10)
plt.title("Analyzed " + str(len(list(graphsByMap.keys()))) + " reactions", fontsize = 9)
plt.grid(color = "lightgray", axis = "y", linestyle = "--", linewidth = 0.7)
plt.tight_layout()
plt.savefig(outputFileNameBoxPlot)
plt.close()


# task message
print("* making summary ...")


# save summary in txt file
dataStr = ""
summaryOutput = ""
# add lines per reaction
for eachReaction in list(graphsByMap.keys()):
    # add reaction smiles
    summaryOutput = summaryOutput + "#" + "\t" + eachReaction + "\n"
    # add number of input maps for reaction
    dataStr = str(len(graphsByMap[eachReaction]))
    summaryOutput = summaryOutput + "maps given" + "\t" + dataStr + "\n"
    # add number of classes obtained when comparing auxiliary graphs
    dataStr = str(max([classR for (classR, idR, mapR, GR, HR, auxR) in resultsAUX[eachReaction]]))
    summaryOutput = summaryOutput + "classes AUX" + "\t" + dataStr + "\n"    
    # add number of classes obtained when comparing CGRs
    dataStr = str(max([classR for (classR, idR, mapR, GR, HR, cgrR) in resultsCGR[eachReaction]]))
    summaryOutput = summaryOutput + "classes ITS" + "\t" + dataStr + "\n"
    # add number of classes obtained when comparing ISOs
    dataStr = str(max([classR for (classR, idR, mapR, GR, HR) in resultsISO[eachReaction]]))
    summaryOutput = summaryOutput + "classes ISO" + "\t" + dataStr + "\n"    
    # add result AUX
    if(max([classR for (classR, idR, mapR, GR, HR, auxR) in resultsAUX[eachReaction]]) == 1):
        dataStr = "equivalent maps"
    else:
        dataStr = "not equivalent maps"        
    summaryOutput = summaryOutput + "result AUX" + "\t" + dataStr + "\n"
    # add result CGR
    if(max([classR for (classR, idR, mapR, GR, HR, cgrR) in resultsCGR[eachReaction]]) == 1):
        dataStr = "equivalent maps"
    else:
        dataStr = "not equivalent maps"        
    summaryOutput = summaryOutput + "result ITS" + "\t" + dataStr + "\n"
    # add result ISO
    if(max([classR for (classR, idR, mapR, GR, HR) in resultsISO[eachReaction]]) == 1):
        dataStr = "equivalent maps"
    else:
        dataStr = "not equivalent maps"
    summaryOutput = summaryOutput + "result ISO" + "\t" + dataStr + "\n"    
    # add time AUX        
    dataStr = str(timeAUX[eachReaction])
    summaryOutput = summaryOutput + "time AUX" + "\t" + dataStr + "\n"
    # add time CGR        
    dataStr = str(timeCGR[eachReaction])
    summaryOutput = summaryOutput + "time ITS" + "\t" + dataStr + "\n"
    # add time ISO        
    dataStr = str(timeISO[eachReaction])
    summaryOutput = summaryOutput + "time ISO" + "\t" + dataStr + "\n"        
outputFile = open(outputFileNameSummary, "w")
outputFile.writelines(summaryOutput)
outputFile.close()


# task message
print("* saving data into pkl files ...")


# save pkl data of AUX
outputFile = open(outputFileNamePklAUX, "wb")
pickle.dump(resultsAUX, outputFile)
outputFile.close()


# save pkl data of CGR
outputFile = open(outputFileNamePklCGR, "wb")
pickle.dump(resultsCGR, outputFile)
outputFile.close()


# save pkl data of ISO
outputFile = open(outputFileNamePklISO, "wb")
pickle.dump(resultsISO, outputFile)
outputFile.close()


# final message
print("\n")
print(">>> Finished!!!")
print("\n")


# End ##########################################################################
################################################################################
