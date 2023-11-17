################################################################################
#                                                                              #
#  > 5_Comparison_of_AAMs                                                      #
#                                                                              #
#  > Run:   python  5_Comparison_of_AAMs.py   [Reactions_aam.smiles]           #
#                                                                              #
#  > Description: receives a tab-separated txt (\t) containing both the simple #
#    and mapped SMILES of a set of reactions, and ids for them as follows:     #
#                                                                              #
#    * e.g. Reactions_aam.smiles                                               #
#                                                                              #
#    #            CCCCCCCCCC>>CCCCCCCCC.C                                      #
#    id_map_1     [CH3:1]...[CH:10]>>[CH3:1]...[CH3:9].[CH4:10]                #
#    id_map_2     [CH3:2]...[CH:4]>>[CH3:1]...[CH3:9].[CH4:2]                  #
#    id_map_3     [CH3:6]...[CH:7]>>[CH3:7]...[CH3:1].[CH4:6]                  #
#    #            CCCCCCCCCC>>CCCC.CCCCCC                                      #
#    id_map_1     [CH3:1]...[CH:10]>>[CH3:1]...[CH3:4].[CH3:5]...[CH4:10]      #
#    id_map_2     [CH3:2]...[CH:4]>>[CH3:1]...[CH3:9].[CH4:2]...[CH4:8]        #
#    id_map_3     [CH3:6]...[CH:7]>>[CH3:7]...[CH3:1].[CH3:8]...[CH3:4]        #
#    ...                                                                       #
#                                                                              #
#    * the line with the simple SMILES should start with "#"                   #
#    * there is no limit to the number of reactions that can be analized       #
#    * the number of maps per rxn doesn't need to be the same for all rxns     #
#                                                                              #
#    Right now the program runs a "sanity check" over the 3 alternative        #
#    methods for evaluating the equivalence aams: (AUX) comparin the auxiliary #
#    graphs, (CGR) comparing the condensed graphs of the reactions and         #
#    (ISO) comparing the isomorphisms associated to the reactions.             #
#    It returns: a pdf with the box plots of the time taken by each method,    #
#    3 pickles each containig the eqivalence relation of the maps per rxn,     #
#    and a txt summary indicating whether the provided maps were equivalent.   #
#                                                                              #
#    * the atom map equivalence classes can be retrieved from the PKL's        #
#                                                                              #
################################################################################


# WARNING: pysmiles may do weird thigs with aromatic rings and the hcount!!!


# Requires #####################################################################
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
from sys import argv
from operator import eq
from copy import deepcopy
from math import factorial, modf, log10
import itertools

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
outputFileNameBoxPlot = inputFileName.replace(".smiles", "_bxps.pdf")
outputFileNamePklAUX = inputFileName.replace(".smiles", "_auxr.pkl")
outputFileNamePklCGR = inputFileName.replace(".smiles", "_cgrr.pkl")
outputFileNamePklISO = inputFileName.replace(".smiles", "_isor.pkl")
outputFileNameSummary = inputFileName.replace(".smiles", "_smry.txt")
outputFile = None
summaryOutput = ""


# data -------------------------------------------------------------------------
# properties of pysmiles atoms e.g.
# vertex v : data = {'charge': 0, 'hcount': 2, 'aromatic': False, 'element': C, N or O, 'class': aamNum}
# properties of pysmiles bonds e.g.
# edge (u, v) : data = {'order': bond-type 1 (-), 2 (=) or 3 (#)}
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


# Functions ####################################################################


# function: print custom progress bar ------------------------------------------
def printProgress(casePercentage):
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
    bar = "- Progress:   0%  [" + finished + "]  100%"
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
    nodeMatch = generic_node_match(["GorH", "element", "hcount"], [0, "C", 0], [eq, eq, eq])
    edgeMatch = generic_edge_match("order", "", eq)
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
    nodeMatch = generic_node_match(["element", "hcount"], ["C", 0], [eq, eq])
    edgeMatch = generic_edge_match("order", "", eq)
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
    nodeMatch = generic_node_match(["element", "hcount"], ["C", 0], [eq, eq])
    edgeMatch = generic_edge_match("order", "", eq)
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
print(">>> 5_Comparison_of_AAMs")

def isa_group_separator(line):
    return line=='\n'


with open(inputFileName, mode='r') as rf:
  for key,group in itertools.groupby(rf,isa_group_separator):
    current_maps = list(group)
    #print(key)
    print(list(group))
    print(current_maps)
    if key:
         print(key)
         continue
    print(len(list(group)))
    currentReaction, smiles, rxn_map, graphormer_map, rdt_map = map(str.strip, current_maps)
    #currentReaction, smiles, rxn_map, graphormer_map, rdt_map = list(group)
    print(currentReaction)
    print(smiles)
    print('RXN MAPPER:')
    print(rxn_map)
    print('Graphormer Map:')
    print(graphormer_map)
    print('RDT Map:')
    print(rdt_map) 
    # working message
    print("\n")
    print("* Received " + str(len(list(mappedSMILES.keys()))) + " reaction SMILES ...")
        
    # turn mapped smiles into networkx graphs G of reactants and H of products
    tempResults = []

    # split reaction
    reactantsSide_rxn = (rxn_map.split(">>"))[0]
    productsSide_rxn = (rxn_map.split(">>"))[1]
    
    reactantsSide_graphormer = (graphormer_map.split(">>"))[0]
    productsSide_graphormer = (graphormer_map.split(">>"))[1]
    
    reactantsSide_rdt = (rdt_map.split(">>"))[0]
    productsSide_rdt = (rdt_map.split(">>"))[1]
    
    # get reactants graph G
    reactantsList_rxn = reactantsSide_rxn.split(".")    
    reactantsList_graphormer = reactantsSide_graphormer.split(".")     
    reactantsList_rdt = reactantsSide_rdt.split(".")     
        
    G_rxn = nx.Graph()
    G_graphormer = nx.Graph()
    G_rdt = nx.Graph()
    for eachReactant in range(len(reactantsList_rxn)):
            # read smiles into networkx graph            
            molecule_rxn = ps.read_smiles(reactantsList_rxn[eachReactant])
            molecule_graphormer = ps.read_smiles(reactantsList_graphormer[eachReactant])
            molecule_rdt = ps.read_smiles(reactantsList_rdt[eachReactant])
            # rename vertices using atom map            
            atomMapping_rxn = nx.get_node_attributes(molecule_rxn, "class")
            atomMapping_graphormer = nx.get_node_attributes(molecule_graphormer, "class")
            atomMapping_rdt = nx.get_node_attributes(molecule_rdt, "class")
            molecule_rxn = nx.relabel_nodes(molecule_rxn, atomMapping_rxn)
            molecule_graphormer = nx.relabel_nodes(molecule_graphormer, atomMapping_graphormer)
            molecule_rdt = nx.relabel_nodes(molecule_rdt, atomMapping_rxn)
            # make (disjoint) union over G          
            try:  
                G_rxn = deepcopy(nx.union(G_rxn, molecule_rxn))
                G_graphormer = deepcopy(nx.union(G_graphormer, molecule_graphormer))
                G_rdt = deepcopy(nx.union(G_rdt, molecule_rdt))
            except:
                print('UNION PROBLEMS!')
                continue
    # get products graph H
    productsList_rxn = productsSide_rxn.split(".")     
    productsList_graphormer = productsSide_graphormer.split(".")        
    productsList_rdt = productsSide_rdt.split(".")           
    H_rxn = nx.Graph()
    H_graphormer = nx.Graph()
    H_rdt = nx.Graph()
    for eachProduct in range(len(productsList_rxn)):
            # read smiles into networkx graph            
            molecule_rxn = ps.read_smiles(productsList_rxn[eachProduct])
            molecule_graphormer = ps.read_smiles(productsList_graphormer[eachProduct])
            molecule_rdt = ps.read_smiles(productsList_rdt[eachProduct])
            # rename vertices using atom map            
            atomMapping_rxn = nx.get_node_attributes(molecule_rxn, "class")
            atomMapping_graphormer = nx.get_node_attributes(molecule_graphormer, "class")
            atomMapping_rdt = nx.get_node_attributes(molecule_rdt, "class")
            molecule_rxn = nx.relabel_nodes(molecule_rxn, atomMapping_rxn)
            molecule_graphormer = nx.relabel_nodes(molecule_graphormer, atomMapping_graphormer)
            molecule_rdt = nx.relabel_nodes(molecule_rdt, atomMapping_rxn)
            # make (disjoint) union over H            
            try:    
                H_rxn = deepcopy(nx.union(H_rxn, molecule_rxn))
                H_graphormer = deepcopy(nx.union(H_graphormer, molecule_graphormer))
                H_rdt = deepcopy(nx.union(H_rdt, molecule_rdt))
            except:
                print('UNION PROBLEMS!')      
    # save map data
    #print(rxn_map)
    tempResults.append(('RXNmap', rxn_map , deepcopy(G_rxn), deepcopy(H_rxn)))
    tempResults.append(('Graphormer', graphormer_map , deepcopy(G_graphormer), deepcopy(H_graphormer)))
    tempResults.append(('RDT', rdt_map , deepcopy(G_rdt), deepcopy(H_rdt)))

    # save reaction data
    graphsByMap[currentReaction] = deepcopy(tempResults)
    #print(graphsByMap)
    tempResults = []

    # working message
    print("* Comparing Auxiliary Graphs of each reaction ...")
    # analyze auxiliary graphs for the maps of each reaction
    for eachReaction in list(graphsByMap.keys()):
        initialTime = time.time()
        resultsAUX[eachReaction] = analyzeAuxiliaryGraphs(graphsByMap[eachReaction])
        finalTime = time.time()    
        timeAUX[eachReaction] = finalTime-initialTime

    # working message
    print("* Comparing Condensed Graphs of each reaction ...")   
    # analyze CGRs for the maps of each reaction
    for eachReaction in list(graphsByMap.keys()):
        initialTime = time.time()
        resultsCGR[eachReaction] = analyzeCGRs(graphsByMap[eachReaction])
        finalTime = time.time()    
        timeCGR[eachReaction] = finalTime-initialTime

    # working message
    print("* Comparing Isomorphisms associated to each reaction ...")
    # analyze ISOs for the maps of each reaction
    for eachReaction in list(graphsByMap.keys()):
        initialTime = time.time()
        resultsISO[eachReaction] = analyzeISOs(graphsByMap[eachReaction])
        finalTime = time.time()    
        timeISO[eachReaction] = finalTime-initialTime
     
# working message
print("* Making plots ...")

    
# define data to plot
logTimeCGR = [log10(t) for t in list(timeCGR.values())]
logTimeISO = [log10(t) for t in list(timeISO.values())]
logTimeAUX = [log10(t) for t in list(timeAUX.values())]
timeData = [logTimeCGR, logTimeISO, logTimeAUX]
timeLabels = ["CGR", "ISO", "AUX"]
whiskers = 1.25
widthsBoxes = 0.4
colorsB = ["steelblue", "tab:green", "tab:red"]
colorsW = ["tab:blue", "tab:blue", "darkgreen", "darkgreen", "firebrick", "firebrick"]
colorsL = ["tab:blue", "darkgreen", "firebrick"]
hatchB = ["...", "///", "oo"]
# hatchB = ["", "", ""]
# hatchB = [".....", ".....", "....."]
# hatchB = [".....", "/////", "\\\\\\\\\\"]
# hatchB = ["//", "////", "//////"]
# make plot 
fig, ax = plt.subplots()
bps = ax.boxplot(timeData, whis = whiskers, widths = widthsBoxes, labels = timeLabels, sym = ".", patch_artist = True)
# set colors
for i in range(len(bps["fliers"])):
    bps["fliers"][i].set(markeredgecolor = "grey", markerfacecolor = "grey")
    # bps["fliers"][i].set(markeredgecolor = colorsL[i], markerfacecolor = colorsB[i])
for i in range(len(bps["caps"])):
    bps["caps"][i].set(color = "dimgrey")
    # bps["caps"][i].set(color = colorsW[i])    
for i in range(len(bps["boxes"])):
    bps["boxes"][i].set(facecolor = "w", edgecolor = "grey", hatch = hatchB[i])
    # bps["boxes"][i].set(facecolor = colorsB[i], edgecolor = colorsL[i], alpha = 0.5)
for i in range(len(bps["medians"])):    
    bps["medians"][i].set(color = "dimgrey")
    # bps["medians"][i].set(color = colorsL[i])        
for i in range(len(bps["whiskers"])):
    bps["whiskers"][i].set(color = "dimgrey")
    # bps["whiskers"][i].set(color = colorsW[i])
# set some attributes and save figures
plt.ylabel("$Log_{10}$ of Running Time\n")
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.title("Analyzed " + str(len(list(graphsByMap.keys()))) + " reactions", fontsize = 8)
plt.grid(color = "lightgray", axis = "y", linestyle = "--", linewidth = 0.7)
plt.tight_layout()
plt.savefig(outputFileNameBoxPlot)
plt.close()


# working message
print("* Making txt summary ...")


# save summary in txt file
dataStr = ""
summaryOutput = ""
# add lines per reaction
for eachReaction in list(graphsByMap.keys()):
    # add reaction smiles
    summaryOutput = summaryOutput + "#" + "\t" + eachReaction + "\n"
    # add number of input maps for reaction
    dataStr = str(len(graphsByMap[eachReaction]))
    summaryOutput = summaryOutput + "Maps Given" + "\t" + dataStr + "\n"
    # add number of classes obtained when comparing auxiliary graphs
    dataStr = str(max([classR for (classR, idR, mapR, GR, HR, auxR) in resultsAUX[eachReaction]]))
    summaryOutput = summaryOutput + "Classes AUX" + "\t" + dataStr + "\n"    
    # add number of classes obtained when comparing CGRs
    dataStr = str(max([classR for (classR, idR, mapR, GR, HR, cgrR) in resultsCGR[eachReaction]]))
    summaryOutput = summaryOutput + "Classes CGR" + "\t" + dataStr + "\n"
    # add number of classes obtained when comparing ISOs
    dataStr = str(max([classR for (classR, idR, mapR, GR, HR) in resultsISO[eachReaction]]))
    summaryOutput = summaryOutput + "Classes ISO" + "\t" + dataStr + "\n"    
    # add result AUX
    if(max([classR for (classR, idR, mapR, GR, HR, auxR) in resultsAUX[eachReaction]]) == 1):
        dataStr = "equivalent maps"
    else:
        dataStr = "not equivalent maps"        
    summaryOutput = summaryOutput + "Result AUX" + "\t" + dataStr + "\n"
    # add result CGR
    if(max([classR for (classR, idR, mapR, GR, HR, cgrR) in resultsCGR[eachReaction]]) == 1):
        dataStr = "equivalent maps"
    else:
        dataStr = "not equivalent maps"        
    summaryOutput = summaryOutput + "Result CGR" + "\t" + dataStr + "\n"
    # add result ISO
    if(max([classR for (classR, idR, mapR, GR, HR) in resultsISO[eachReaction]]) == 1):
        dataStr = "equivalent maps"
    else:
        dataStr = "not equivalent maps"
    summaryOutput = summaryOutput + "Result ISO" + "\t" + dataStr + "\n"    
    # add time AUX        
    dataStr = str(timeAUX[eachReaction])
    summaryOutput = summaryOutput + "Time AUX" + "\t" + dataStr + "\n"
    # add time CGR        
    dataStr = str(timeCGR[eachReaction])
    summaryOutput = summaryOutput + "Time CGR" + "\t" + dataStr + "\n"
    # add time ISO        
    dataStr = str(timeISO[eachReaction])
    summaryOutput = summaryOutput + "Time ISO" + "\t" + dataStr + "\n"        
outputFile = open(outputFileNameSummary, "w")
outputFile.writelines(summaryOutput)
outputFile.close()


# working message
print("* Saving data in PKL's ...")


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
