################################################################################
#                                                                              #
#  README - Program: 2_Symmetries.py                                           #
#                                                                              #
#  - Paper: [tbd]                                                              #
#                                                                              #
#  - Github: [tbd]                                                             #
#                                                                              #
#  - Contributor(s):                                                           #
#  * @MarcosLaffitte - Marcos E. González Laffitte                             #
#                                                                              #
#  - Description: this script receives a list of reaction-SMILES and           #
#    for each reaction G>>H it determines the number of symmetries of the      #
#    simple graph G whose connected components are the reactants of G>>H,      #
#    as well as the symmetries of the simple graph H whose connected           #
#    components are the products of the reaction G>>H. Here symmetries are     #
#    to be understood as node-label and edge-label preserving automorphisms    #
#    of G and H, i.e., all possible "combinatorial symmetries" of such graphs, #
#    which includes the more usual "geometric symmetries" of molecules.        #
#                                                                              #
#  - Input: plain-text file with *.smiles extension whose lines are each a     #
#    single reaction-SMILES (mapped or unmapped).                              #
#                                                                              #
#  - Output: (1) plain-text file with *.txt extension containing the number of #
#    symmetries of the reactants, number of symmetries of products and the     #
#    original reaction smiles corresponding to each given reaction, these 3    #
#    values separated by tabs. (2) *.pdf with barplots showing the number of   #
#    reactions whose reactants/products had a specific number of symmetries.   #
#                                                                              #
#  - Run with:                                                                 #
#                    python  2_Symmetries.py  [myFile.smiles]                  #
#                                                                              #
#  - Expected output:                                                          #
#    (1)   myFile_summary.txt                                                  #
#    (2)   myFile_barplot.pdf                                                  #
#                                                                              #
#  - Date: 12 January 2023                                                     #
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
from sys import argv
from operator import eq
from copy import deepcopy
from math import factorial, modf


# turn off warnings and matplotlib to run in a server --------------------------
import warnings
matplotlib.use("Agg")
warnings.filterwarnings("ignore")


# variables ####################################################################


# input ------------------------------------------------------------------------
inputSMILES = []
inputFileName = argv[1]
inputFile = None


# output -----------------------------------------------------------------------
outputFileSummary = inputFileName.replace(".smiles", "_summary.txt")
outputFileBarplot = inputFileName.replace(".smiles", "_barplot.pdf")
outputFile = None
newLine = ""
theSummary = ""
summaryLines = []


# data -------------------------------------------------------------------------
percentage = 0
eachReaction = ""
productsSMILES = ""
reactantsSMILES = ""
G = None
H = None
resultFormat = ()
results = []
sRvals = []
sPvals = []
XR = []
XP = []
tempList = []


# Functions ####################################################################


# function: print custom progress bar ------------------------------------------
def printProgress(casePercentage, caseNumber, totCases):
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
    bar = "- progress:   0%  [" + finished + "]  100%"
    # add case number to bar
    bar = bar + "  ;  done:  " + str(caseNumber) + " / " + str(totCases)
    # message    
    print(bar + tail, end = "\r")


# function: obtain graph from molecules-SMILES ---------------------------------
def moleculesIntoGraph(someMolecules):
    # local variables
    eachMolecule = ""
    moleculeList = []
    theGraph = nx.Graph()
    newMolecule = nx.Graph()            
    # get graph whose connected components are molecules
    moleculeList = someMolecules.split(".")        
    for eachMolecule in moleculeList:
        # read smiles into networkx graph
        newMolecule = ps.read_smiles(eachMolecule)
        # make disjoint union over theGraph (presrving smiles attributes)
        theGraph = nx.disjoint_union(theGraph, newMolecule)
    # end of function
    return(theGraph)


# function: obtain label-preserving automorphisms ------------------------------
def obtainSymmetries(someGraph):
    # local variables
    totSymmetries = 0
    nodeMatch = None
    edgeMatch = None
    nodeLabelNames = []
    nodeLabelDefault = []
    nodeLabelOperator = []
    graphMatcher = None
    # determine labels to preserv during automorphism check (as defined by pysmiles)
    # * "class" usually denotes atom-to-atom map, but can be accesed as follows if needed:
    # * nodeLabelNames = ["element", "aromatic", "hcount", "charge", "class"]
    # * nodeLabelDefault = ["*", False, 0, 0, 0]
    # * nodeLabelOperator = [eq, eq, eq, eq, eq]    
    nodeLabelNames = ["element", "aromatic", "hcount", "charge"]
    nodeLabelDefault = ["*", False, 0, 0]
    nodeLabelOperator = [eq, eq, eq, eq]
    nodeMatch = generic_node_match(nodeLabelNames,
                                   nodeLabelDefault,
                                   nodeLabelOperator)
    edgeMatch = generic_edge_match("order", 1, eq)
    # obtain number of symmetries
    graphMatcher = isomorphism.GraphMatcher(someGraph, someGraph,
                                            node_match = nodeMatch,
                                            edge_match = edgeMatch)
    totSymmetries = len(list(graphMatcher.isomorphisms_iter()))
    # end of function
    return(totSymmetries)


# Main #########################################################################


# initial message
print("\n")
print(">>> 2_Symmetries")
print("\n")


# task message
print("* accessing input file ...")


# open input file
inputFile = open(inputFileName, "r")
inputSMILES = inputFile.read().splitlines()
inputFile.close()


# task message
print("* turning SMILES into graphs and obtaining symmetries ...")


# turn reaction-SMILES into graphs and obtain symmetries
for i in range(len(inputSMILES)):
    eachReaction = inputSMILES[i]
    reactantsSMILES = eachReaction.split(">>")[0]
    productsSMILES = eachReaction.split(">>")[1]
    # obtain reactants graph
    G = moleculesIntoGraph(reactantsSMILES)
    # obtain products graph
    H = moleculesIntoGraph(productsSMILES)
    # obtain symmetries of reactants
    symmsG = obtainSymmetries(G)
    # obtain symmetries of products
    symmsH = obtainSymmetries(H)
    # save results
    resultFormat = (symmsG, symmsH, eachReaction)
    results.append(resultFormat)
    # print progress
    percentage = round((i+1) * 100 / len(inputSMILES), 2)
    printProgress(percentage, i+1, len(inputSMILES))
print("\n")


# task message
print("* making txt summary file ...")


# make lines for summary file
for (sR, sP, eachReaction) in results:
    newLine = str(sR) + "\t" + str(sP) + "\t" + eachReaction + "\n"
    summaryLines.append(newLine)


# print summary file    
theSummary = "".join(summaryLines)
outputFile = open(outputFileSummary, "w")
outputFile.write(theSummary)
outputFile.close()


# task message
print("* making pdf barplot file ...")


# get proportion of reactions with each value of symmetries and make barplot
sRvals = [sR for (sR, sP, eachReaction) in results]
sPvals = [sP for (sR, sP, eachReaction) in results]
# get X positions and H heights of data
XR = list(set(sRvals))
XP = list(set(sPvals))
XR.sort()
XP.sort()
HR = []
HP = []
Xticks = list(set(XR + XP))
Xticks.sort()
# get heights indexed by values of symmetries
for i in range(len(Xticks)):
    tempList = [val for val in sRvals if val == Xticks[i]]
    if(len(tempList) == 0):
        HR.append(0)
    else:
        HR.append(len(tempList) * 100 / len(sRvals))
for i in range(len(Xticks)):
    tempList = [val for val in sPvals if val == Xticks[i]]
    if(len(tempList) == 0):
        HP.append(0)
    else:
        HP.append(len(tempList) * 100 / len(sPvals))
# mak bar plots
positionsR = [val - 0.15 for val in range(len(Xticks))]
positionsP = [val + 0.15 for val in range(len(Xticks))]
plt.bar(x = positionsR, height = HR, width = 0.3, color = "w", edgecolor = "darkgrey", zorder = 3, label = "Reactants Graph", hatch = "......")
plt.bar(x = positionsP, height = HP, width = 0.3, color = "w", edgecolor = "dimgrey", zorder = 3, label = "Products Graph", hatch = "//////")
# change image attributes
plt.legend(fontsize = 8)
plt.grid(visible = True, color = "lightgray", linestyle = "--", linewidth = 0.4)
plt.xticks(ticks = range(len(Xticks)), labels = Xticks, fontsize = 7, rotation = 45)
plt.yticks(fontsize = 8)
plt.xlabel("Number of symmetries computationally detected", fontsize = 10)
plt.ylabel("Percentage of reactions\n", fontsize = 10)
plt.title("Analyzed " + str(len(inputSMILES)) + " reactions", fontsize = 9)
plt.tight_layout()
plt.savefig(outputFileBarplot)
plt.close()


# final message
print("\n")
print(">>> Finished")
print("\n")


# End ##########################################################################
################################################################################
