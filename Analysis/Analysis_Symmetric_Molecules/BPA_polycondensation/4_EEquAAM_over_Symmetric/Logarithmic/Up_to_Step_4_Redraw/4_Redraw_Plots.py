################################################################################
#                                                                              #
# - Redraw some plots                                                          #
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
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf


# already in python ------------------------------------------------------------
import pickle
from sys import argv
from math import log10


# turn off warnings and matplotlib to run in the server ------------------------
import warnings
matplotlib.use("Agg")
warnings.filterwarnings("ignore")


# Main #########################################################################


# initial message
print("\n")
print(">>> Redraw")


# save original data
repetitions  = 5
inputFile = open(argv[1], "rb")
[timesReactionAUX, timesReactionISO, timesReactionCGR] = pickle.load(inputFile)
inputFile.close()


# transform log scale and obtain average and std of time data
avgReactionAUX = []
stdReactionAUX = []
avgReactionISO = []
stdReactionISO = []
avgReactionCGR = []
stdReactionCGR = []
for someL in [0, 1, 2, 3]:
    timesReactionAUX[1 + someL] = [log10(t) for t in timesReactionAUX[1 + someL]]
    timesReactionISO[1 + someL] = [log10(t) for t in timesReactionISO[1 + someL]]
    timesReactionCGR[1 + someL] = [log10(t) for t in timesReactionCGR[1 + someL]]
for someL in [0, 1, 2, 3]:
    avgReactionAUX.append(np.mean(timesReactionAUX[1 + someL]))
    stdReactionAUX.append(np.std(timesReactionAUX[1 + someL]))
    avgReactionISO.append(np.mean(timesReactionISO[1 + someL]))
    stdReactionISO.append(np.std(timesReactionISO[1 + someL]))
    avgReactionCGR.append(np.mean(timesReactionCGR[1 + someL]))
    stdReactionCGR.append(np.std(timesReactionCGR[1 + someL]))


# task message
print("\n")
print("* making plots ...")


# plot average and standard deviation
fig, ax = plt.subplots()
X = list(timesReactionAUX.keys())
ax.errorbar(X, avgReactionISO, yerr = stdReactionISO, color = "k", fmt = "-.o", label = r"ISO-$\equiv$",
             capsize = 3, elinewidth = 0.1, capthick = 0.3, linewidth = 1, markersize = 2)
ax.errorbar(X, avgReactionAUX, yerr = stdReactionAUX, color = "dimgrey", fmt = "--o", label = r"AUX-$\Gamma$",
             capsize = 3, elinewidth = 0.1, capthick = 0.3, linewidth = 1, markersize = 2)
ax.errorbar(X, avgReactionCGR, yerr = stdReactionCGR, color = "grey", fmt = ":o", label = r"ITS-$\Upsilon$",
             capsize = 3, elinewidth = 0.1, capthick = 0.3, linewidth = 1, markersize = 2)
ax.legend(fontsize = 13)
plt.yticks(fontsize = 13)
plt.xticks(ticks = X, labels = ["Step 1", "Step 2", "Step 3", "Step 4"], fontsize = 12, rotation = 45)
plt.ylabel("Average of  " + r"$Log_{10}$" + "  of running time [s]\n" + "over " + str(repetitions) + " iterations\n", fontsize = 14)
plt.grid(color = "lightgray", axis = "y", linestyle = "--", linewidth = 0.7)
plt.tight_layout()
plt.savefig("real_symmetric_suitable_aam_bxp.pdf")
plt.close()

    
# define data to plot
logTimeISO = []
logTimeAUX = []
logTimeCGR = []
for someL in [0, 1, 2, 3]:
    logTimeISO = logTimeISO + timesReactionISO[1 + someL]
    logTimeAUX = logTimeAUX + timesReactionAUX[1 + someL]
    logTimeCGR = logTimeCGR + timesReactionCGR[1 + someL]
timeData = [logTimeISO, logTimeAUX, logTimeCGR]
timeLabels = [r"ISO-$\equiv$", r"AUX-$\Gamma$", r"ITS-$\Upsilon$"]
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
plt.ylabel(r"$Log_{10}$" + " of running times [s]\n" + "over " + str(repetitions) + " iterations\n", fontsize = 16)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 14)
plt.title("Analyzed 4 reactions", fontsize = 14)
plt.grid(color = "lightgray", axis = "y", linestyle = "--", linewidth = 0.7)
plt.tight_layout()
plt.savefig("real_symmetric_suitable_aam_lnp.pdf")
plt.close()


# final message
print("\n")
print(">>> Finished")
print("\n")


# End ##########################################################################
################################################################################
