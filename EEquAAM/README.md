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
#      * default ITS   python EEquAAM.py [myFile_aam.smiles]                   #
#      * or            python EEquAAM.py --sanity-check [myFile_aam.smiles]    #
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