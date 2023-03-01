- Description: this script receives a plain-text file with *_aam.smiles format containing COMPLETELY mapped reaction SMILES and identifiers for them, as shown in the following toy-example:

#,reaction_1
RXNmap,[CH3:1][...][CH:8]>>[CH3:1][...][CH3:9].[CH4:8]                    
RDTmap,[CH3:2][...][CH:4]>>[CH3:1][...][CH3:9].[CH4:2]                    
CHYmap,[CH3:6][...][CH:7]>>[CH3:7][...][CH3:1].[CH4:6]                    
#,reaction_2                                                              
RXNmap,[CH3:1][...][CH:8]>>[CH3:1][...][CH3:4].[CH3:5][...][CH4:8]        
RDTmap,[CH3:2][...][CH:4]>>[CH3:1][...][CH3:9].[CH4:2][...][CH4:8]        
CHYmap,[CH3:6][...][CH:7]>>[CH3:7][...][CH3:1].[CH3:8][...][CH3:4]        
    ...                                                                  
                                                                              
- Input: plain-text csv *_aam.smiles file as specified above.

- Output: (1) pdf with running time of the methods. (2) pkl with the auxiliary graphs and equivalence classes found by AUX. (3) pkl with the condensed graphs of the reaction and the equivalence classes found by the CGR method. (4) pkl with the reactants graph G, the products graph H and the equivalence classes found by ISO. (5) a txt file with a summary indicating if the atom maps where equivalent or not under each method, as well as the number of equivalence classes found and the time (seconds) taken by each method for each reaction. Note that the equivalence classes of the atom maps can be retrieved from the PKL's. The time boxplots show Log_10 of time; linear time is always included in the summary.            

- Expected output:                                                          
    (1) myFile_aam_times.pdf                                                  
    (2) myFile_aam_aux.pkl        (only with --sanity-check option)             
    (3) myFile_aam_its.pkl                                                    
    (4) myFile_aam_iso.pkl        (only with --sanity-check option)             
    (5) myFile_aam_summary.txt    (tab separated file, may be changed in-code)  