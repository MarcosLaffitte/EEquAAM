- How to run the toy example (after activating eequaam conda environment)<br/>

```
python  EEquAAM.py  ToyExample_id_aam.smiles
```

- General Input Format: plain-text file [*_aam.smiles] with the original identifiers and the mapped SMILES as follows (this is also the output of MappingTool.py)<br/>

```
#,reaction_1
RXNmap,[CH3:1][...][CH:10]>>[CH3:1][...][CH3:9].[CH4:10]
RDTmap,[CH3:2][...][CH:4]>>[CH3:1][...][CH3:9].[CH4:2]
CHYmap,[CH3:6][...][CH:7]>>[CH3:7][...][CH3:1].[CH4:6]
#,reaction_2
RXNmap,[CH3:1][...][CH:10]>>[CH3:1][...][CH3:4].[CH3:5][...][CH4:10]
RDTmap,[CH3:2][...][CH:4]>>[CH3:1][...][CH3:9].[CH4:2][...][CH4:8]
CHYmap,[CH3:6][...][CH:7]>>[CH3:7][...][CH3:1].[CH3:8][...][CH3:4]
...
```

- Output:<br/>

a) [*_aam_summary.txt]: total of reactions with equivalent and non-equivalent mappings.<br/>

b) [*_aam_times.pdf]: PDF boxplots of the time taken by each method to compare the provided atom maps.<br/>

c) [*_aam_out_alleq.smiles]: reactions with all equivalent maps, printed out in the same format as the input file.<br/>

d) [*_aam_out_noneq.smiles]: reactions with at least one pair of non-equivalent maps, printed out in the same format as the input file, but with the following additional information:<br/>

```
#,reaction_1
(1),RXNmap,[CH3:1][...][CH:8]>>[CH3:1][...][CH3:9].[CH4:8]
(2),RDTmap,[CH3:2][...][CH:4]>>[CH3:1][...][CH3:9].[CH4:2]
(1),CHYmap,[CH3:6][...][CH:7]>>[CH3:7][...][CH3:1].[CH4:6]
#,reaction_2
(1),RXNmap,[CH3:1][...][CH:8]>>[CH3:1][...][CH3:4].[CH3:5][...][CH4:8]
(2),RDTmap,[CH3:2][...][CH:4]>>[CH3:1][...][CH3:9].[CH4:2][...][CH4:8]
(3),CHYmap,[CH3:6][...][CH:7]>>[CH3:7][...][CH3:1].[CH3:8][...][CH3:4]
...
```

where the numer between parenthesis in the lines with the maps depicts the equivalence class of such map, meaning, for example, that for reaction_1, the maps RXNmap and CHYmap are equivalent to each other but not with RDTmap, while for reaction_2 the three given maps are all non-equivalent. Note that maps of each reaction in file (_aam_out_alleq) do not need an specification like this because they all are in the same equivalence class.<br/>

- Further input options:

If you only want to run the ITS method described in <a href="https://match.pmf.kg.ac.rs/issues/m90n1/m90n1_75-102.html">[1]</a>
```
python EEquAAM.py [myFile_aam.smiles]
```

If you want to run the 3 methods ITS, AUX and ISO described in <a href="https://match.pmf.kg.ac.rs/issues/m90n1/m90n1_75-102.html">[1]</a><br/>
(warning: ISO is more time consuming due to its mathematical properties)<br/>

```
python EEquAAM.py --sanity-check [myFile_aam.smiles]
```

- Mind the *.smiles extension when runnig these programs.<br/>


> you can find a more detailed explanation inside <a href="./EEquAAM.py">EEquAAM.py</a>
