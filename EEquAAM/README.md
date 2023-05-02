- How to run the toy example (after activating eequaam conda environment)<br/>

```
python  EEquAAM.py  ToyExample_id_aam.smiles
```

- This last command will produce (in particular) the ToyExample_id_aam_summary.txt and the ToyExample_id_aam_times.pdf files.<br/>


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

b) [*_aam_out_alleq.txt] and [*_aam_out_noneq.txt]: reactions with all equivalent maps and with some non-equivalent maps.<br/>

c) [*_aam_times.pdf]: PDF boxplots of the time taken by each method to compare the provided atom maps.<br/>


- Input options:

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
