- How to run the toy example (after activating eequaam conda environment)<br/>

```
python  MappingTool.py  ToyExample_id.smiles
```

- This last command will produce (in particular) the ToyExample_id_aam.smiles file when blanaced and suitable reactions are found.<br/>


- General Input Format: plain-text file with reaction SMILES and identifiers as follows (this is also the output of NumberingTool.py)<br/>

```
#,reaction_1
CCCCCCCCCC>>CCCCCCCCC.C
#,reaction_2
CCCCCCCCCC>>CCCC.CCCCCC
...
```

- General Output Format: plain-text file with the original identifiers and the mapped SMILES as follows<br/>

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

> you can find a more detailed explanation inside <a href="./MappingTool.py">MappingTool.py</a>