- How to run the toy example (after activating eequaam conda environment)<br/>

```
python  NumberingTool.py  ToyExample.smiles
```

- This last command will produce the ToyExample_id.smiles file.<br/>


- General Input Format: plain text file with reaction SMILES as follows<br/>

```
CCCCCCCCCC>>CCCCCCCCC.C
CCCCCCCCCC>>CCCC.CCCCCC
...
```

- General Output Format: plain text file with reaction SMILES and identifiers as follows<br/>

```
#,reaction_1
CCCCCCCCCC>>CCCCCCCCC.C
#,reaction_2
CCCCCCCCCC>>CCCC.CCCCCC
...
```

- Mind the *.smiles extension when runnig these programs.<br/>


> you can find a more detailed explanation inside <a href="./NumberingTool.py">NumberingTool.py</a>