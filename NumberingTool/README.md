- How to run the toy example (after activating eequaam conda environment)<br/>

       	   python  NumberingTool.py  ToyExample.smiles<br/>#


- This last command will produce the ToyExample_id.smiles file.


- Input format: plain text file with reaction SMILES as in the following example<br/>

CCCCCCCCCC>>CCCCCCCCC.C<br/>
CCCCCCCCCC>>CCCC.CCCCCC<br/>
...<br/>

- Output format: plain text file with reaction SMILES and identifiers as in the following example<br/>

#,reaction_1<br/>
CCCCCCCCCC>>CCCCCCCCC.C<br/>
#,reaction_2<br/>
CCCCCCCCCC>>CCCC.CCCCCC<br/>
...<br/>