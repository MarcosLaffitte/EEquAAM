# EEquAAM - Welcome to the Repo


<p align="center">
<img src="./ReadmeLogo/EEquAAM_logo.png" width="460"/>
</p>


<div align="center">
Cheminformatics software for the automatic<br/>
<strong>E</strong>valuation of the <strong>Eq</strong>uivalence of <strong>A</strong>tom-to-<strong>A</strong>tom <strong>M</strong>aps
</div>


## Cite as


**[1]**   M. E. González Laffitte, N. Beier, N. Domschke, P. F. Stadler, Comparison of Atom Maps. *MATCH Commun. Math. Comput. Chem.* **90** (2023) 75–102.
> **Link:** https://match.pmf.kg.ac.rs/issues/m90n1/m90n1_75-102.html


## Instructions


### Create Anaconda Environment
###### 1) First create the eequaam environment while adding some of the required dependencies:
```
conda create -n eequaam python=3.9.13 networkx=2.8.4 matplotlib=3.5.2 numpy=1.21.5 beautifulsoup4=4.11.1
```
###### 2) Activate the eequaam conda environment:
```
conda activate eequaam
```
###### 3) Install pip dependencies inside the eequaam conda environment:
```
pip install pysmiles==1.0.2 rdkit==2022.9.3 chytorch-rxnmap==1.3 rxnmapper==0.2.4 && pip install --no-warn-conflicts transformers==4.24.0
```


### Run NumberingTool
```
python  NumberingTool.py  [myFile.smiles]
```


### Run MappingTool
```
python  MappingTool.py  [myFile_id.smiles]
```


### Run EEquAAM
###### There are two options when running EEquAAM. The first evaluates the equivalence of atom maps only through the ITS method. This is the fastest method and therefore the default option:
```
python  EEquAAM.py  [myFile_aam.smiles]
```
###### The second option evaluates 
```
python  EEquAAM.py  --sanity-check  [myFile_aam.smiles]
```


### Deactivate or Remove Anaconda Environment
###### You can deactivate the eequaam environment after using it:
```
conda deactivate
```
###### or remove it if not needed anymore:
```
conda env remove -n eequaam
```