# EEquAAM - Welcome to the Repo


<p align="center">
<img src="./ReadmeLogo/EEquAAM_logo.png" width="400"/>
</p>


<p align="center">
Cheminformatics Software for the Automatic Evaluation of the Equivalence of Atom-to-Atom Maps
</p>


## Instructions


### Create Anaconda Environment
1) First create the eequaam environment while adding some of the required dependencies:
```
conda create -n eequaam python=3.9.13 networkx=2.8.4 matplotlib=3.5.2 numpy=1.21.5 beautifulsoup4=4.11.1
```
```
conda activate eequaam
```
```
pip install pysmiles==1.0.2 rdkit==2022.9.3 chytorch-rxnmap==1.3 rxnmapper==0.2.4 && pip install --no-warn-conflicts transformers==4.24.0
```


### Run MappingTool
```
python  MappingTool.py  [myFile.smiles]
```


### Run EEquAAM
```
python EEquAAM.py [myFile_aam.smiles]
```


### Deactivate or Remove Anaconda Environment
```
conda deactivate
```
```
conda env remove -n eequaam
```
