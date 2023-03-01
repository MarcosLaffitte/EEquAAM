# EEquAAM - Welcome to the Repo


<p align="center">
<img src="./ReadmeLogo/EEquAAM_logo.png" width="420"/>
</p>


<div align="center">
Cheminformatics software for the automatic<br/>
<strong>E</strong>valuation of the <strong>Eq</strong>uivalence of <strong>A</strong>tom-to-<strong>A</strong>tom <strong>M</strong>aps
</div>


## Developed by

- Marcos E. González Laffitte<br/>
  Bioinformatics Group and ScaDS.AI, Leipzig University<br/>
  marcoslaffitte@gmail.com<br/>
  marcos@bioinf.uni-leipzig.de<br/>

- Nora Beier<br/>
  Bioinformatics Group, Leipzig University<br/>
  nora@bioinf.uni-leipzig.de<br/>

- Nico Domschke<br/>
  Bioinformatics Group, Leipzig University<br/>
  dnico@bioinf.uni-leipzig.de<br/>

- Prof. Dr. Peter F. Stadler<br/>
  Bioinformatics Group and ScaDS.AI and Interdisciplinary Center for Bioinformatics, Leipzig University<br/>
  studla@bioinf.uni-leipzig.de<br/>
  

## Description

<div align="justify">
Here we provide a self-contained pipeline for the analysis and comparison of atom-to-atom maps explaining the mechanisms behind chemical reactions. You will find 3 main python tools: <a href="./EEquAAM">EEquAAM</a> will compare atom maps in an automatic manner. If you still don't have annotated reaction SMILES describing atom maps, you can use our <a href="./MappingTool">MappingTool</a>, which implements the 3 representative atom mapping tools RXN, RDT and Graphormer, see information about them below. Finally, to make use of EEquAAM and of our MappingTool you will require a list of reaction SMILES displaying names / identifiers for each reaction. If you have a list of such SMILES but don't want to name them one-by-one, you may use our auxiliary <a href="./NumberingTool">NumberingTool</a>.
</div>


## Cite as

This repository was developed as part of the contribution:

**[1]**   M. E. González Laffitte, N. Beier, N. Domschke, P. F. Stadler, Comparison of Atom Maps. *MATCH Commun. Math. Comput. Chem.* **90** (2023) 75–102.
> **Link:** https://match.pmf.kg.ac.rs/issues/m90n1/m90n1_75-102.html

<div align="justify">
There you can find detailed information on the algorithms implemented here. This work was developed for research purposes. Please cite as above if you find this information useful for your own research.
</div>


## Instructions

###### In order to run these programs you will require some python packages, which can be installed inside an anaconda environment. After this you only need a list of unannotated reaction SMILES inside a plain-text file, over which you can apply the pipeline NumberingTool > MappingTool > EEquAAM, meaning that the output of one program is the input for the following. Below you can find the python commands to run each program.


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
###### If your initial list of unannotated reaction SMILES doesn't have identifiers for each reaction, you can make use of this script to add them to your list. In the following command [myFile.smiles] is a plain-text file whose lines are each a reaction SMILES:
```
python  NumberingTool.py  [myFile.smiles]
```
###### The output will be a plain-text file [myFile_id.smiles] whose lines alternate between the reaction SMILES and their identifiers. See <a href="./NumberingTool">NumberingTool</a> for an example of the input and output formats.<br/>



### Run MappingTool
###### If you already have a list of unannotated reaction SMILES with identifiers and you want to produce atom maps for these reactions, you can make use of this script. This wraps the three atom mapping tools RXN, RDT and Graphormer.
```
python  MappingTool.py  [myFile_id.smiles]
```
###### The output will include a plain-text file [myFile_aam.smiles] containing the annotated SMILES of those suitable and balanced reactions that were completely mapped by the three mappers. See <a href="./MappingTool">MappingTool</a> for an example of the input and output formats.<br/>



### Run EEquAAM
###### There are two options when running EEquAAM. The first evaluates the equivalence of atom maps only through the ITS method. This is the fastest method and therefore the default option:
```
python  EEquAAM.py  [myFile_aam.smiles]
```
###### The second option evaluates the equivalence of atom maps with the three methods ITS, AUX and ISO described in [1]. These are mathematically equivalent, but the ISO method may be more time consuming depending on the number of symmetries of the involved molecules. To run and compare these methods use the following command:
```
python  EEquAAM.py  --sanity-check  [myFile_aam.smiles]
```
###### The output will include a summary indicating for which reactions the given maps were all equivalent, and for which the given maps were non-equivalent. See <a href="./EEquAAM">EEquAAM</a> for an example of the input and output formats.<br/>



### Deactivate or Remove Anaconda Environment
###### Finally you can deactivate the eequaam environment after using it:
```
conda deactivate
```
###### or remove it if not needed anymore:
```
conda env remove -n eequaam
```


## Additional Information


