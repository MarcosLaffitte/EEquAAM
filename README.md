# EEquAAM - Welcome to the Repo


<p align="center">
<img src="./ReadmeLogo/EEquAAM_logo.png" width="400"/>
</p>


<div align="center">
Cheminformatics software for the automatic<br/>
<strong>E</strong>valuation of the <strong>Equ</strong>ivalence of <strong>A</strong>tom-to-<strong>A</strong>tom <strong>M</strong>aps
</div>


## Institutions

> Center for Scalable Data Analytics and Artificial Intelligence, Leipzig / Dresden, Germany. See <a href="https://scads.ai/">ScaDS.AI</a>.<br/>

> Bioinformatics Group, Department of Computer Science, Leipzig University, Germany. See <a href="https://www.bioinf.uni-leipzig.de/home.html">Bioinf</a>.<br/>

> Interdisciplinary Center for Bioinformatics, Leipzig University, Germany. See <a href="https://www.izbi.uni-leipzig.de/">IZBI</a>.<br/>



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
Here we provide a self-contained pipeline for the analysis and comparison of atom-to-atom maps explaining the mechanisms behind chemical reactions. You will find 3 main python tools: <a href="./EEquAAM">EEquAAM</a> will compare atom maps in an automatic manner. If you still don't have annotated reaction SMILES describing atom maps, you can use our <a href="./MappingTool">MappingTool</a>, which implements the 3 representative atom mapping tools RXN, RDT and Graphormer; find more information about them below. Finally, to make use of EEquAAM and of our MappingTool you will require a list of reaction SMILES displaying names / identifiers for each reaction. If you have a list of such SMILES but don't want to name them one-by-one, you may use our auxiliary <a href="./NumberingTool">NumberingTool</a>.
</div>

<br/>

**Note:** the atom maps to be compared need to be complete, i.e., *bijective*. An incomplete atom map is that in which there are atoms on one side of the reaction not being mapped to atoms on the other side, in other words, the atom map does not represent a bijection. This in turn requires for the analyzed reactions to be stoichiometrically balanced. For more information see the following reference.<br/>


## Cite as

This repository was developed as part of the contribution:

**[1]**   M. E. González Laffitte, N. Beier, N. Domschke, P. F. Stadler, Comparison of Atom Maps. *MATCH Commun. Math. Comput. Chem.* **90** (2023) 75–102.
> **Link:** https://match.pmf.kg.ac.rs/issues/m90n1/m90n1_75-102.html

<div align="justify">
There you can find detailed information on the algorithms implemented here. This work was developed for research purposes. Please cite as above if you find this work or these programs useful for your own research.
</div>


## Instructions

###### In order to run these programs you will require some python packages, below we show how to install them directly into an anaconda environment. In particular, if you are runnig MappingTool you will need to install Java or OpenJDK in your computer (required to run the RDT mapper). After this you only need a list of unannotated reaction SMILES inside a plain-text file, over which you can apply the pipeline NumberingTool > MappingTool > EEquAAM, meaning that the output of one program is the input for the following. Below you can find the python commands to run each program as well.


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
###### Remember to always activate the eequaam conda environment before using the programs in this repository.


### Run NumberingTool
###### If your initial list of unannotated reaction SMILES doesn't have identifiers for each reaction, you can make use of this script to add them to your list. In the following command [myFile.smiles] is a plain-text file whose lines are each a reaction SMILES. Please mind the *.smiles extension:
```
python  NumberingTool.py  [myFile.smiles]
```
###### The output will be a plain-text file [myFile_id.smiles] whose lines alternate between the reaction SMILES and their identifiers. See <a href="./NumberingTool">NumberingTool</a> for an example of the input and output formats.<br/>



### Run MappingTool
###### If you already have a list of unannotated reaction SMILES with identifiers and you want to produce atom maps for these reactions, you can make use of this script. This wraps the three atom mapping tools RXN, RDT and Graphormer.
```
python  MappingTool.py  [myFile_id.smiles]
```
###### The output will include a plain-text file [myFile_aam.smiles] containing the annotated SMILES of those suitable and balanced reactions that were completely mapped by the three mappers. The extension of these files is again *.smiles. Nevertheless they have different formats. See <a href="./MappingTool">MappingTool</a> for an example of the input and output formats.<br/>



### Run EEquAAM
###### There are two options when running EEquAAM. The first evaluates the equivalence of atom maps only through the ITS method. This is the fastest method and therefore the default option:
```
python  EEquAAM.py  [myFile_aam.smiles]
```
###### The second option evaluates the equivalence of atom maps with the three methods ITS, AUX and ISO described in <a href="https://match.pmf.kg.ac.rs/issues/m90n1/m90n1_75-102.html">[1]</a>. These are mathematically equivalent, but the ISO method may be more time consuming depending on the number of symmetries of the involved molecules. To run and compare these methods use the following command:
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

- RXNMapper [RXN]<br/>
  http://rxnmapper.ai/<br/>
  https://pypi.org/project/rxnmapper/<br/>

- Reaction Decoder Tool [RDT]<br/>
  https://academic.oup.com/bioinformatics/article/32/13/2065/1743096<br/>
  https://github.com/asad/ReactionDecoder<br/>
  https://github.com/asad/ReactionDecoder/releases<br/>
  
- GraphormerMapper [Graphormer]<br/>
  https://pubs.acs.org/doi/10.1021/acs.jcim.2c00344<br/>
  https://github.com/chython/chytorch-rxnmap<br/>
  https://pypi.org/project/chytorch-rxnmap/<br/>

- General specifications on both anotated and unannotated SMILES<br/>
  http://opensmiles.org/opensmiles.html<br/>



## LICENSE

The programs in this repository are part of the work published in<br/>
https://match.pmf.kg.ac.rs/issues/m90n1/m90n1_75-102.html<br/>
and are released under<br/>
<strong>MIT License Copyright (c) 2023 Marcos E. González Laffitte</strong><br/>
See <a href="./LICENSE">LICENSE</a> file for full license details.