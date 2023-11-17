#!/usr/bin/env python3

from rxnmapper import RXNMapper
from rdkit import Chem
import sys
import itertools
from transformers import logging
from chython import smiles
import os

logging.set_verbosity_error()

reactions = sys.argv[1]
outputmapped = sys.argv[2] 


def isa_group_separator(line):
    return line=='\n'

with open(reactions, mode='r') as rf:
 with open(outputmapped, 'w') as omf:
   for key,group in itertools.groupby(rf,isa_group_separator):
      #print(key, list(group))
      if key:
         continue
      id_line, ec_line, name_reaction, smiles_reaction = map(str.strip, list(group))
      try:
          # RXN MAPPER
          rxn_mapper = RXNMapper()
          results_rxn = rxn_mapper.get_attention_guided_atom_maps([smiles_reaction], canonicalize_rxns=False)   
          
          # Graphormer
          print(smiles_reaction)
          r = smiles(smiles_reaction)
          r.reset_mapping()
          result_graphor = format(r, 'm')    
          
          
          # RDT mapper (in silet mode)  
          os.system("nohup java -jar RDT_2.4.1.jar -Q SMI -q " + "\"" + smiles_reaction + "\"" + " -c -j AAM -f TEXT >/dev/null 2>&1")
          inputFile = open("ECBLAST_smiles_AAM.txt", "r")
          results_rdt = (inputFile.read().splitlines())[3]
          inputFile.close()
          os.system("rm *.txt")
          os.system("rm *.rxn")
         
          # Wirte results
          smiles_sides = smiles_reaction.split('>>')
          left_string = ".".join( [ Chem.MolToSmiles(Chem.MolFromSmiles(s), allHsExplicit=True, canonical=False) for s in smiles_sides[0].split('.')  ] )
          right_string = ".".join( [ Chem.MolToSmiles(Chem.MolFromSmiles(s), allHsExplicit=True, canonical=False) for s in smiles_sides[1].split('.')  ] )    

          print(name_reaction, file=omf)
          print(left_string+">>"+right_string, file=omf)
          print(results_rxn[0]['mapped_rxn'], file=omf)
          print(result_graphor, file=omf)
          print(results_rdt, file=omf)
                
      except RuntimeError:
         print("skipped reaction of length", len(smiles_reaction), file=sys.stderr)
         print(id_line, file=sys.stderr)
