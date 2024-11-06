# -*- coding: utf-8 -*-
"""
Input folder is the alphafold output folder
    change the directory "dirparent" and name of "folder"
User can choose to use the protein position as the label in the output file (NAMECHOICENUMBER = True)
    or the name of the protein in the file (NAMECHOICENUMBER = False)
    The code assumes that the protein names are separated by a "_" and that there are no extra "_" in the names
The script needs to find the .pdb file based on the .json file name
    If you get the error "name 'i_pdb' is not defined", adjust lines 119 and 120 to find the file name.

Output is that for every protein combination prediction, each pairwise comparison PAE plot is checked if it passes threshold
    each pair is a new line in the output excel file
    for ease of reading, each line is labelled as "self" if it is the protein against itself (AA) or "intermolecular" if not (AB)
    if the protein name is used, a third option is "homodimer" if the proteins have the same name but are in different positions (e.g. the first and second in AAB)
"""

# install/import dependencies and define functions
import numpy as np
import json
import os
from collections import defaultdict
import pandas as pd

# User choice on how to label the output proteins:
# "True" the position of the protein in the prediction (1, 2, 3, etc) - will always work
# "False" if the file is named as "protein1_protein2_protein3", this code will parse the name based on "_"
NAMECHOICENUMBER = True

# file paths
dirparent = r'Z:\Group_Data\Timothy_Soh\Documents\Collaborations\FLI_Axel-Karger\2024-11-05\test'
folder = 'results'
path = os.path.join(dirparent, folder)
filepath = os.path.join(dirparent, f'Scores_pairwise_{folder}.xlsx')

def parse_atm_record(line):
    '''Get the atm record
    '''
    record = defaultdict()
    record['name'] = line[0:6].strip()
    record['atm_no'] = int(line[6:11])
    record['atm_name'] = line[12:16].strip()
    record['atm_alt'] = line[17]
    record['res_name'] = line[17:20].strip()
    record['chain'] = line[21]
    record['res_no'] = int(line[22:26])
    record['insert'] = line[26].strip()
    record['resid'] = line[22:29]
    record['x'] = float(line[30:38])
    record['y'] = float(line[38:46])
    record['z'] = float(line[46:54])
    record['occ'] = float(line[54:60])
    record['B'] = float(line[60:66])

    return record

def read_pdb(pdbfile):
    '''Read a pdb file predicted with AF and rewritten to contain all chains
    '''

    chain_coords, chain_plddt = {}, {}
    with open(pdbfile, 'r') as file:
        for line in file:
            if not line.startswith('ATOM'):
                continue
            record = parse_atm_record(line)
            #Get CB - CA for GLY
            if record['atm_name']=='CB' or (record['atm_name']=='CA' and record['res_name']=='GLY'):
                if record['chain'] in [*chain_coords.keys()]:
                    chain_coords[record['chain']].append([record['x'],record['y'],record['z']])
                    chain_plddt[record['chain']].append(record['B'])
                else:
                    chain_coords[record['chain']] = [[record['x'],record['y'],record['z']]]
                    chain_plddt[record['chain']] = [record['B']]


    #Convert to arrays
    for chain in chain_coords:
        chain_coords[chain] = np.array(chain_coords[chain])
        chain_plddt[chain] = np.array(chain_plddt[chain])

    return chain_coords, chain_plddt

# list to store files
respdb = []
resjson = []

# make list of pdb files in dir
for file in os.listdir(path):
    # check only text files
    if file.endswith(("000.pdb")):
        respdb.append(file)
respdbmissing = respdb.copy()

# make list of json files in dir
for file in os.listdir(path):
    # check only text files
    if file.endswith(("000.json")):
        resjson.append(file)

# extract and calculate PAEs
scores = [[0,0,0,0,0,0,0,0,0]]
scores = np.delete(scores,0,0)
print("files found:  ", len(resjson))
for i in range(0,len(resjson)):
  if i%10==0:
      print("currently at count: "+ str(i) + ", file: " + resjson[i])
    
  # write file name into scores list
  i_filename=resjson[i][0:-14]

  # Reading metadatafile for a single prediction model
  f = open(path+'/'+resjson[i])
  metadata = json.load(f)
  pae=np.array(metadata['pae'])
  plddt=np.array(metadata['plddt'])

  # find pdb
  delim1 = [k for k,m in enumerate(resjson[i]) if m=='_']
  proteincombname = resjson[i][0:delim1[len(delim1)-10]+1]+"u"
  proteincombmodel = resjson[i][delim1[len(delim1)-4]:-5]
  proteincombfiles = []
  
  for j in range(0,len(respdb)):
      if proteincombname in respdb[j]:
          proteincombfiles.append(respdb[j])
  if not len(proteincombfiles)==5:
      print(f'Problem: {resjson[i]} \n{proteincombname} \n{proteincombfiles}')
          
  for j in range(0,len(proteincombfiles)):
      if proteincombmodel in proteincombfiles[j]:
          i_pdb=proteincombfiles[j]

  respdbmissing.remove(i_pdb)

  # read pdb
  pdbpath=(path+'/'+i_pdb)
  chain_coords, chain_plddt = read_pdb(pdbpath)

  ## Read chain lengths out of pdb file because I cannot find the individual chain lengths in the PAE json files 
  chainlength=np.empty(len(chain_coords), dtype="int")
  numchains=0;
  for j in chain_coords:
      chainlength[numchains]=chain_coords[j].shape[0]
      numchains=numchains+1

  ## make arrays of intermolecular PAE values
  paechainall=np.empty([numchains,numchains], dtype="object")
  xstart=0
  for j in range(0,numchains):
      ystart=0
      for k in range(0,numchains):
          paetemp=pae[xstart:xstart+chainlength[j], ystart:ystart+chainlength[k]]
          if k>=j:
              paechainall[j,k]=paetemp.flatten()
          else:
              paechainall[k,j]=np.append(paechainall[k,j],paetemp.flatten())
          ystart=ystart+chainlength[k]
      xstart=xstart+chainlength[j]

  #grab pTM and piTM
  ptm=np.array(metadata['ptm'])
  iptm=np.array(metadata['iptm'])
  
  ## calculate from PAE scores max value and sum over background
  paescores=np.empty([1, 8+1], dtype="object")
  for j in range(0,numchains):
      for k in range(0,numchains):
        if k>=j:
          PAEmin = np.amin(paechainall[j,k])
          PAEmedian = np.median(paechainall[j,k])
          PAEmean = np.mean(paechainall[j,k])
          PAEstd = np.std(paechainall[j,k])
          
          paescores[0,0] = f'{i_filename} Chain {j} with Chain {k}' #ptm
          paescores[0,1] = np.round(np.max(ptm),2) #ptm
          paescores[0,2] = np.round(np.max(iptm),2) #iptm
          paescores[0,3] = 30-PAEmin #strongest signal
          paescores[0,4] = PAEstd/PAEmean #how broad is the distribution to identify muliple populations/size        
          if paescores[0,3]>24.84 and paescores[0,4]>0.2131:
              paescores[0,5] = "Yes"
          else:
              paescores[0,5] = "No"

          if NAMECHOICENUMBER:
              paescores[0,6] = j
              paescores[0,7] = k
          else:
              if j == 0:
                  startj = 0
              else:
                  startj = delim1[j-1]+1
              endj = delim1[j]
           
              if k == 0:
                  startk = 0
              else:
                  startk = delim1[k-1]+1
              endk = delim1[k]
           
              paescores[0,6] = resjson[i][startj:endj]
              paescores[0,7] = resjson[i][startk:endk]

          if j==k:
              paescores[0,8] = 'self'
          elif paescores[0,6] == paescores[0,7]:
              paescores[0,8] = 'homodimer'
          else:
              paescores[0,8] = 'intermolecular'
          scores = np.append(scores,paescores,0)

## save scores
scores = pd.DataFrame(scores)
with pd.ExcelWriter(filepath) as writer:
    scores.to_excel(writer, header=["file", "pTM", "ipTM", "30-PAEmin", "PAEstd/PAEmean", "interaction?", "Protein 1", "Protein 2", "self or intermolecular?"])


print("\nfiles found:  ", len(resjson))
print('done calculating and saved in: \n'+filepath)

print('\n' + str(len(respdb)-len(resjson)) + " pdb file(s) are missing a .json file")
if (len(respdb)-len(resjson)) != 0:
    print(respdbmissing)
