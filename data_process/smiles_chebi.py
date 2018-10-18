# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 21:35:42 2018

@author: Administrator
"""

with open('E:/IGEM/data/ChEBI_complete.sdf') as f:
    f = f.readlines()
smiles_list = {}

for i,line in enumerate(f):
    if line.startswith('  Marvin  '):
        line = line.strip().split()
        marvin = line[1]
        smiles_list[marvin] = {}
    elif line.startswith('> <SMILES>'):
        smiles_list[marvin]['smiles'] = f[i+1].strip()
    elif line.startswith('> <KEGG COMPOUND'):
        smiles_list[marvin]['kegg_id'] = f[i+1].strip()
    if i%100 == 0:
        print(i)


smiles_list_1 = [[smiles_list[x]['kegg_id'], smiles_list[x]['smiles']] for x in smiles_list.keys() if 'smiles' in smiles_list[x] and 'kegg_id' in smiles_list[x]]
kegg_dict = {x[0]:x[1] for x in smiles_list_1}


with open('E:/IGEM/kegg_smiles.txt') as f:
    f = f.readlines()

result = ''

n=0
m=0
for line in f:
    line = line.strip().split()
    if not line[0] in kegg_dict:
        result += line[0] + '\t' + line[1] + '\n'
        n+=1
    else:
        result += line[0] + '\t' + kegg_dict[line[0]] + '\n'
        m+=1

l=0
for key in kegg_dict:
    #if not 'C' in key:
    #    print(key)
    s = key + '\t' +kegg_dict[key] + '\n'
    if not s in result:
        result += s
        l+=1

with open('E:/IGEM/kegg_smiles2.txt','w') as f:
    f.write(result)
        