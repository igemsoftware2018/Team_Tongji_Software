# -*- coding: utf-8 -*-
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem.Fingerprints import FingerprintMols
import csv
from rdkit import Chem
from rdkit.Chem import Draw


compound_dict = dict()
with open("data//compound_name_pair.csv") as c_file:
    c_reader = csv.reader(c_file)
    for c in c_reader:
        compound_dict.update({c[0]: c[1:]})
    compound_dict.update({'':''})

def save_img(smiles, img_path, x=300, y=300):
    mol = Chem.MolFromSmiles(smiles)
    Draw.MolToFile(mol, img_path, size=(x,y))

def get_similar_compound(condon):
    com=condon['smiles']

    save_img(com,'static//compound_img//smiles_img.png',300,300)

    output_num=condon['MaxLength']
    smiles_file_path = 'data//kegg_smiles2.txt'
    
    with open(smiles_file_path) as file:
        f = file.readlines()
        
    smiles_list = [x.split()[1] for x in f]
    output_num = min(output_num, len(smiles_list))
    top_idx = [0] * output_num
    top_score = [0] * output_num
    
    
    mol1 = Chem.MolFromSmiles(com)
    if mol1 is None:
        print('input smiles not exist')
        return []
    mol1 = AllChem.AddHs(mol1)
    fps1 = AllChem.GetMorganFingerprint(mol1, 2)
    
    for i,item in enumerate(smiles_list):
        mol2 = Chem.MolFromSmiles(item)
        if mol2 is None:
            continue
        mol2 = AllChem.AddHs(mol2)
        fps2 = AllChem.GetMorganFingerprint(mol2, 2)

        score = DataStructs.DiceSimilarity(fps1, fps2)
        score = round(score,2)
        
        if score > min(top_score):
            min_idx = top_score.index(min(top_score))
            top_idx[min_idx] = i
            top_score[min_idx] = score
        
    top_keggid = [f[i].split()[0] for i in top_idx]
    top_smiles = [f[i].split()[1] for i in top_idx]
    result=sorted(zip(top_keggid, top_smiles, top_score), key=lambda x:x[2], reverse=True)
    for i in range(len(result)):
        result[i]=list(result[i])
        result[i].insert(1,compound_dict[result[i][0]][0])
    return result
if __name__=='__main__':
    print(get_similar_compound('CCCO',10))