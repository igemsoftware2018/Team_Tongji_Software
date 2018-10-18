with open('E:/IGEM/atom/data/compounds.dat','rb') as file:
    f = file.readlines()

compounds = []
for line in f:
    if compounds == [] and not line.startswith(b'UNIQUE-ID - '):
        continue
    if line.startswith(b'UNIQUE-ID - '):
        compounds.append({'id' : line.strip()[12:]})
    elif line.startswith(b'COMMON-NAME - '):
        compounds[-1]['name'] = line.strip()[14:]
    elif line.startswith(b'SMILES - '):
        compounds[-1]['smiles'] = line.strip()[9:]
    elif line.startswith(b'CHEMICAL-FORMULA - '):
        if 'formula' in compounds[-1].keys():
            compounds[-1]['formula'].append(line.strip()[18:])
        else:
            compounds[-1]['formula'] = [line.strip()[18:]]



smiles_dict = {c['id']: c['smiles'] for c in compounds if 'smiles' in c.keys()}
smiles_dict.update({c['name']: c['smiles'] for c in compounds if 'name' in c.keys() and 'smiles' in c.keys()})

smiles_list = [x.decode() for x in list(smiles_dict.values())]
atom_list = []
for s in smiles_list:
    s = s.replace('DNA','').replace('RNA','')
    for i,w in enumerate(s):
        if w.isupper():
            if i == len(s)-1:
                atom_list.append(w)
            elif s[i+1].islower():
                atom_list.append(s[i:i+2])
            else:
                atom_list.append(w)
atom_list = list(set(atom_list))
atom_list.remove('H')



import os
mol_dict = {}
mol_path = 'E:/IGEM/data/atom/data/MetaCyc-MOLfiles/'
path_list = os.listdir(mol_path)
for path in path_list:
    with open(mol_path + path,'rb') as file:
        f = file.readlines()
    mol = []
    for line in f:
        line = line.split()
        if len(line) >= 4 and line[3].strip()[:1].isupper() and (line[3].strip().isalpha() or line[3].strip() == b'R#'):
            atom = line[3].strip().decode()
            if len(atom)<3 and atom != 'H':
                mol.append(atom)
    mol_dict[path[:-4]] = mol

'''
def comma(l):
    for i, item in enumerate(l):
        list[i] = b'"' + list[i] + b'"'
    return list


result_csv = b''
for item in result:
    result_csv += b','.join(comma(item)) + b'\n'
with open('C:/Users/Lenovo/Desktop/igem/22.0/atom/compounds.csv','wb') as file:
    file.write(result_csv)

'''

with open('E:/IGEM/data/atom/data/reactions.dat', 'rb') as file:
    f = file.readlines()

reaction=[]
for i in f:
    if i.startswith(b'UNIQUE-ID - '):
        reaction.append([i.strip()[12:]])
    elif i.startswith(b'ATOM-MAPPINGS - '):
        reaction[-1].append(i.strip()[16:])
    #elif i.startswith(b'EC-NUMBER - EC- '):
    #    reaction[-1].append(i.strip()[16:])
    elif i.startswith(b'LEFT - '):
        reaction[-1].append(b'left-' + i.strip()[7:])
    elif i.startswith(b'RIGHT - '):
        reaction[-1].append(b'right-' + i.strip()[8:])

reaction = [x for x in reaction if b'ENCODING' in x[1]] #NO ENCODING



# metabolic-reactions
import re

with open('E:/IGEM/atom/data/metabolic-reactions.xml', 'rb') as file:
    f = file.readlines()
biocyc_kegg = []
for line in f:
    if b'BIOCYC:' in line:
        biocyc_kegg.append([re.findall(b'BIOCYC: (.*?)</p>',line)[0]])
    #elif b'PUBCHEM:' in line:
    #    biocyc_kegg[-1].append(re.findall(b'PUBCHEM: (\d*?)</p>',line)[0])
    elif b'KEGG:' in line:
        biocyc_kegg[-1].append(re.findall(b'KEGG: ([CDGcR]?\d*?)</p>',line)[0])
    # 13096 / 36954
biocyc_kegg = [x for x in biocyc_kegg if len(x)==2]
biocyc_kegg = list(set([tuple(x) for x in biocyc_kegg]))
biocyc_kegg = [(x[0].decode(), x[1].decode()) for x in biocyc_kegg]
# 9627
# R:4546 / CDG:5081
biocyc_kegg_r = [x for x in biocyc_kegg if 'R' in x[1]]
biocyc_kegg_c = [x for x in biocyc_kegg if not 'R' in x[1]]

valid_reaction = [x for x in reaction if x[0].decode() in [y[0] for y in biocyc_kegg_r]]

'''
with open('E:/IGEM/atom/compound_raw.txt', 'rb') as file:
    f = file.readlines()
kegg_pubchem = []
for i,line in enumerate(f):
    if line.startswith(b'ENTRY  '):
        kegg_pubchem.append([re.findall(b'C\d{5}',line)[0]])
    elif b'PubChem:' in line:
        kegg_pubchem[-1].append(re.findall(b'PubChem: (\d*)',line)[0])
    # 18261 / 18294
kegg_pubchem = [x for x in kegg_pubchem if len(x) == 2]
k_id = [int(x[1]) for x in kegg_pubchem]


biocyc_kegg = []
pubid_in_kegg = [x[1] for x in kegg_pubchem]
pubid_in_biocyc = [x[1] for x in biocyc_pubchem]
for i,pubid in enumerate(pubid_in_biocyc):
    if pubid in pubid_in_kegg:
        i = i 
        j = pubid_in_kegg.index(pubid)
        biocyc_kegg.append([biocyc_pubchem[i][0], kegg_pubchem[j][0], pubid])
'''

def smiles_purify(smiles):
    purified = []
    for i,w in enumerate(smiles):
        if w.isupper():
            if i==len(smiles)-1:
                if w in atom_list:
                    purified.append(w)
            else:
                if not smiles[i+1].islower():
                    if w in atom_list:
                        purified.append(w)    
                elif smiles[i:i+2] in atom_list:
                    purified.append(smiles[i:i+2])   
    purified.sort()
    return purified


def same_idx(a1,a2,b1,b2):
    min_idx = max(a1,b1)
    max_idx = min(a2,b2)
    return [min_idx,max_idx]


import re

def reaction_process(encoding_text):
    encoding = re.findall(b'NO-HYDROGEN-ENCODING \(.*?\)', encoding_text)
    encoding_text = encoding_text[(len(encoding[0])+2):]
    
    msg_list = re.findall(b'(\([^\(]*?\s\d*?\s\d*?\))', encoding_text)
    msg_list = [encoding[0][21:]] + msg_list
    
    msg_list = [x[1:-1] for x in msg_list]
    return msg_list

a_c = []
error_list = []
for i,item in enumerate(reaction):
    com = [[],[]]
    com_smiles = [[],[]]
    com_idx = [[],[]]
    
    rea_process = reaction_process(item[1])
    left_idx = [int(x) for x in rea_process[0].decode().split()]
    right_idx = list(range(max(left_idx)+1))
    
    error = False
    for c in rea_process[1:]:
        c = c.split()
        #example: reaction[2] |CoM| / reaction[764] "CPD-415"
        if c[0][:1] + c[0][-1:] == b'||' or c[0][:1] + c[0][-1:] == b'""':
            c[0] = c[0][1:-1]
        c_name = c[0] if len(c) == 3 else b' '.join(c[:2])
        
        if c[-2] == b'0':
            if com[0] == []:
                j = 0
            else:
                j = 1
        else:
            if com[1] == []:
                j = 0
            else:
                j = 1
                
        com[j].append(c_name)
        com_idx[j].append(list(range(int(c[-2]),int(c[-1])+1)))
        try:
            com_smiles[j].append(mol_dict[c[0].decode()])
        except KeyError as e:
            error_list.append([i,e.args[0]])
            error = True
            break
        '''
        # NAD / NADH / 1124 DNA-containing-a-Apyrimidinic-Sites / 6926 CPD0-1028
        if b'left-'+c[0] in item or (c[0]==b'NAD' and b'left-Oxidized-ferredoxins' in item)\
                                 or (c[0]==b'NADH' and b'left-Reduced-ferredoxins' in item)\
                                 or (c[0]==b'CPD0-1028' and b'left-GERANYLGERANYL-PP' in item):
            com[0].append(c_name)
            com_idx[0].append(list(range(int(c[-2]),int(c[-1])+1)))
            try:
                com_smiles[0].append(smiles_purify(smiles_dict[c[0]].decode()))
            except KeyError as e:
                error_list.append([i,e.args[0]])
                error = True
                break
        elif b'right-'+c[0] in item or (c[0]==b'NAD' and b'right-Oxidized-ferredoxins' in item)\
                                    or (c[0]==b'NADH' and b'right-Reduced-ferredoxins' in item)\
                                    or (c[0]==b'DNA-containing-a-Apyrimidinic-Sites' and b'right-DNA-containing-aPurinic-Sites' in item)\
                                    or (c[0]==b'CPD0-1028' and b'right-GERANYLGERANYL-PP' in item):
            com[1].append(c_name)
            com_idx[1].append(list(range(int(c[-2]),int(c[-1])+1)))
            try:
                com_smiles[1].append(smiles_purify(smiles_dict[c[0]].decode()))
            except KeyError as e:
                error_list.append([i,e.args[0]])
                error = True
                break
        else:
            error_list.append([i,'NotInList: %s' % c[0].decode()])
            error = True
            break
        '''
    if error:
        continue

    a_c.append([i,com,com_smiles,com_idx,left_idx,right_idx])


atom_idx_list = {reaction[x[0]][0].decode():x[3] for x in a_c}


unbalanced_list = []
for i,item in enumerate(a_c):
    error_message = [i]
    com_smiles = item[2]
    com_idx = item[3]
    smiles_len = [[],[]]
    idx_len = [[],[]]
    for j in range(2):
        smiles_len[j] = [len(x) for x in com_smiles[j]]
        idx_len[j] = [len(x) for x in com_idx[j]]
    if sum(smiles_len[0]) != sum(smiles_len[1]):
        error_message.append('smiles_unbalanced')
    if sum(idx_len[0]) != sum(idx_len[1]):
        error_message.append('idx_unbalanced')
    if smiles_len != idx_len:
        error_message.append('len_unbalanced')
        
    if error_message != [i]:
        unbalanced_list.append(error_message)
# num: 1


import copy
num=0
atom_transfer = {}
for item in a_c:
    biocyc_id = reaction[item[0]][0].decode()
    
    rea_atom = item[2][0]
    pro_atom = item[2][1]
    rea_idx = item[3][0]
    pro_idx = item[3][1]
    rea_list = item[4]
    pro_list = item[5]
    transfer_list = [copy.deepcopy(rea_idx), copy.deepcopy(pro_idx)]
    
    flag = True
    for i_r,com_r in enumerate(rea_atom):
        for j_r,atom_r in enumerate(com_r):
            idx_r = rea_idx[i_r][j_r]
            idx_p = pro_list[rea_list.index(idx_r)]
            for i_p in range(len(pro_idx)):
                if idx_p in pro_idx[i_p]:
                    j_p = pro_idx[i_p].index(idx_p)
                    break
            if atom_r != pro_atom[i_p][j_p]:
                print(biocyc_id,'Mismatched: %i-%i / %i-%i' % (i_r,j_r,i_p,j_p))
                flag = False
                break
            transfer_list[0][i_r][j_r] = i_p
            transfer_list[1][i_p][j_p] = i_r
    if flag:
        com_list = [[x.decode() for x in item[1][0]],[x.decode() for x in item[1][1]]]
        atom_transfer[biocyc_id] = [transfer_list, com_list]


corresponding_com = []
with open('E:/IGEM/biocyc-kegg-compound-2.txt') as f:
    f = f.readlines()
for line in f:
    corresponding_com.append(line.strip().split())

corresponding_rea = []
with open('E:/IGEM/biocyc-kegg-reaction-2.txt') as f:
    f = f.readlines()
for line in f:
    corresponding_rea.append(line.strip().split())



biocyc_compound = list(set([x[0] for x in corresponding_com]))
kegg_compound = list(set([x[1] for x in corresponding_com]))
biocyc_reaction = list(set([x[0] for x in corresponding_rea if x[0] in atom_transfer.keys()]))
kegg_reaction = list(set([x[1] for x in corresponding_rea if x[0] in atom_transfer.keys()]))

#corresponding_rea_dict = {b:k for b,k in corresponding_rea if b in atom_transfer.keys()}
#atom_conservation = {corresponding_rea_dict[b]:atom_transfer[b] for b in corresponding_rea_dict.keys()}


def get_corresponding(item, data_type): #0:compound b2k  1:compound k2b  2:reaction b2k  3:reaction k2b
    result_list = []
    if data_type < 2:
        if data_type == 0:
            for com in corresponding_com:
                if com[0] == item:
                    result_list.append(com[1])
        if data_type == 1:
            for com in corresponding_com:
                if com[1] == item:
                    result_list.append(com[0])
    else:
        if data_type == 2:
            for rea in corresponding_rea:
                if rea[0] == item:
                    result_list.append(rea[1])
        if data_type == 3:
            for rea in corresponding_rea:
                if rea[1] == item:
                    result_list.append(rea[0])
    return result_list
    
'''
import pickle
output = open('E:/IGEM/atom/atom.pkl','wb')
pickle.dump(atom_conservation, output)
output.close()
'''

'''
pkl_file = open('E:/IGEM/atom/atom.pkl','rb')
atom_conservation = pickle.load(pkl_file)
'''




import pickle
pkl_file = open('E:/IGEM/web/software-test4/data/reaction_dictionary.pk','rb')
reaction_pair = pickle.load(pkl_file)


result = []
for pair in reaction_pair.keys():
    k_rea_gruop = reaction_pair[pair]
    for k_rea_id in k_rea_gruop:
        if k_rea_id in kegg_reaction:
            com1, com2 = pair.split('_')
            if com1 in kegg_compound and com2 in kegg_compound:
                for b_rea_id in get_corresponding(k_rea_id, 3):
                    if b_rea_id in atom_transfer.keys():
                        transfer_list, com_list = atom_transfer[b_rea_id]
                        exist1 = True
                        exist2 = True
                        for b_com1_id in get_corresponding(com1, 1):
                            if b_com1_id in com_list[0]:
                                i_1, j_1 = 0, com_list[0].index(b_com1_id)
                                break
                            elif b_com1_id in com_list[1]:
                                i_1, j_1 = 1, com_list[1].index(b_com1_id)
                                break
                            exist1 = False
                        for b_com2_id in get_corresponding(com2, 1):
                            if b_com2_id in com_list[0]:
                                i_2, j_2 = 0, com_list[0].index(b_com2_id)
                                break
                            elif b_com2_id in com_list[1]:
                                i_2, j_2 = 1, com_list[1].index(b_com2_id)
                                break
                            exist2 = False
                        if exist1 and exist2:
                            atom_idx = atom_idx_list[b_rea_id]
                            atom_idx_data1 = atom_idx[i_1][j_1]
                            atom_idx_data2 = atom_idx[i_2][j_2]
                            transfer_data1 = [x if x in atom_idx_data2 else -1 for x in atom_idx_data1]
                            transfer_data2 = [x if x in atom_idx_data1 else -1 for x in atom_idx_data2]
                            
                            if len(transfer_data1)-transfer_data1.count(-1) == len(transfer_data2)-transfer_data2.count(-1):
                                result.append([pair,k_rea_id,transfer_data1,transfer_data2])
                                

'''
output = open('E:/IGEM/atom/atom_transfer.pkl','wb')
pickle.dump(result, output)
output.close()
'''

        
        






'''
result = ''
for item in bio_result:
    if item in atom_transfer:
        result += item
        result += '/'
        for side in atom_transfer[item]:
            for com in side:
                result += ','.join([str(x) for x in com]) + ' '
            result += '/'
        result += '\n'
with open('E:/IGEM/atom_conservation.txt','w') as f:
    f.write(result)
'''