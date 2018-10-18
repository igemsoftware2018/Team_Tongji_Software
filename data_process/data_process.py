# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 15:41:30 2018

@author: Administrator
"""

import re

corresponding_com = []
corresponding_rea = []

# directly corresponding

with open('E:/IGEM/data/atom/data/metabolic-reactions.xml', 'rb') as file:
    f = file.readlines()
biocyc = []
for line in f:
    if b'BIOCYC:' in line:
        name = re.findall(b'BIOCYC: (.*?)</p>',line)[0]
        if name[0] == b'|' and name[-1] == b'|':
            name = name[1:-1]
        biocyc.append([name])
    #elif b'PUBCHEM:' in line:
    #    biocyc_pubchem[-1].append(re.findall(b'PUBCHEM: (\d*?)</p>',line)[0])
    elif b'KEGG:' in line:
        biocyc[-1].append(re.findall(b'KEGG: ([CDGcR]?\d*?)</p>',line)[0])
    # 13096 / 36954
biocyc_kegg = [x for x in biocyc if len(x)==2]
biocyc_rest = [x[0] for x in biocyc if len(x) == 1]

com_list = [[item.decode() for item in x] for x in biocyc_kegg if not b'R' in x[1]]
for item in com_list:
    if not item in corresponding_com:
        corresponding_com.append(item)
# 5081
rea_list = [[item.decode() for item in x] for x in biocyc_kegg if b'R' in x[1]]
for item in rea_list:
    if not item in corresponding_rea:
        corresponding_rea.append(item)
# 4546
# no repeating
corresponding_com = [tuple(x) for x in corresponding_com]
corresponding_rea = [tuple(x) for x in corresponding_rea]


# compound 
#   biocyc

with open('E:/IGEM/atom/data/compounds.dat','rb') as file:
    f = file.readlines()
biocyc_com = []
for line in f:
    if biocyc_com == [] and not line.startswith(b'UNIQUE-ID - '):
        continue
    if line.startswith(b'UNIQUE-ID - '):
        biocyc_com.append({'id' : line.strip()[12:]})
    elif line.startswith(b'COMMON-NAME - '):
        biocyc_com[-1]['name'] = line.strip()[14:]
    elif line.startswith(b'SMILES - '):
        biocyc_com[-1]['smiles'] = line.strip()[9:]
    elif line.startswith(b'CHEMICAL-FORMULA - '):
        if 'formula' in biocyc_com[-1].keys():
            biocyc_com[-1]['formula'].append(line.strip()[18:])
        else:
            biocyc_com[-1]['formula'] = [line.strip()[18:]]
biocyc_com = {c['id']:c for c in biocyc_com}
biocyc_com1 = {c['id']:c for c in biocyc_com.values() if 'name' in c}
# 14513/14520
b_name1 = [x['name'] for x in biocyc_com.values() if 'name' in x]
b_name1 = [x.lower() for x in b_name1]





#   kegg

with open('E:/IGEM/atom/compound_raw.txt','rb') as file:
    f = file.readlines()
kegg_com = []
for i,line in enumerate(f):
    if kegg_com == [] and not line.startswith(b'ENTRY '):
        continue
    if line.startswith(b'ENTRY '):
        kegg_com.append({'id' : line.split()[1]})
    elif line.startswith(b'NAME '):
        name = [b' '.join(line.strip().split()[1:])]
        while name[-1][-1:] == b';':
            name[-1] = name[-1][:-1]
            i = i+1
            name.append(f[i].strip())
        kegg_com[-1]['name'] = name
    elif line.startswith(b'FORMULA '):
        kegg_com[-1]['formula'] = line.strip().split()[1]
kegg_com = {c['id']:c for c in kegg_com}
# 18293/18293
k_name1 = [x['name'] for x in kegg_com.values() if 'name' in x]
k_name1 = [[y.lower() for y in x] for x in k_name1]


with open('E:/IGEM/data/kegg_compound_name.txt') as file:
    f = file.readlines()
k_name_all = {x[4:10]:x.strip().split('\t')[1].split('; ') for x in f}


# lower, <i>, &alpha;
def name_process_1(name):
    name = re.sub('<.*?>','',name.lower())
    name = re.sub('\(-\)','l',name)
    name = re.sub('\(+\)','d',name)
    name = re.sub('amp','adenylate',name)
    name = re.split('[-\s,]', name)
    if '' in name:
        name.remove('')
    if name.count('l')>1:
        for i in range(name.count('l')-1):
            name.remove('l')
    if name.count('d')>1:
        for i in range(name.count('d')-1):
            name.remove('d')
    #if 'cis' in name:
    #    name.remove('cis')
    #if 'trans' in name:
    #    name.remove('trans')
    name = [re.sub('\W','',x) for x in name]
    name = [re.sub('_','',x) for x in name]
    #name = re.sub('&alpha;','alpha',name)
    #name = re.sub('&beta;','beta',name)
    return ''.join(name)

def name_process_2(name):
    name = re.sub('<.*?>','',name.lower())
    name = re.sub('\(-\)','l',name)
    name = re.sub('\(+\)','d',name)
    name = re.sub('amp','adenylate',name)
    name = re.split('[-\s,]', name)
    if '' in name:
        name.remove('')
    if 'l' in name:
        name.remove('l')
    if 'd' in name:
        name.remove('d')
    if 'cis' in name:
        name.remove('cis')
    if 'trans' in name:
        name.remove('trans')
    name = [re.sub('\W','',x) for x in name]
    name = [re.sub('_','',x) for x in name]
    return ''.join(name)

same_name_1 = []
same_name_2 = []
i = 0
for bn in biocyc_com1.keys():
    bname = biocyc_com1[bn]['name'].decode('gb2312')
    bname_1 = name_process_1(bname)
    for kn in k_name_all.keys():
        for kname in k_name_all[kn]:
            kname_1 = name_process_1(kname)
            if bname_1 == kname_1:
                same_name_1.append([bn.decode(),kn])
            else:
                bname_2 = name_process_2(bname)
                kname_2 = name_process_2(kname)
                if bname_2 == kname_2:
                    same_name_2.append([bn.decode(),kn])
    i += 1
    print(i)
    


corresponding_com = [tuple(x) for x in corresponding_com]
same_name_1 = [tuple(x) for x in same_name_1]
same_name_2 = [tuple(x) for x in same_name_2]

all_com = list(set(same_name_1).union(set(corresponding_com)))
# 7012
mismatch_com = list(set(corresponding_com).difference(set(same_name_1)))
# 1481
addition_com_1 = list(set(same_name_1).difference(set(corresponding_com)))
# 1931
addition_com_2 = list(set(same_name_2).difference(set(all_com)))
# 409


            
def print_com_name(i,com_list,end = 0):
    if end == 0:
        end = len(com_list)-1
    com_pair = com_list[i:end+1]
    for pair in com_pair:
        if not pair[0].encode() in biocyc_com or not pair[1].encode() in kegg_com:
            continue
            #print('not exist') 
        else:
            b_name = biocyc_com[pair[0].encode()]['name']
            k_name = kegg_com[pair[1].encode()]['name']
            print(pair)
            print('biocyc name: %s\nkegg_name:' % b_name, k_name , '\n')
            
            
 




with open('E:/IGEM/biocyc-kegg-compound-2.txt') as f:
    f = f.readlines()
com_list = [tuple(line.strip().split()) for line in f]
 


'''
correct_mismatch_com = []
for com_pair in mismatch_com1:     
    b_name = name_process(biocyc_com[com_pair[0].encode()]['name'])
    k_name = []
    for item in kegg_com[com_pair[1].encode()]['name']:
        k_name.append(name_process(item))
    
    if b_name in k_name:
        correct_mismatch_com.append(com_pair)
# 979 / 1933

mismatch_com2 = []
for item in mismatch_com1:
    if not item in correct_mismatch_com:
        mismatch_com2.append(item)
# 954
        


b_name3 = [x['name'] for x in biocyc_com.values() if 'name' in x]
b_name3 = [name_process(x) for x in b_name3]
k_name3 = [x['name'] for x in kegg_com.values() if 'name' in x]
k_name3 = [[name_process(y) for y in x] for x in k_name3]
b_list3 = [x for x in biocyc_com.values() if 'name' in x]
k_list3 = [x for x in kegg_com.values() if 'name' in x]

same_name3 = []
for com_pair in k_name3:
    for name in com_pair:
        if name in b_name3:
            same_name3.append(name)
            break
same_name3 = list(set(same_name3))
# same name: 5273
# no repeating: 5258
corresponding_com3 = []
for i,name in enumerate(same_name3):
    corresponding_com3.append([b_list[b_name3.index(name)]['id'].decode()])
    for item in k_list:
        for n in item['name']:
            if name_process(n) == name:
                corresponding_com3[-1].append(item['id'].decode())
                break
# 5258
                
new_com3 = []
corresponding_com3_line1 = [x[0] for x in corresponding_com3] # no repeating
for i,item in enumerate(corresponding_com3_line1):
    if item in corresponding_com1_line1:
        continue
    new_com3.append(corresponding_com3[i])
# 1723



final_com1 = corresponding_com + new_com3

def save_corresponding(path = 'E:/IGEM/biocyc-kegg-compound-1.txt'):
    with open(path,'w') as file:
        result = ''
        for item in final_com1:
            result += ' '.join(item) + '\n'
        file.write(result)


def save_smiles(path = 'E:/IGEM/biocyc-kegg-compound-smiles.txt'):
    with open(path,'w') as file:
        result = ''
        for item in final_com1:
            if item[0].encode() in biocyc_com and 'smiles' in biocyc_com[item[0].encode()]:
                result += item[0] + '\t' +  biocyc_com[item[0].encode()]['smiles'].decode() + '\n'
        file.write(result)







biocyc_id = [x[0] for x in final_com1]
kegg_id = [x[1] for x in final_com1]

'''






# reaction 
#   biocyc

with open('E:/IGEM/atom/data/reactions.dat', 'rb') as file:
    f = file.readlines()

biocyc_reaction=[]
for line in f:
    if line.startswith(b'UNIQUE-ID - '):
        biocyc_reaction.append([line.strip()[12:].decode()])
    elif line.startswith(b'EC-NUMBER'):
        ec_number = line.strip()[12:].decode()
        if ec_number[0] == '|' and ec_number[-1] == '|':
            ec_number = ec_number[1:-1]
        biocyc_reaction[-1].append(ec_number[3:])
    elif line.startswith(b'LEFT - '):
        biocyc_reaction[-1].append('left-' + line.strip()[7:].decode())
    elif line.startswith(b'RIGHT - '):
        biocyc_reaction[-1].append('right-' + line.strip()[8:].decode())


biocyc_without_enzyme = [x for x in biocyc_reaction if not x[1][0].isnumeric()]
biocyc_with_enzyme = [x for x in biocyc_reaction if x[1][0].isnumeric()]

biocyc_id2enz = {x[0]:x[1] for x in biocyc_with_enzyme}
biocyc_enz2id = {}
for item in biocyc_id2enz.keys():
    enz = biocyc_id2enz[item]
    if not enz in biocyc_enz2id:
        biocyc_enz2id[enz] = [item]
    else:
        biocyc_enz2id[enz].append(item)
        



# kegg


with open('E:/IGEM/atom/reaction_raw.txt','rb') as file:
    f = file.readlines()

kegg_reaction = []
for i,line in enumerate(f):
    if kegg_reaction == [] and not line.startswith(b'ENTRY '):
        continue
    if line.startswith(b'ENTRY '):
        kegg_reaction.append([line.split()[1].decode()])
    elif line.startswith(b'ENZYME'):
        kegg_reaction[-1].append(line.strip().split()[1:])
    elif line.startswith(b'EQUATION '):
        kegg_reaction[-1].append(b' '.join(line.split()[1:]))


kegg_without_enzyme = [x for x in kegg_reaction if len(x)==2]
kegg_with_enzyme = [x for x in kegg_reaction if len(x)==3]

kegg_id2enz = {x[0]:b' '.join(x[2]).decode() for x in kegg_with_enzyme}
kegg_enz2id = {}
for item in kegg_id2enz:
    enz_group = kegg_id2enz[item].split()
    for enz in enz_group:
        if not enz in kegg_enz2id:
            kegg_enz2id[enz] = [item]
        else:
            kegg_enz2id[enz].append(item)



inter_enzyme = sorted(set([enz for item in kegg_id2enz.values() for enz in item.split()])\
                      .intersection(set(biocyc_id2enz.values())))

common_enzyme_reaction = {enz:[biocyc_enz2id[enz], kegg_enz2id[enz]] for enz in inter_enzyme}



# get corresponding compounds

with open('E:/IGEM/biocyc-kegg-compound-1.txt') as f:
    f = f.readlines()
corresponding_com = [tuple(line.strip().split()) for line in f]
biocyc_com_id = [x[0] for x in corresponding_com]
kegg_com_id = [x[1] for x in corresponding_com]





biocyc_reaction_idx = {}

for i,item in enumerate(biocyc_reaction):
    com_list = [[],[]]
    for com in item[1:]:
        if com.startswith('left-'):
            com = com[5:]
            com_list[0].append(com)
        elif com.startswith('right-'):
            com = com[6:]
            com_list[1].append(com)

    biocyc_reaction_idx[item[0]] = com_list



kegg_reaction_idx = {}

for i,item in enumerate(kegg_reaction):
    com_list = [[],[]]
    com_rea, com_pro = item[1].decode().split('<=>')
    
    com_list[0] = re.findall('\D\d{5}',com_rea)
    com_list[1] = re.findall('\D\d{5}',com_pro)
    
    kegg_reaction_idx[item[0]] = com_list



def reaction_compare(b_id, k_id, reverse=False):
    if not b_id in biocyc_reaction_idx:
        #print('%s not in biocyc' % b_id)
        return [2]
    if not k_id in kegg_reaction_idx:
        #print('%s not in kegg' % k_id)
        return [3]

    b_rea_com, b_pro_com = biocyc_reaction_idx[b_id]
    if not reverse:
        k_rea_com, k_pro_com = kegg_reaction_idx[k_id]
    else:
        k_pro_com, k_rea_com = kegg_reaction_idx[k_id]
    # number compare
    if len(b_rea_com) != len(k_rea_com) or len(b_pro_com) != len(k_pro_com):
        # except PROTON / C00080
        b_rea_com = [x for x in b_rea_com if not x == 'PROTON']
        b_pro_com = [x for x in b_pro_com if not x == 'PROTON']        
        k_rea_com = [x for x in k_rea_com if not x == 'C00080']
        k_pro_com = [x for x in k_pro_com if not x == 'C00080']
        
        if len(b_rea_com) != len(k_rea_com) or len(b_pro_com) != len(k_pro_com):
            #print('number unequal')
            return [4]
        
    # compound compare
    b_result = [[-1]*len(b_rea_com), [-1]*len(b_pro_com)] # compared:1 / else:0
    k_result = [[-1]*len(k_rea_com), [-1]*len(k_pro_com)] # compared:1 / else:0
    
    for k_rea in k_rea_com:
        for b_rea in b_rea_com:
            if (b_rea, k_rea) in corresponding_com:
                if not b_rea_com.index(b_rea) in k_result[0] and not k_rea_com.index(k_rea) in b_result[0]:
                    k_result[0][k_rea_com.index(k_rea)] = b_rea_com.index(b_rea)
                    b_result[0][b_rea_com.index(b_rea)] = k_rea_com.index(k_rea)
                    break
    for k_pro in k_pro_com:
        for b_pro in b_pro_com:
            if (b_pro, k_pro) in corresponding_com:
                if not b_pro_com.index(b_pro) in k_result[1] and not k_pro_com.index(k_pro) in b_result[1]:
                    k_result[1][k_pro_com.index(k_pro)] = b_pro_com.index(b_pro)
                    b_result[1][b_pro_com.index(b_pro)] = k_pro_com.index(k_pro)
                    break
    
    b_result = [[1 if x>=0 else 0 for x in com_group] for com_group in b_result]
    k_result = [[1 if x>=0 else 0 for x in com_group] for com_group in k_result]
    return [1, b_result, k_result, [b_rea_com,b_pro_com], [k_rea_com,k_pro_com]]



# check corresponding reaction
def check_reaction(corresponding_rea):
    append_com = []
    undetermined_list = []
    rest_list =[]
    for item in corresponding_rea:
        b_id, k_id = item
        result = reaction_compare(b_id, k_id)
        if result[0] == 4 or result[0] == 1 and not 1 in [x for item in result[1] for x in item]:
            result = reaction_compare(b_id, k_id, reverse = True)
        
        if result[0] != 1 or not 0 in [x for item in result[1] for x in item]:
            continue
        
        if 0 in [x for item in result[1] for x in item]:
            b_rea_com, b_pro_com = result[3]
            k_rea_com, k_pro_com = result[4]
                
            if result[1][0].count(0) == 1:
                b_idx = result[1][0].index(0)
                k_idx = result[2][0].index(0)
                compound_pair = [b_rea_com[b_idx], k_rea_com[k_idx]]
                if not tuple(compound_pair) in append_com:
                    append_com.append(tuple(compound_pair))
            elif not item in undetermined_list:
                undetermined_list.append(item)
            if result[1][1].count(0) == 1:
                b_idx = result[1][1].index(0)
                k_idx = result[2][1].index(0)
                compound_pair = [b_pro_com[b_idx], k_pro_com[k_idx]]
                if not tuple(compound_pair) in append_com:
                    append_com.append(tuple(compound_pair))
            elif not item in undetermined_list:
                undetermined_list.append(item)
                
        if not item in append_com and not item in undetermined_list:
            rest_list.append(item)
    return [append_com, undetermined_list, rest_list]

i = 0
while True:
    flag = False
    append_com, undetermined_list, rest_list = check_reaction(corresponding_rea)
    for item in append_com:
        if item not in corresponding_com:
            corresponding_com.append(item)
            flag = True
    if not flag:
        break
    i += 1
    print('iter:',i)



'''
Triacylglycerols                C00422
DIACYLGLYCEROL                  C00165
3-Beta-Hydroxysterols           C02945
3-Beta-Hydroxysterol-Esters     C03587
Lysidine-tRNA-Ile2              C19723
GLN                             C00064
GLT                             C00025
CPD-9517                        C00251
DNA-Guanines                    C11475
Protein-S-methyl-L-cysteine     C03800
Ubiquitin-activating-protein-E1-L-cys       C00496
Ubiquitin-C-Terminal-Glycine    C02188
CPD-16184                       C07278
Glycoprotein-L-serine-or-L-threonine        C02189
Nucleoside-Monophosphates       C00215
Deoxy-Ribonucleosides           C02269
Thyroglobulin-triiodothyronines C16730
Thyroglobulin-aminoacrylates    C16734
Thyroglobulin-L-thyroxines      C16729
SUCROSE                         G00370
CPD-12462                       G01504
2-KETOGLUTARATE                 C00026
LysW-L-lysine                   C19889
CO-A                            C00040
ACYL-COA                        C00010
'''

'''
with open('E:/IGEM/biocyc-kegg-compound-append.txt') as f:
    f = f.readlines()
f = [tuple(x.strip().split()) for x in f]
corresponding_com = list(set(corresponding_com).union(set(f)))
'''



matched_reaction = {}
undetermined_list = {}
one_list = []
# corresponding reaction
for enz in common_enzyme_reaction.keys():
    matched_reaction[enz] = []
    undetermined_list[enz] = []
    biocyc_reaction_list, kegg_reaction_list = common_enzyme_reaction[enz]
    if len(biocyc_reaction_list) == 1 and len(kegg_reaction_list) == 1:
        one_list.append([biocyc_reaction_list[0], kegg_reaction_list[0]])
    for br in biocyc_reaction_list:
        for kr in kegg_reaction_list:
            result = reaction_compare(br, kr)
            if result[0] != 1:
                continue
            if not 1 in [x for y in result[1] for x in y]:
                result = reaction_compare(br, kr, reverse=True)
                if result[0] == 1 and not 0 in [x for y in result[1] for x in y]:
                    if not (br, kr) in matched_reaction.values():
                        matched_reaction[enz].append((br, kr))
                elif result[0] == 1 and 1 in [x for y in result[1] for x in y]:
                    if not (br, kr) in undetermined_list.values():
                        undetermined_list[enz].append((br, kr))
                else:
                    continue
            elif not 0 in [x for y in result[1] for x in y]:
                if not (br, kr) in matched_reaction.values():
                    matched_reaction[enz].append((br, kr))
            elif 1 in [x for y in result[1] for x in y]:
                if not (br, kr) in undetermined_list.values():
                    undetermined_list[enz].append((br, kr))

corresponding_rea2 = [tuple(rea) for item in list(matched_reaction.values()) for rea in item]
corresponding_rea = list(set([tuple(x) for x in corresponding_rea + corresponding_rea2]))



corresponding_rea = list(set(corresponding_rea).union(set([tuple(x) for x in one_list])))





def check_duplicate(rea_com_list, path = ''):
    result = ''
    line_1 = [x[0] for x in rea_com_list]
    line_2 = [x[1] for x in rea_com_list]
    for name in set(line_1):
        if line_1.count(name) > 1:
            if path == '':
                print([x for x in rea_com_list if x[0] == name])
            else:
                result += '\t'.join([x[0] for x in rea_com_list if x[0] == name])
    for name in line_2:
        if line_2.count(name) > 1:
            if path == '':
                print([x for x in rea_com_list if x[1] == name])
    
'''
result = ''
for item in corresponding_rea:
    result +=' '.join(item) + '\n'
with open('E:/IGEM/biocyc-kegg-reaction-2.txt','w') as f:
    f.write(result)
'''
'''
result = ''
for item in corresponding_com:
    result +=' '.join(item) + '\n'
with open('E:/IGEM/biocyc-kegg-compound-2.txt','w') as f:
    f.write(result)
'''


'''
result = []
with open('E:/IGEM/biocyc-kegg-reaction-2.txt') as f:
    f = f.readlines()
for line in f:
    result.append(line.strip())
'''



# search example
'''
import pickle
file = open('E:/IGEM/data/atom/atom_transfer.pkl','rb')
f = pickle.load(file)


num = 0
for item in f:
    if len(item[2])<15 and len(item[3])<15 and len(item[2])>5 and len(item[3])>5:
        num +=1
    if num>N:
        break
print(item)
N += 1
'''