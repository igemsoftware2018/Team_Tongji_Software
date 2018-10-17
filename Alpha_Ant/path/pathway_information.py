# -*- coding: utf-8 -*-

import pickle
import itertools

#import os
#os.chdir('E:\Igem\web-software\software-加速版')


rea_dict_file = ('data//reaction_dictionary.pk')
with open(rea_dict_file, 'rb') as f:
    rea_dict = pickle.load(f)

with open('data//Gibbs.pkl', 'rb') as g_f:
    gibbs_dict = pickle.load(g_f)

eco_dic={}
with open('data//ECO.pkl','rb') as f:
    eco_dic=pickle.load(f)

name_dic={}
with open('data//compound_name_pair.csv','r') as f:
    for line in f:
        line=line.strip().split(',')
        name_dic[line[0]]=line[1]
        
        
# com_list: split by →
def com_to_rea(com_list):
    # 将反应物链的所有可能列出
    path = com_list.split('→')
    name=''
    for com in path:
        name+=name_dic[com]+'|'
    name=name[:-1]
    full_path=[]

    temp = []
    num = len(path)
    for i in range(0, num - 1):
        temp.append(str(path[i]) + '_' + str(path[i + 1]))
    # temp: C1_C2 C2_C3....
    rea_set = []
    for item in temp:
        rea_set.append(rea_dict[item])
    first, *second = rea_set
    temp_rea = itertools.product(first, *second)
    full_path=[temp, list(temp_rea)]
    full_path.append(name)
     
    return full_path

def rank_list(full_path):
    #[['C00002_C00008', 'C00008_C00020'], [('R00002', 'R00122'), ('R00002', 'R00127'), ('R00002', 'R00157')
    rank_result = []
    name=full_path[-1]
    full_path=full_path[:-1]
    com_list=full_path[0]
    
    for rec_list in full_path[1]:  #('R00002', 'R00122'
        weight = 0.0
        gb = 0.0
        for i in range(len(rec_list)): # R00002
            weight += float(eco_dic[(rec_list[i], com_list[i])])

            gb += gibbs_dict[(rec_list[i], com_list[i])]
            gb=round(gb,2)
            weight=round(weight,2)
        rank_result.append((rec_list,name,gb,weight))

    sort_result = sorted(rank_result, key=lambda w: w[3], reverse=True)
    return sort_result

def main(com_lsit):
    full_path=com_to_rea(com_lsit)
    result=list(rank_list(full_path))
    for i in range(len(result)):
        result[i]=list(result[i])
        r=''
        for j in range(len(result[i][0])):
            r+=result[i][0][j]+'→'
        r=r[:-1]
        result[i][0]=r
    return result

if __name__=='__main__':
    a=main('C00082→C00022→C00024→C00223→C06561→C00509')