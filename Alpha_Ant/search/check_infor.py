# -*- coding: utf-8 -*-


# 检查表单输入的信息是否符合本地数据库标准

import os
os.chdir('C:\\Users\\27364\\Desktop\\software 3')

import pickle


c_list=[] #compound编号信息列表
c_name=[] #compound中文名称信息表
organ=[]  #物种信息以及缩写全称信息表


with open('data//nodes3.csv','r') as f:
    c_list=[line.strip().split(',')[0] for line in f]
    
with open('data//compound_name_pair.csv','r') as f:
    for line in f:
        c_name+=line.strip().split(',')[1:]

with open('data//org_short.pickle','rb') as f:
    organ_short=pickle.load(f)
        
for key in organ_short:
    organ.append(key)
    short_list=organ_short[key].split(',')
    organ+=short_list

 
