# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 19:44:47 2018

@author: Administrator
"""
import pickle
'''
with open('E:/IGEM/data/new_bac_enzyme.pickle','rb') as f:
    org_dict = pickle.load(f)

with open('E:/IGEM/data/short_bact.csv') as f:
    org_short = f.readlines()

org_short = [x.strip().split(',') for x in org_short]


full2short = {}
for item in org_short:
    if item[1] in org_dict:
        if not item[1] in full2short:
            full2short[item[1]] = item[0]
        else:
            full2short[item[1]] += ',' + item[0]
    

with open('E:/IGEM/data/org_short.pickle','wb') as f:
    pickle.dump(full2short, f)
'''

MAX = 5

with open('data//org_short.pickle','rb') as f:
    full2short = pickle.load(f)
short2full = {}
for full in full2short.keys():
    for short in full2short[full].split(','):
        short2full[short] = full


def start_with(pattern, text):
    pattern = pattern.lower()
    text = text.strip().lower()
    
    for item in text.split():
        if item.startswith(pattern):
            return 1
    return 0


def search_org(pattern):
    result = ''
    num = 0
    if len(pattern)<=3:
        for short_name in short2full.keys():
            if start_with(pattern, short_name):
                result += short_name + ': ' + short2full[short_name] + '$'
                num += 1
                if num == MAX:
                    return result[:-1]
    for full_name in full2short.keys():
        if start_with(pattern, full_name):
            result += full_name + '$'
            num += 1
            if num == MAX:
                return result[:-1]
    return result[:-1]



'''
import time
start = time.time()
name = 'coli'
for i in range(1,len(name)+1):
    print(search_org(name[:i]))
end = time.time()
print(end-start)
'''