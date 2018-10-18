# -*- coding: utf-8 -*-

import pickle
import numpy as np

main_path = 'E:/IGEM/'

with open(main_path + 'data/frequency/R_E_C_G_S.csv') as file:
    f = file.readlines()
enz_file = open(main_path + 'data/frequency/bac_enzyme.pkl','rb')
org2enz = pickle.load(enz_file)

enz2pair = {}
org2pair = {}
pair_freq = {}


non_essential_pair = []
non_existent_enz = []


for line in f:
    enz_group, pairs = line.split(',')[1:3]
    enz_group = enz_group.strip().split()
    pairs = pairs.strip().split()
    
    if enz_group == []:
        non_essential_pair += pairs
        
    for enz in enz_group:
        for pair in pairs:
            if enz in enz2pair:
                enz2pair[enz].append(pair)
            else:
                enz2pair[enz] = [pair]


for org in org2enz:
    pair_group = []
    for enz in org2enz[org]:
        if enz in enz2pair:
            pair_group += enz2pair[enz]
        else:
            non_existent_enz.append(enz)
    pair_group = list(set(pair_group))
    
    org2pair[org] = pair_group
        



for pair_group in org2pair.values():
    for pair in pair_group:
        if pair in pair_freq:
            pair_freq[pair] += 1
        else:
            pair_freq[pair] = 1





# norm
pair_freq = sorted(pair_freq.items(), key=lambda x:x[1])
freq_list = np.array([int(x[1]) for x in pair_freq])
freq_list = np.log(freq_list)
std = freq_list.std()
avg = freq_list.mean()


freq_list = (freq_list-max(freq_list))/std * 100
pair_freq = {item[0]:item[1] for item in zip([x[0] for x in pair_freq],freq_list)}


file = open(main_path + 'data/frequency/frequency.pkl', 'wb')
pickle.dump(pair_freq, file)








import pickle
with open(main_path + 'data/frequency/frequency.pkl', 'rb') as f:
    pair_freq = pickle.load(f)

with open(main_path + 'data/frequency/frequency.pkl', 'wb') as f:
    pickle.dump(pair_freq, f)
    
    
    
    