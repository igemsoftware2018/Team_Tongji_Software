# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 15:56:17 2018

@author: Administrator
"""

import re
import pickle


atom_file = 'data//atom_transfer.pkl'
file = open(atom_file,'rb')
atom_transfer_list = pickle.load(file)


def get_transfer_data(com_pair, k_rea_id, atom_transfer_list):
    for item in atom_transfer_list:
        if item[0] == com_pair and item[1] == k_rea_id:
            return [item[2],item[3]]


def compare_transfer_data(data1, data2):
    if len(data1[1]) != len(data2[0]):
        return []
    new_data2 = [[],[]]
    new_data2[0] = [x if data1[1][i] > 0 else -1 for i,x in enumerate(data2[0])]
    new_data2[1] = [x if x in new_data2[0] else -1 for x in data2[1]]
    return new_data2


def atom_conservation_score(compound_list, reaction_list):
    compound_list = compound_list.split('→')
    compound_list = [x+'_'+compound_list[i+1] for i,x in enumerate(compound_list) if i!= len(compound_list)-1]
    reaction_list = reaction_list.split('→')
    print(compound_list)
    print(reaction_list)
    
    
    atom_c_list = list(set([x[0] for x in atom_transfer_list]))
    atom_r_list = list(set([x[1] for x in atom_transfer_list]))
    for coms in compound_list:
        if not coms in atom_c_list:
            print('compound error', coms)
            return -1
    for rea in reaction_list:
        if not rea in atom_r_list:
            print('reaction error', rea)
            return -1
    
    for i,coms in enumerate(compound_list):
        if i == len(compound_list)-1:
            break
        if coms.split('_')[1] != compound_list[i+1].split('_')[0]:
            print('error in %s and %s' % (compound_list[i], compound_list[i+1]))
            return -1

    transfer_list = []
    for com_pair, k_rea_id in zip(compound_list, reaction_list):
        transfer_data = get_transfer_data(com_pair, k_rea_id, atom_transfer_list)
        if transfer_list != []:
            compare_data = compare_transfer_data(transfer_list[-1], transfer_data)
            if compare_data == []:
                print('mismatch')
            transfer_list.append(compare_data)
        else:
            transfer_list.append(transfer_data)
    score = len([x for x in transfer_list[-1][1] if x != -1]) / len(transfer_list[0][0])
    return score

        