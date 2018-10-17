import pickle
import csv
from results.data_processing import data_clean
from results.data_processing import trans_compound
from results.atom_conservation import atom_conservation_score
#import os
#os.chdir('E:\Igem\web-software\software-加速版')

visited_nodes = dict()
neighbor_dict = dict()
enzyme_dict = dict()
compound_dict = dict()
gibbs_dict = dict()
toxicity_dict = dict()
max_eco = dict()
freq_dict = dict()

with open("data//nodes3.csv") as nodes_file:
    node_reader = csv.reader(nodes_file)
    for node in node_reader:
        visited_nodes.update({node[0]: False})

with open('data//adj_map2.pk', 'rb') as f:
    neighbor_dict = pickle.load(f)

with open('data//reaction_dictionary.pk', 'rb') as f:
    reaction_dict = pickle.load(f)

with open('data//max_eco.pkl', 'rb') as f:
    max_eco = pickle.load(f)

with open('data//frequency.pkl', 'rb') as f:
    freq_dict = pickle.load(f)

with open("data//reaction_enzyme_pair.csv") as enzyme_file:
    emzyme_reader = csv.reader(enzyme_file)
    for e in emzyme_reader:
        enzyme_dict.update({e[0]: e[1:]})

with open("data//compound_name_pair.csv") as c_file:
    c_reader = csv.reader(c_file)
    for c in c_reader:
        compound_dict.update({c[0]: c[1:]})
    compound_dict.update({'':''})

with open('data//Gibbs.pkl', 'rb') as g_f:
    gibbs_dict = pickle.load(g_f)
    
with open('data//eco_Toxicity.csv','r') as f:
    for line in f:
        line=line.strip().split(',')
        toxicity_dict.update({line[0]:round(float(line[1]),2)})


def dfs(start_compound, target_compound, depth=10):

    path_list = []

    cur_depth = 0
    path_stack = []

    visited_nodes[start_compound] = True
    path_stack.append([start_compound, 0])
    while len(path_stack) > 0:
        if cur_depth >= depth:
            temp_top = path_stack.pop()
            cur_depth -= 1
            visited_nodes[temp_top[0]] = False
            continue
        cur_compound = path_stack[-1]
        temp_adj = neighbor_dict[cur_compound[0]]
        # for i in graph.neighbors(cur_compound[0]):
        #     temp_adj.append(i)

        length = len(temp_adj)
        if cur_compound[1] >= length:
            temp_top = path_stack.pop()
            visited_nodes[temp_top[0]] = False
            cur_depth -= 1
            continue
        else:
            next_compound = temp_adj[cur_compound[1]]
            path_stack[-1][1] += 1

            if next_compound == target_compound:
                path_stack.append([next_compound, 0])
                path_list.append(path_stack[:])
                path_stack.pop()
            else:
                if not visited_nodes[next_compound]:
                    visited_nodes[next_compound] = True
                    path_stack.append([next_compound, 0])
                    cur_depth += 1
    path_result = []
    for path in path_list:
        temp = []
        for item in path:
            temp.append(item[0])
        path_result.append(temp)

    return path_result


def all_dfs(start_compound, depth=3):

    path_list = []

    cur_depth = 0
    path_stack = []

    visited_nodes[start_compound] = True
    path_stack.append([start_compound, 0])
    while len(path_stack) > 0:

        if cur_depth >= depth:
            temp_top = path_stack.pop()
            cur_depth -= 1
            visited_nodes[temp_top[0]] = False
            continue
        cur_compound = path_stack[-1]
        temp_adj = neighbor_dict[cur_compound[0]]

        length = len(temp_adj)
        if cur_compound[1] >= length:
            temp_top = path_stack.pop()
            visited_nodes[temp_top[0]] = False
            cur_depth -= 1
            continue
        else:
            next_compound = temp_adj[cur_compound[1]]
            path_stack[-1][1] += 1

            if not visited_nodes[next_compound]:
                visited_nodes[next_compound] = True
                path_stack.append([next_compound, 0])
                path_list.append(path_stack[:])
                cur_depth += 1

    path_result = []
    for path in path_list:
        temp = []
        for item in path:
            temp.append(item[0])
        path_result.append(temp)

    return path_result

def rank_list(path_result, w_gibbs = 1, w_toxicity = 1, w_frequency = 1):
    #[['C00002', 'C00008'], ['C00002', 'C00020', 'C00008'], ['C00002', 'C00020', 'C00498', 'C00008']]
    rank_result = []
    for com_list in path_result:
        num = len(com_list)
        weight = 0.0
        for i in range(0, num - 1):
            tempstr = com_list[i] + "_" + com_list[i + 1]
            weight += w_gibbs*float(max_eco[tempstr])
            if freq_dict.get(tempstr) != None:
                weight += w_frequency*float(freq_dict[tempstr])
            else:
                weight+=w_frequency*(-400)
        for item in com_list:
            if toxicity_dict.get(item) != None:
                weight += w_toxicity*toxicity_dict[item]
        weight = round(weight, 2)
        rank_result.append([com_list,weight])

    sort_result = sorted(rank_result,key=lambda w: w[1], reverse=True)

    return sort_result
    '''
    if conservation >= len(rank_result):
        return sort_result
    else:
        return sort_result[:conservation]
    '''


def re_rank_list(path_result, w_gibbs = 1, w_toxicity = 1, w_frequency = 1):

    rank_result = []
    for com_list in path_result:
        num = len(com_list)
        print(com_list)
        com_list = list(reversed(com_list))
        print(com_list)
        weight = 0.0
        for i in range(0, num - 1):
            tempstr = com_list[i] + "_" + com_list[i + 1]
            weight += w_gibbs*float(max_eco[tempstr])
            if freq_dict.get(tempstr) != None:
                weight += w_frequency*float(freq_dict[tempstr])
            else:
                weight+=w_frequency*(-400)
        for item in com_list:
            if toxicity_dict.get(item) != None:
                weight += w_toxicity*toxicity_dict[item]
        weight = round(weight, 2)
        rank_result.append([com_list,weight])

    sort_result = sorted(rank_result,key=lambda w: w[1], reverse=True)

    return sort_result
    '''
    if conservation >= len(rank_result):
        return sort_result
    else:
        return sort_result[:conservation]
    '''


def attach_inform(sort_result):
    inform_result = []
    for row in sort_result:  # (('R01417', 'R02421'), ['C00002_C00008', 'C00008_C00498'], 0.0)
        cnames = []
        for com in row[0]:
            cnames.append(compound_dict[com][0])

        row.append(cnames)
        inform_result.append(row)

    return inform_result


def reverse_info(inform):
    #[('R01417', 'R02421'), ['C00002_C00008', 'C00008_C00498'], 0.0, [['2.7.3.10'], ['2.4.1.21', '2.4.1.242']], [['ATP', 'ADP'], ['ADP', 'ADP-glucose']]]
    reversed = []
    for row in inform:
        temp_re = []
        for i in range(len(row[0])-1,-1,-1):
            # print(len(row[0]))
            temp_re.append(row[0][i])
        temp_comp = []
        for i in range(len(row[1])-1,-1,-1):
            s = row[1][i].split('_')
            temp_comp.append(s[1] + "_" + s[0])
        temp_enzyme = []
        for i in range(len(row[4])-1,-1,-1):
            temp_enzyme.append(row[4][i])
        temp_name = []
        for i in range(len(row[5])-1,-1,-1):
            temp_name.append([row[5][i][1],row[5][i][0]])
        reversed.append([temp_re,temp_comp,row[2],row[3],temp_enzyme,temp_name])

    return reversed;

def simple_path(beg, end, depth,w_gibbs = 1, w_toxicity = 1, w_frequency = 1):

    

    compound_result = dfs(beg, end, int(depth))

    sort_result = rank_list(compound_result, w_gibbs, w_toxicity, w_frequency)
    print(sort_result)
    final_result = attach_inform(sort_result)
    print(final_result)
    return final_result


def all_path(comp, depth, w_gibbs = 1, w_toxicity = 1, w_frequency = 1):

    
    compound_result = all_dfs(comp, int(depth))

    sort_result = rank_list(compound_result, w_gibbs, w_toxicity, w_frequency)
    print(sort_result)
    final_result = attach_inform(sort_result)
    print(final_result)
    return final_result


def reverse_all_path(comp, depth, w_gibbs = 1, w_toxicity = 1, w_frequency = 1):
    
    compound_result = all_dfs(comp, int(depth))

    sort_result = re_rank_list(compound_result, w_gibbs, w_toxicity, w_frequency)
    print(sort_result)
    final_result = attach_inform(sort_result)
    print(final_result)
    return final_result

def trans_C(c):
    #转中文名为编号
    for key,val in compound_dict.items():
        if key==c:
            return c
        elif c in val:
            return key
        
def trans_list(l):
    for i in range(len(l)):
        l[i]=trans_C(l[i])
    return l

def main(condon):
    start_compound = trans_C(condon['Input'])    #起始化合物
    target_compound = trans_C(condon['Output'])  #终止化合物
    
    depth = condon['MaxLength']                        #搜索深度
    conservation = int(condon['result_conservation'])  #结果保留量

    #  取权值
    w_gibbs = float(condon['Gibbs'])
    w_toxicity = float(condon['Toxicity'])
    w_frequency = float(condon['Frequency'])

    #  结果筛选
    required_compound=trans_list(condon['requrired'].strip().split(','))
    not_required_compound=trans_list(condon['not_requrired'].strip().split(','))
    
    result=[]
    if start_compound!='' and target_compound!='':
        result = simple_path(start_compound, target_compound, depth,w_gibbs,w_toxicity,w_frequency)
    elif start_compound!='' and target_compound=='':
        result = all_path(start_compound,depth,w_gibbs,w_toxicity,w_frequency)
    elif start_compound=='' and target_compound!='':
        result = reverse_all_path(target_compound, depth,w_gibbs,w_toxicity,w_frequency)

    result=data_clean(result,required_compound,not_required_compound)

    if conservation < len(result):
        result = result[:conservation]
    # ['C00082→C00022→C00024→C00223', 'L-Tyrosine | Pyruvate | Acetyl-CoA | p-Coumaroyl-CoA', 361.4]

    '''
    if conservation < len(result):
        result = result[:conservation]
    '''

    return result



if __name__ == '__main__':
    data=main(condon)

