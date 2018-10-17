import pickle
import csv
import itertools
import os

visited_nodes = dict()
neighbor_dict = dict()
enzyme_dict = dict()
compound_dict = dict()
gibbs_dict = dict()
species_list = os.listdir('data//all_oragn_score')
species_score = []
full_path = []


with open("data//nodes3.csv") as nodes_file:
    node_reader = csv.reader(nodes_file)
    for node in node_reader:
        visited_nodes.update({node[0]: False})

with open('data//adj_map2.pk', 'rb') as f:
    neighbor_dict = pickle.load(f)

with open('data//reaction_dictionary.pk', 'rb') as f:
    reaction_dict = pickle.load(f)

with open('data//ECO.pkl', 'rb') as f:
    eco = pickle.load(f)

with open("data//reaction_enzyme_pair.csv") as enzyme_file:
    emzyme_reader = csv.reader(enzyme_file)
    for e in emzyme_reader:
        enzyme_dict.update({e[0]: e[1:]})

with open("data//compound_name_pair.csv") as c_file:
    c_reader = csv.reader(c_file)
    for c in c_reader:
        compound_dict.update({c[0]: c[1]})

with open('data//Gibbs.pkl', 'rb') as g_f:
    gibbs_dict = pickle.load(g_f)


def load_data():

    for spe in species_list:
        global species_score
        with open('data//all_oragn_score//'+spe, 'rb') as f:
            species_score.append(pickle.load(f))


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

    return path_list


def get_compound_pair_list(path_list):

    # path_list :  C1 C2 C3...
    for path in path_list:
        temp = []
        num = len(path)
        for i in range(0, num - 1):
            temp.append(str(path[i][0]) + '_' + str(path[i + 1][0]))
        # temp: C1_C2 C2_C3....
        rea_set = []
        for item in temp:
            rea_set.append(reaction_dict[item])
        first, *second = rea_set
        temp_rea = itertools.product(first, *second)
        global full_path
        full_path.append([temp, list(temp_rea)])

    return full_path


def rank_list(conservation, species, w_gibbs = 0, w_in_ex = 0, w_posion = 0):
    #[['C00002_C00008', 'C00008_C00020'], [('R00002', 'R00122'), ('R00002', 'R00127'), ('R00002', 'R00157')
    rank_result = []
    global full_path
    for com_list in full_path:
        for rec_list in com_list[1]:  #('R00002', 'R00122'
            weight = 0
            gb = 0.0
            for i in range(len(rec_list)): # R00002
                weight += float(species_score[species][(rec_list[i], com_list[0][i])])
                # gb += gibbs_dict[(rec_list[i], com_list[0][i])]
            rank_result.append(weight)

    sort_result = sorted(rank_result, reverse=True)
    if conservation >= len(rank_result):
        return average(sort_result)
    else:
        return average(sort_result[:conservation])


def attach_inform(sort_result):
    inform_result = []
    for row in sort_result:  # (('R01417', 'R02421'), ['C00002_C00008', 'C00008_C00498'], 0.0)
        enzymes = []
        for rection in row[0]:
            enzymes.append(enzyme_dict[rection])
        cnames = []
        for com in row[1]:
            temp = com.split('_')
            cnames.append([compound_dict[temp[0]],compound_dict[temp[1]]])

        # gb = 0.0
        # for i in range(0,len(row[0])):
        #     gb += gibbs_dict[row[0][i],row[1][i]]

        row = list(row)
        row.append(enzymes)
        row.append(cnames)
        # row.append(gb)
        inform_result.append(row)

    return inform_result


def simple_path(beg, end, depth, conserv):
    load_data()
    result = dict()
    compound_result = dfs(beg, end, int(depth))

    print('search done!')
    #路径格式转换
    temp_ans = get_compound_pair_list(compound_result)

    print('convert done!')
    #搜索排序
    for i in range(len(species_list)):
        result.update({species_list[i].replace('.pkl',''): round(rank_list(int(conserv), i),2)})
    print('sort done!')

    return result


def average(list):
    avg = 0.0
    avg = sum(list)/(len(list)*1.0)
    return avg


def main(condon):
    global full_path
    full_path = []
    start_compound = condon['Input']
    target_compound = condon['Output']
    depth = condon['MaxLength']
    conservation = condon['result_conservation']


    result = simple_path(start_compound, target_compound, depth, conservation)
    result=sorted(result.items(),key=lambda x:x[1],reverse=True)
    
    '''
    for item in result.keys():
        data.append([item,result[item]])
    '''
    return result
    '''
    for item in result.keys():
        print(item)
        print(result[item])
    '''


if __name__ == '__main__':
    main()
