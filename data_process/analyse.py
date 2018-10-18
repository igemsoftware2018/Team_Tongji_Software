# -*- coding: utf-8 -*-
"""
Created on Sat May  5 20:57:24 2018

@author: Administrator
"""
import re
import requests
'''import matplotlib.pyplot as plt'''


#by_product ---- H2O, ATP, NAD+, NADH, NADPH, NADP+, O2
#           ---- ADP, P, CoA, CO2, PPi, NH3, UDP
#           ---- FAD, H2O2, Acceptor, Reduced acceptor, GDP , GTP, e-
#           ---- H+, H2, Reduced NADPH, Oxidized NADPH
by_product = ['C00001','C00002','C00003','C00004','C00005','C00006','C00007',\
              'C00008','C00009','C00010','C00011','C00013','C00014','C00015',\
              'C00016','C00027','C00028','C00030','C00035','C00044','C05359',\
              'C00080','C00282','C03024','C03161']


def get_num(i):
    if isinstance(i,str):
        return i[1:]
    if i<10:
        result = '0000' + str(i)
    elif i<100:
        result = '000' + str(i)
    elif i<1000:
        result = '00' + str(i)
    elif i<10000:
        result = '0' + str(i)
    else:
        result = str(i)
    return result


def get_r_info(i):
    r = requests.get('http://rest.kegg.jp/get/R' + get_num(i))
    if r.text == '':
        return []
    info = {'id':'R' + get_num(i)}
    r = r.text.split('\n')
    for i,item in enumerate(r):
        if item.startswith('EQUATION'):
            info['EQUATION'] = ' '.join(item.split()[1:]).replace('<=>','=')
        elif item.startswith('RCLASS'):
            info['RCLASS'] = item.split()[2:]
            while r[i+1].startswith('  '):
                info['RCLASS'] += r[i+1].split()[1:]
                i += 1
        elif item.startswith('ENZYME'):
            info['ENZYME'] = item.split()[1:]
    return info


def get_c_info(i):
    r = requests.get('http://rest.kegg.jp/get/C' + get_num(i))
    if r.text == '':
        return []
    info = {'id':'C' + get_num(i)}
    for item in r.text.split('\n'):
        if item.startswith('FORMULA'):
            info['FORMULA'] = item.split()[-1]
        elif item.startswith('NAME'):
            info['NAME'] = ' '.join(item.split()[1:])
        elif item.startswith('MOL_WEIGHT'):
            info['WEIGHT'] = item.split()[-1]
    return info


def enz_from_file(line):
    pattern = '\d\.\d'
    line = line.strip().split(',')
    enz = []
    for item in line:
        if re.findall(pattern,item) != []:
            enz.append(item)
    return enz


def com_from_file(line):
    pattern = 'C\d{5}'
    line = line.strip().split(',')
    for item in line:
        compounds = re.findall(pattern,item)
        if compounds != []:
            return compounds
    return [None]


def g_from_file(line):
    pattern = 'G\d{5}'
    line = line.strip().split(',')
    for item in line:
        compounds = re.findall(pattern,item)
        if compounds != []:
            return compounds
    return [None]
'''
g = {}
for line in f:
    result = g_from_file(line)
    if result != [None]:
        for item in result:
            if item not in g:
                g[item] = 1
            else:
                g[item] += 1
'''
    

def eq_to_pairs(eq):
    pattern = 'C\d{5}'
    Rea, Pro = eq.split('=')
    rea = re.findall(pattern, Rea)
    pro = re.findall(pattern, Pro)
    
    pairs = []
    for r in rea:
        for p in pro:
            pairs.append(r+'_'+p)
    
    coe = []    
    for r in rea:
        Rea = re.sub(r,'',Rea)
    for p in pro:
        Pro = re.sub(p,'',Pro)
    coe.append([int(x) if x.strip()!= '' else 1 for x in Rea.split('+')])
    coe.append([int(x) if x.strip()!= '' else 1 for x in Pro.split('+')])
    
    return pairs,coe


def r_main_info(a,b):
    result = ''
    for i in range(a,b):
        line = ''
        info = get_r_info(i)
        if info == []:
            continue
        if 'id' in info:
            line += info['id'] + ', '
        if 'EQUATION' in info:
            line += info['EQUATION'] + ', '
        if 'ENZYME' in info:
            line += ', '.join(info['ENZYME']) + ', '
        if 'RCLASS'  in info:
            line += ', '.join(info['RCLASS'])
        else:
            line = '*' + line + ', '.join(eq_to_pairs(info['EQUATION'])[0])
        result += line + '\n'
    return result



def get_polymer():
    polymer = []
    pattern = '\(.*?\)n'
    for i in range(1,21832):
        c_info = get_c_info(i)
        if 'FORMULA' in c_info and re.findall(pattern,c_info['FORMULA']) != []:
            print(i,c_info)
            polymer.append(c_info['id'])
    out = open('polymer.txt','w')
    out.write('\n'.join(polymer))
    out.close()
'''
get_polymer()
'''
    
#34  38  3184  6925  12570  13427
 

def polymer_list():
    file = open('E://IGEM//polymer.txt')
    f = file.readlines()
    file.close()
    
    p_list = []
    for line in f:
        p_list.append(line.strip())
    return p_list



def del_byproduct():
    pattern = 'C\d{5}_C\d{5}'
    result = ''
    file = open('reaction.csv')
    f = file.readlines()
    file.close()
    for line in f:
        if re.findall(pattern,line) == []:
            result += line
            continue
        if 'RCLASS' in get_r_info(int(line[1:6])):
            result += line
            continue
        
        line = line.strip().split(',')
        for i,item in enumerate(line):
            if re.findall(pattern,item) != []:
                pairs = line[i:]
                result += ','.join(line[:i])
                break
        for pair in pairs:
            if pair == '':
                break
            pair = pair.split('_')
            if pair[0] not in by_product and pair[1] not in by_product:
                result += ',' + '_'.join(pair)
            else:
                print('del',line[0])
        result += '\n'
    out = open('reaction_main_pairs.csv','w')
    out.write(result)
    out.close()
'''
del_byproduct()                
'''



def get_file_info(file_name, info_type):
    file = open(file_name)
    f = file.readlines()
    enzyme = {}
    compound = {}
    g_compound = {}
    g_reaction = []
    polymer_reaction = []
    
    for line in f:
        e = enz_from_file(line)
        c = com_from_file(line)
        g = g_from_file(line)
        for item in e:
            if item not in enzyme:
                enzyme[item] = 1
            else:
                enzyme[item] += 1
        for item in c:
            if item not in compound:
                compound[item] = 1
            else:
                compound[item]+= 1
            
        if g != [None]:
            g_reaction.append(line.split(',')[0])
        for item in g:
            if item not in g_compound:
                g_compound[item] = 1
            else:
                g_compound[item]+= 1
        
    if info_type == 1:
        return enzyme,compound
    elif info_type == 2:
        return g_compound,g_reaction
        
'''
enzyme,compound = get_file_info('reaction.csv',1)

sort_enzyme = sorted(enzyme.items(),key=lambda x:x[1],reverse=True)
sort_compound = sorted(compound.items(),key=lambda x:x[1],reverse=True)

print('enzyme num: %i,  max: %i' % (len(enzyme), max(enzyme.items(),key=lambda x:x[1])[1]))
print('compound num: %i,  max: %i' % (len(compound), max(compound.items(),key=lambda x:x[1])[1]))


e_y = [0] * 172
for item in list(enzyme.items()):
    e_y[item[1]-1] += 1
plt.ylim(0,50)
plt.plot(range(1,173),e_y)

c_y = [0] * 3807
for item in list(compound.items()):
    c_y[item[1]-1] += 1
plt.ylim(0,50)
plt.plot(range(1,3808),c_y)
'''



def exist_main_weight(rea_cw,pro_cw):
    W1 = max(rea_cw)
    W2 = max(pro_cw)
    if len([w for w in rea_cw if w>W1/2]) == 1 and len([w for w in pro_cw if w>W2/2]) == 1 and abs(W1-W2)<W1/3:
        return True
    
    
def find_main_pair(rea,pro,coe):
    main_pair = [''] * 2
    main_weight = [0] * 2
    message = ''
    rea_w = [get_c_info(c) for c in rea]
    rea_w = [float(c['WEIGHT']) if 'WEIGHT' in c else -1 for c in rea_w]
    rea_cw = [c*w for c,w in zip(coe[0],rea_w)]
    pro_w = [get_c_info(c) for c in pro]
    pro_w = [float(c['WEIGHT']) if 'WEIGHT' in c else -1 for c in pro_w]    
    pro_cw = [c*w for c,w in zip(coe[1],pro_w)]
    print(rea_w, pro_w, coe, rea_cw, pro_cw, end = '\t')
    
    if len(rea) == 1 and len(pro) == 1:
        message = 'only one'
        main_pair = [rea[0],pro[0]]
        main_weight = [rea_cw[0],pro_cw[0]]
        print(main_pair,main_weight, message)
        return main_pair, main_weight, message
    
    if len(set(rea) - set(by_product)) == 1:
        main_pair[0] = list(set(rea) - set(by_product))[0]
        main_weight[0] = rea_cw[rea.index(main_pair[0])]
    if len(set(pro) - set(by_product)) == 1:
        main_pair[1] = list(set(pro) - set(by_product))[0]
        main_weight[1] = pro_cw[pro.index(main_pair[1])]
    
    if main_pair[0] != '' and main_pair[1] != '':
        message = 'curtain r & p'
        print(main_pair,main_weight, message)
        return main_pair, main_weight, message
    elif main_pair[0] != '' and main_pair[1] == '':
        min_D = 1000
        closest_w = 0
        closest_c = ''
        for i,w in enumerate(pro_cw):
            if abs(main_weight[0]-w) < min_D:
                min_D = abs(main_weight[0]-w)
                closest_w = w
                closest_c = pro[i]
        main_pair[1] = closest_c
        main_weight[1] = closest_w
        message = 'curtain r'
    elif main_pair[1] != '' and main_pair[0] == '':
        min_D = 1000
        closest_w = 0
        closest_c = ''
        for i,w in enumerate(rea_cw):
            if abs(main_weight[1]-w) < min_D:
                min_D = abs(main_weight[1]-w)
                closest_w = w
                closest_c = rea[i]
        main_pair[0] = closest_c
        main_weight[0] = closest_w
        message = 'curtain p'
    else:
        if -1 in rea + pro:
            message = 'uncertain weight'
        elif exist_main_weight(rea_cw,pro_cw):
            main_weight[0] = max(rea_cw)
            main_weight[1] = max(pro_cw)
            main_pair[0] = rea[rea_cw.index(main_weight[0])]
            main_pair[1] = pro[pro_cw.index(main_weight[1])]
            message = 'main weight of: %0.2f & %0.2f' % (main_weight[0],main_weight[1]) 
        else:
            message = 'no main weight: ' + ' '.join([str(x) for x in rea_cw]) + ' = ' + ' '.join([str(x) for x in pro_cw])
    print(main_pair,main_weight, message)
    return main_pair, main_weight, message
    
    
    
#none on rclass
def reaction_process(file_name,out_name,n):
    file = open(file_name)
    f = file.readlines()
    file.close()
    
    if n == -1:
        n = len(f)
    
    not_rc = []
    processed_id = []
    for line in f[:n]:
        if line[0] == '*':
            not_rc.append(line[1:7])
    result = []
    for R in not_rc:
        eq = get_r_info(R)['EQUATION']
        if 'n' in eq:
            print(R + '\n' + 'n in eq')
            result.append(R + ',n in eq')
            continue
        if 'G' in eq:
            print(R + '\n' + 'G in eq')
            result.append(R + ',G in eq')
            continue
        processed_id.append(R)
        pattern = 'C\d{5}'
        rea, pro = eq.split('=')
        rea = re.findall(pattern, rea)
        pro = re.findall(pattern, pro)
        coe = eq_to_pairs(eq)[1]
        print(R + '\n' + eq)
        pair, weight, message = find_main_pair(rea,pro,coe)
        result.append(R + ',' + '_'.join(pair) + ',' + message + ',' + ','.join([str(x) for x in weight]))
    
    print(len(processed_id), len(result))
    out = open(out_name,'w')
    out.write('\n'.join(result))
    out.close()


def save_file(n):
    result = ''
    for i in range(1,n):
        print(i)
        line = ''
        info = get_r_info(i)
        if info == []:
            continue
        if 'RCLASS' not in info:
            line += '*'
        line += info['id']
        line += ',' + info['EQUATION']
        if 'ENZYME' in info:
            for e in info['ENZYME']:
                line += ',' + e
        if line[0] != '*':
            for r in info['RCLASS']:
                line += ',' + r
        result += line + '\n'
    return result
    

    
    