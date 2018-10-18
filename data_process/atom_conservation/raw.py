import requests
import numpy as np


def get_compounds():
    c_list = []
    for i in range(1,21832):
        id = str(i)
        r = requests.get('http://rest.kegg.jp/get/C' + '0'*(5-len(id)) + id)
        c_list.append(r.text)
        if i%100 == 0:
            print(i)
    np.save('C:/raw_c',np.array(c_list))


def get_reactions():
    r_list = []
    for i in range(1,11991):
        id = str(i)
        r = requests.get('http://rest.kegg.jp/get/R' + '0'*(5-len(id)) + id)
        r_list.append(r.text)
        if i%100 == 0:
            print(i)
    np.save('C:/raw_r',np.array(r_list))