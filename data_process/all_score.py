# -*- coding: utf-8 -*-


# 给定物种输出分数

import pickle
import math
import csv

#读取各物种内源酶文件
with open('bac_enzyme.pkl','rb') as f:
    bact_enzyme=pickle.load(f)

#读取各个反应吉布斯分数初步处理数据
data=[]
with open('R_E_C_G_S.csv','r') as f:
    for line in f:
        line=line.strip().split(',')
        data.append(line)    

def if_enzyme(e_list,enzyme):
    # native reactions or native reactions
    
    for e in e_list:
        if e in enzyme:
            return True
    return False

def c_score(c,direction,enzyme):
    weight_score=0
    for line in data:
       
        e_list=line[1].strip().split(',')
        if if_enzyme(e_list,enzyme):
            c_pair=line[2].strip().split(' ')
            if direction=='Positive':
                for C in c_pair:
                    C=C.split('_')[0]
                    if C==c:
                        score=float(line[4])
                        if score>=700:
                            score=700
                        elif score<=-700:
                            score=-700
                        e_exp=math.exp(-score)
                        weight_score+=e_exp
            elif direction=='Reverse':
                for C in c_pair:
                    C=C.split('_')[1]
                    if C==c:
                        score=float(line[4])
                        if score>=700:
                            score=700
                        elif score<=-700:
                            score=-700
                        e_exp=math.exp(score)
                        weight_score+=e_exp
    return weight_score

def trans_C(C):
    c=C.split('_')
    output=c[1]+'_'+c[0]
    return output

def main(bact):
    bact_name=bact
    bact=bact_enzyme[bact]
    data1=[]
    for line in data:
        e_list=line[1].strip().split(',')
        c_pair=line[2].strip().split(' ')
        score=float(line[4])
        weight_score=''
        if score>=700:
            score=700
        elif score<=-700:
            score=-700
        e_exp=math.exp(-score)
        if if_enzyme(e_list,bact):
            for c in c_pair:
                c=c.split('_')[0]
                foreign_score=c_score(c,direction='Positive',enzyme=bact)
                if e_exp/(1+foreign_score) ==0:
                    weight_score+=str(-1000)+' '
                else:
                    weight_score+=str(math.log(e_exp/(1+foreign_score)))+' '
            e_exp=1/e_exp
            for c in c_pair:
                c=c.split('_')[1]
                foreign_score=c_score(c,direction='Reverse',enzyme=bact)
                if e_exp/(1+foreign_score) ==0:
                    weight_score+=str(-1000)+' '
                else:
                    weight_score+=str(math.log(e_exp/(1+foreign_score)))+' '
        else:
            for c in c_pair:
                c=c.split('_')[0]
                foreign_score=c_score(c,direction='Positive',enzyme=bact)
                if e_exp/(1+foreign_score+e_exp) ==0:
                    weight_score+=str(-1000)+' '
                else:
                    weight_score+=str(math.log(e_exp/(1+foreign_score+e_exp)))+' '
            e_exp=1/e_exp
            for c in c_pair:
                c=c.split('_')[1]
                foreign_score=c_score(c,direction='Reverse',enzyme=bact)
                if e_exp/(1+foreign_score+e_exp) ==0:
                    weight_score+=str(-1000)+' '
                else:
                    weight_score+=str(math.log(e_exp/(1+foreign_score+e_exp)))+' '
        line.append(weight_score)
        data1.append(line)
        print(len(data1),line[5])
        
    for i in range(len(data1)):
        data1[i][5]=data1[i][5].strip().replace(' ','_')
    
    cond={}
    for line in data1:
        R=line[0]
        C_list=line[2].strip().split(' ')
        S_list=line[5].split('_')
        l=len(C_list)
        for i in range(l):
            C=C_list[i]
            cond[(R,C)]=round(float(S_list[i]),2)
            cond[(R,trans_C(C))]=round(float(S_list[i+l]),2)
            
    with open('pkl\\'+bact_name+'.pkl','wb') as f:
        pickle.dump(cond,f)
        
    with open('EXCEL\\'+bact_name+'.csv','w',newline='') as f:
        csv_writer=csv.writer(f)
        csv_writer.writerows(data1)
        
if __name__=='__main__':
    for key in bact_enzyme:
        main('Acaryochloris marina')
    