# -*- coding: utf-8 -*-


import pickle
#import os

#os.chdir('E:\Igem\web-software\software\data')

with open('data//enzyme_info.pkl','rb') as f:
    enzyme_infor=pickle.load(f)

R_E={}
with open('data//R_E.csv','r') as f:
    for line in f:
        line=line.strip().split(',')
        R_E[line[0]]=line[1].split(' ')[0]
        
def enzyme_information(r):
    enzyme=R_E[r]
    print(enzyme)
    infor=enzyme_infor[enzyme]
    infor=infor.split('\n')
    result='<p>'+'<b>'+'Reaction ID: '+r+'</br>'
   
    result+="<a href='/download/gene/'>"+infor[0]+'</a>'+'&nbsp;&nbsp;&nbsp;(click and acquire gene information)'+'</br>'
    result+=infor[1]+'</br>'
    #result+=infor[2]+'</br>'
    result+='</b>+</p><p>'
    infor=infor[2:]
    for line in infor:
        if '#' in line:
            line=line.replace('#','')
            #line.replace('#','&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;')
            result+='</br>'+'<b>'+line+'</b>'+'</br>'
        elif '#' not in line:
            line.replace('#','&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;')
            result+=line+'</br>'
    result+='</p>'
    '''
    infor=infor.replace('\n','</br>')
    infor=infor.replace('\t','&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;')
    '''
    return result

if __name__=='__main__':
    print(enzyme_information('R00005'))