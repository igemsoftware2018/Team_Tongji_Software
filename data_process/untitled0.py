# -*- coding: utf-8 -*-
"""
Created on Sun Aug 19 20:41:32 2018

@author: Administrator
"""
import re
path = 'D:/Documents/WeChat Files/blizzardsky7/Files/'

with open(path + 'full_database_1.xml','rb') as f:
    f = f.readlines()
    
drugtype = []
text = b''
name = b''
for line in f:
    if line.startswith(b'<drug type='):
        if name != b'':
            drugtype.append([name,text])
        text = b''
        name = re.findall(b'<drug type="(.*?)"', line)
        if len(name) > 1:
            print('error')
            break
        name = name[0]
    else:
        text += line
        
drugtype.append([name,text])



result = b''
for item in drugtype:
    result += item[0] + b'\n'
    result += item[1]
    
with open('E:/drug.txt','wb') as f:
    f.write(result)