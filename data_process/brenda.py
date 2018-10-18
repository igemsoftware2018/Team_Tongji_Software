import re
import pickle

main_path = 'E:/IGEM/'
with open(main_path + 'data/brenda/brenda_download.txt', 'rb') as f:
    f = f.readlines()

head_list = []

for line in f:
    if line[:2].isalpha() and line[:2].isupper() and line[2] == 9:
        head = line[:2]
        if not head in head_list:
            head_list.append(head)



brenda_dict = {}

for i,line in enumerate(f):
    if line.startswith(b'ID'):
        if brenda_dict != {}:
            brenda_dict[enz]['ORG'] = org_list
            brenda_dict[enz]['RN'] = rn_list
            brenda_dict[enz]['CF'] = cf_list
            brenda_dict[enz]['KM'] = km_list
            brenda_dict[enz]['PHR'] = phr_list
            brenda_dict[enz]['PHO'] = pho_list
            brenda_dict[enz]['TO'] = to_list
            brenda_dict[enz]['TR'] = tr_list
            brenda_dict[enz]['LO'] = lo_list
            
        enz = line.strip().split()[1].decode()
        brenda_dict[enz] = {}
        org_list = []
        rn_list = []
        cf_list = []
        km_list =[]
        phr_list = []
        pho_list = []
        to_list = []
        tr_list = []
        lo_list = []            
    
    
    elif line.startswith(b'PR\t'):
        while(not f[i+1].startswith(b'PR\t') and f[i+1] != b'\n'):
            line = line[:-1] + f[i+1]
            i+=1
        org = re.sub('\(.*?\)', '', line.decode()).strip()
        org = re.findall('(#[\d\t,]*?#.*?)<', org)[0]
        org_list.append(org)
        
        
    elif line.startswith(b'RN\t'):
        while(not f[i+1].startswith(b'PR\t') and f[i+1] != b'\n'):
            line = line[:-1] + f[i+1]
            i+=1
        rn = re.sub('RN\t', '', line.decode()).strip()
        rn_list.append(rn)
        
        
    elif line.startswith(b'CF\t'):
        while(not f[i+1].startswith(b'CF\t') and f[i+1] != b'\n'):
            line = line[:-1] + f[i+1]
            i+=1
        #
        org = re.findall('#([\d\t,]*?)#.*?', line.decode())[0].strip()
        org = re.split('[\t,]', org)
        cf = re.findall('#[\d\t,]*?#(.*?)[<\(]', line.decode())[0].strip()
        cf_list = [org, cf]
        
        
    elif line.startswith(b'KM\t'):
        while(not f[i+1].startswith(b'KM\t') and f[i+1] != b'\n'):
            line = line[:-1] + f[i+1]
            i+=1
        km = re.findall('(#[\d\t,]*?#.*?})', line.decode())[0].strip()
        km_list.append(km)
        
        
    elif line.startswith(b'PHR\t'):
        while(not f[i+1].startswith(b'PHR\t') and f[i+1] != b'\n'):
            line = line[:-1] + f[i+1]
            i+=1
        phr = re.sub('\(.*?\)', '', line.decode())
        phr = re.findall('(#[\d\t,]*?#.*?)<', phr)[0].strip()
        phr_list.append(phr)
        
        
    elif line.startswith(b'PHO\t'):
        while(not f[i+1].startswith(b'PHO\t') and f[i+1] != b'\n'):
            line = line[:-1] + f[i+1]
            i+=1
        
        pho = re.sub('\(.*?\)', '', line.decode())
        pho = re.findall('(#[\d\t,]*?#.*?)<', pho)[0].strip()
        pho_list.append(pho)
        
        
    elif line.startswith(b'TO\t'):
        while(not f[i+1].startswith(b'TO\t') and f[i+1] != b'\n'):
            line = line[:-1] + f[i+1]
            i+=1
        to = re.sub('\(.*?\)', '', line.decode())
        to = re.findall('(#[\d\t,]*?#.*?)<', to)[0].strip()
        to_list.append(to)


    elif line.startswith(b'TR\t'):
        while(not f[i+1].startswith(b'TR\t') and f[i+1] != b'\n'):
            line = line[:-1] + f[i+1]
            i+=1
        tr = re.sub('\(.*?\)', '', line.decode())
        tr = re.findall('(#[\d\t,]*?#.*?)<', tr)[0].strip()
        tr_list.append(tr)   


    elif line.startswith(b'LO\t'):
        while(not f[i+1].startswith(b'LO\t') and f[i+1] != b'\n'):
            line = line[:-1] + f[i+1]
            i+=1
        lo = re.sub('\(.*?\)', '', line.decode())
        lo = re.findall('(#[\d\t,]*?#.*?)<', lo)[0].strip()
        lo_list.append(lo)   
        
        
brenda_dict[enz]['ORG'] = org_list
brenda_dict[enz]['RN'] = rn_list
brenda_dict[enz]['CF'] = cf_list
brenda_dict[enz]['KM'] = km_list
brenda_dict[enz]['PHR'] = phr_list
brenda_dict[enz]['PHO'] = pho_list
brenda_dict[enz]['TO'] = to_list
brenda_dict[enz]['TR'] = tr_list
brenda_dict[enz]['LO'] = lo_list



save_path = main_path + 'data/brenda/brenda.pkl'
with open(save_path,'wb') as save_file:
    pickle.dump(brenda_dict, save_file)
    
    
    

'''
import pickle
with open('E:/IGEM/data/brenda/brenda.pkl','rb') as save_file:
    brenda_dict = pickle.load(save_file)
'''



brenda_info = {}

for enz in brenda_dict.keys():
    info_text = 'Enzyme ID: %s\nEnzyme name: %s\n' % (enz, brenda_dict[enz]['RN'][0])
    
    org_dict = {}
    for org in brenda_dict[enz]['ORG']:
        org_idx = org.split()[0][1:-1]
        org_name = ' '.join(org.strip().split()[1:])
        org_dict[org_idx] = org_name
    
    
    info_text += '#Localization:\n'
    lo_info = []
    for lo in brenda_dict[enz]['LO']:
        for org in org_dict:
            if org in lo.strip().split()[0]:
                lo_info.append(org_dict[org] + ': ' + ' '.join(lo.strip().split()[1:]) + '\n')
    info_text += ''.join(sorted(lo_info))


    info_text += '#PH Range:\n'
    phr_info = []
    for phr in brenda_dict[enz]['PHR']:
        for org in org_dict:
            if org in phr.strip().split()[0]:
                phr_info.append(org_dict[org] + ': ' + ' '.join(phr.strip().split()[1:]) + '\n')
    info_text += ''.join(sorted(phr_info))
                
    
    info_text += '#PH Optimum:\n'
    pho_info = []
    for pho in brenda_dict[enz]['PHO']:
        for org in org_dict:
            if org in pho.strip().split()[0]:
                pho_info.append(org_dict[org] + ': ' + ' '.join(pho.strip().split()[1:]) + '\n')
    info_text += ''.join(sorted(pho_info))
    
    
    info_text += '#Temperature Range:\n'
    tr_info = []
    for tr in brenda_dict[enz]['TR']:
        for org in org_dict:
            if org in tr.strip().split()[0]:
                tr_info.append(org_dict[org] + ': ' + ' '.join(tr.strip().split()[1:]) + '\n')
    info_text += ''.join(sorted(tr_info))
                
                
    info_text += '#Temperature Optimum:\n'
    to_info = []
    for to in brenda_dict[enz]['TO']:
        for org in org_dict:
            if org in to.strip().split()[0]:
                to_info.append(org_dict[org] + ': ' + ' '.join(to.strip().split()[1:]) + '\n')
    info_text += ''.join(sorted(to_info))
                
    
    info_text += '#KM:\n'       
    km_info = []     
    for km in brenda_dict[enz]['KM']:
        for org in org_dict:
            if org in km.strip().split()[0]:
                km_info.append(org_dict[org] + ': ' + ' '.join(km.strip().split()[1:]) + '\n')
    info_text += ''.join(sorted(km_info))
    
    brenda_info[enz] = info_text




with open('E:/IGEM/data/brenda/enzyme_info.pkl','wb') as f:
    pickle.dump(brenda_info, f)






