# -*- coding: utf-8 -*-


def tran_route(c_list):
    result=''
    for c in c_list:
        result+=c+'→'
    return result[:-1]

def data_clean(data,need_list,not_need_list):
    
    #exp:data[0]
    # ['C00082', 'C00022', 'C00024', 'C00223'],361.4,['L-Tyrosine','Pyruvate','Acetyl-CoA','p-Coumaroyl-CoA']

    result=[]
    
    for item in data:
        
        new_rout=item[0] # ['C00082', 'C00022', 'C00024', 'C00223']
        compound_name=item[2] # ['L-Tyrosine','Pyruvate','Acetyl-CoA','p-Coumaroyl-CoA']
        new_compound_name=''
        for name in compound_name:
            new_compound_name+=name+' | '
        item[2]=new_compound_name.strip()[0:-1].strip()
        
        item[2],item[1]=item[1],item[2]
        #print(item)
        if need_list==['']:
            for c in new_rout:
                if c in not_need_list:
                    break
                if c==new_rout[-1]:
                    item[0]=tran_route(item[0])
                    result.append(item)
        elif need_list!='':
            
            for c in new_rout:
                if c in not_need_list:
                    break
                if c in need_list:
                    need_list.remove(c)
                if c==new_rout[-1] and len(need_list)==0:
                    item[0]=tran_route(item[0])
                    result.append(item)
    return result

def trans_compound(com_list):
    result=[]
    for c_pair in com_list:
        c_list=c_pair.split('_')
        for c in c_list:
            if c not in result:
                result.append(c)
    return result
    # ['C00082→C00022→C00024→C00223', 'L-Tyrosine | Pyruvate | Acetyl-CoA | p-Coumaroyl-CoA', 361.4]
if __name__=='__main__':
    data=[[['C00082', 'C00022', 'C00024', 'C00223'],361.4,['L-Tyrosine','Pyruvate','Acetyl-CoA','p-Coumaroyl-CoA']]]
    print(data_clean(data,[''],['']))