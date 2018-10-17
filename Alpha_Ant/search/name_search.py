import re
import pickle

with open('data//kegg_compound_name.pkl','rb') as f:
    name_dict = pickle.load(f)
with open('data//kegg_to_compound_name.pkl','rb') as f:
    id_dict = pickle.load(f)


MAX_NUM = 10

def start_with(pattern, text):
    text = text.lower()
    pattern = pattern.lower()
    
    text = re.sub('\W', '', text)
    text = re.sub('\d', '', text)
    pattern = re.sub('\W', '', pattern)
    pattern = re.sub('\d', '', pattern)
    
    if pattern != '' and text.startswith(pattern):
        return 1
    return 0
    
    
def is_kegg_id(pattern):
    if pattern[0].lower() == 'c' and pattern[1:6].isdigit() and len(pattern) == 6:
        if pattern.upper() in id_dict:
            return 1
    return 0

    
def search_name_id(pattern):
    num = 0
    result=''
    if pattern=='':
        return '$'
    else:
        for com in name_dict.keys():
            if is_kegg_id(pattern):
                result=id_dict[pattern.upper()].split(';')[0]+'$'+result
                return result[:-1]
                break
            if start_with(pattern, com):
                result+=com+'$'
                num += 1
                if num >= MAX_NUM:
                    return result[:-1]
                    break;
        return result

