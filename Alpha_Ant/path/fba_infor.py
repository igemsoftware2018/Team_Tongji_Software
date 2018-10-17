from __future__ import print_function
import cobra
import cobra.test
import pandas
import pickle
from time import time
from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)
#import os
#os.chdir('C:\\Users\\27364\\Desktop\\software 3')

with open('data//bigg_reaction.pkl','rb') as f:
    bigg_reactions=pickle.load(f)
with open('data//fba_compound.pkl','rb') as f:
    fba_compound=pickle.load(f)
model = cobra.test.create_test_model("ecoli")

def FDA(kegg_reactions):    ###kegg_reactions = ['R00509']
    model_reactions = []
    model_metabolites = []
    for item1 in model.reactions:
        model_reactions.append(item1.id)
    for item2 in model.metabolites:
        model_metabolites.append(item2.id)
    for k_id in kegg_reactions:
        if k_id in bigg_reactions.keys():
            if bigg_reactions[k_id]['id'] not in model_reactions:
               reaction = Reaction(bigg_reactions[k_id]['id'])
               print(reaction.id)
               reaction.name = ''
               reaction.subsystem = ''
               reaction.lower_bound = bigg_reactions[k_id]['lower_bound']  # This is the default
               reaction.upper_bound = bigg_reactions[k_id]['upper_bound']# This is the default
               for key in bigg_reactions[k_id]['metabolites']:
                   if key in model_metabolites:
                       reaction.add_metabolites({
                           model.metabolites.get_by_id(key): bigg_reactions[k_id]['metabolites'][key]
                       })
                   else :
                       a = Metabolite(key,formula='', name='',compartment='')
                       reaction.add_metabolites({
                           a: bigg_reactions[k_id]['metabolites'][key]
                       })
               reaction.gene_reaction_rule = bigg_reactions[k_id]['gene_reaction_rule']
               model.add_reactions([reaction])
###指定一个或一系列目标反应

###下面是执行的一些功能

def summary_data(a):
    #处理总结数据
    a=a.replace('-','').replace(')',')\n\n').replace('CONSUMING','\n\nCONSUMING').replace('REACTION\n      ','REACTION\n').replace(' ','#')  #a的每一行代表一行数据
    fba_list=[]
    for key in fba_compound:
        if key in a:
            fba_list.append(key)
    fba_list.sort(key=lambda x:len(x),reverse=True)
    for fba in fba_list:
        if fba in a:
            a=a.replace(fba,fba_compound[fba])
    a=a.replace('\n','</br>').replace('#','&nbsp;')
    

    return a

def fba(reaction):
    reaction_list = reaction.split('→')
    FDA(reaction_list)
    '''
    result=[]
    for reaction in reaction_list:
        if reaction in bigg_reactions:
            ID=bigg_reactions[reaction]['id']
            result.append([reaction,ID,model.optimize().fluxes[ID]])
            '''
    s=''
    # if len(result)!=0:
      
    #   for line in result:
    #     s+=str(line[0])+'   :   '+str(line[1])+'   :   '+str(line[2])+'</br>'
    #     s='<p>'+s+'</p>'
    
      #s+='<p>'+model.summary().replace('\n','</br>')+'</p>'
    s+='<p>'+summary_data(model.metabolites.nadh_c.summary())+'</p>'
    s+='<p>'+'</br></br>'+summary_data(model.metabolites.atp_c.summary())+'</p>'
    return s
if __name__=='__main__':
    fba('R00728→R00210→R08767→R01613→R02446')
'''
model.optimize().fluxes
model.summary()
model.metabolites.nadh_c.summary()
'''

