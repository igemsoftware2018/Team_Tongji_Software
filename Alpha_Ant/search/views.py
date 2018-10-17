from django.shortcuts import render
from django.shortcuts import render_to_response
from django.shortcuts import HttpResponse
from search.forms import search_form
from django.views.decorators.csrf import csrf_exempt
from search.name_search import search_name_id
from search.org_name_search import search_org
#from software.compound_path import main

def search(request):
    form=search_form()
    return render_to_response('search_form.html',{'form':form})

def main(request):
    return render_to_response('main.html')

#模糊匹配
'''
with open('data//nodes3.csv','r') as f:
    com_number_list=[line.strip().split(',')[0] for line in f]

com_name_list=[]
with open('data//compound_name_pair.csv','r') as f:
    for line in f:
        line=line.strip().split(',')[1:]
        for name in line:
            com_name_list.append(name)
            
all_list=com_name_list+com_number_list

def match(input_data):
    input_data=input_data.upper()
    match_list=[]
    for name in all_list:
        if input_data==name.upper()[0:len(input_data)]:
            match_list.insert(0,name)
        elif input_data in name.upper():
            match_list.append(name)
    if len(match_list)>=10:
        match_list=match_list[0:10]
    result=''
    for match in match_list:
    	result+=match+'$'
    return result[:-1]
'''
@csrf_exempt
def input_ajax(request):
	#input表单的实时提示
	 if request.method=='POST':
		  input_data=request.POST.get('input')  #用户输入的化合物
		  result=search_name_id(input_data)
		  return HttpResponse(result)
@csrf_exempt
def organ_ajax(request):
    if request.method=='POST':
        input_data=request.POST.get('input')  #用户输入的化合物
        result=search_org(input_data)
        return HttpResponse(result)