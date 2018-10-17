from django.shortcuts import HttpResponse
from django.utils.safestring import mark_safe
from django.shortcuts import render_to_response
from django.views.decorators.csrf import csrf_exempt

from django.http import HttpResponseRedirect

from results.compound_path import main
from search.forms import search_form
import json
import pickle

# 读取化合物信息
com_infor={}
with open('data//compound_data.pkl','rb') as f:
  com_infor=pickle.load(f)

def results(request):
    if request.method=='GET':
        form =search_form(request.GET)
        
        if form.is_valid():
            data=form.cleaned_data    
            #result=main(data,data1)
            #return render_to_response('results.html',{'result':result})
            return render_to_response('results.html',{'data':mark_safe(data)})
    else:
        return render_to_response('main.html')

@csrf_exempt    
def ajax_test(request):
    if request.method=='POST':
        data=json.dumps(request.POST)
        mydict=json.loads(data)
        result=main(mydict)
        jarray=trans_graph_data(result)
        output=''
        for item in result:
            path=''
            for row in item:
                path+=str(row)+'$$'
            output+=path[:-2]+'&&'
        output=output[:-2]+'@@@'+jarray
        return HttpResponse(output)

def trans_graph_data(data):
    if 10<len(data):  #图数据最多展示10条
        data=data[:10]
    top_end=data[0][0].split('→')[0]+'-'+data[0][0].split('→')[-1]
    jarray=''
    for item in data:
        c_list=item[0].split('→')
        for i in range(len(c_list)-1):
            c=c_list[i]+'-'+c_list[i+1]
            if c not in jarray:

                jarray+=c+'$$'
    jarray=top_end+'$$'+jarray
    return jarray[:-2]
   
def graph_data(request):
    if request.method=='GET':
        jarray = [
  {"source": "C00002", "target": "C00004", "type": "resolved"},
  {"source": "C00002", "target": "C00003", "type": "resolved"},
  {"source": "C00002", "target": "C00005", "type": "resolved"},
  {"source": "C00002", "target": "C00006", "type": "resolved"},
  {"source": "C00002", "target": "C00007", "type": "resolved"},
  {"source": "C00002", "target": "C00008", "type": "resolved"},
  {"source": "C00002", "target": "C00009", "type": "resolved"},
  {"source": "C00003", "target": "C00010", "type": "resolved"},
  {"source": "C00005", "target": "C00010", "type": "resolved"},
  {"source": "C00011", "target": "C00010", "type": "resolved"}
] 
    print(json.dumps(jarray))
    return HttpResponse(json.dumps(jarray), content_type="application/json")

@csrf_exempt 
def get_compound_infor(request):
  #点击球球，出现化合物的详细信息
  if request.method=="POST":
    c_number=request.POST.get('number')
    infor=com_infor[c_number]
    result=''
    for l in infor:
      result+=str(l)+'$'
    return HttpResponse(result[:-1])