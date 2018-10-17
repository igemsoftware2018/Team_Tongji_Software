from django.shortcuts import render
from django.shortcuts import render_to_response
from django.shortcuts import HttpResponse
from search.forms import search_form
from path.pathway_information import main
from path.get_enzyme_information import enzyme_information
from path.fba_infor import fba
from django.views.decorators.csrf import csrf_exempt
from path.atom_score import atom_conservation_score

#from software.compound_path import main

def FDA_result(request,reaction):
    data=main(reaction)
    return render_to_response('fda.html',{'result':data})

def path_information(request,c_list):
    data=main(c_list)
    result=[]
    for i in range(len(data)):   
        s=round(atom_conservation_score(c_list,data[i][0]),2)
        if s!=-1:
        	s=str(s*100)+'%'
        d = data[i][:]
        d.insert(2,s)
        result.append(d)
    print(result)
    c_list=c_list.split('→')
    graph_data=''
    for line in data:
        r_list=line[0].split('→')
        for i in range(len(r_list)):
            l=c_list[i]+'$'+c_list[i+1]+'$'+r_list[i]
            if l not in graph_data:
                graph_data+=l+'&'
    graph_data=graph_data[:-1]
    return render_to_response('path.html',{'result':result,'graph_data':graph_data})

@csrf_exempt
def reaction_information(request):
    if request.method=="POST":
        r=request.POST.get('reaction')
        result=enzyme_information(r)
        return HttpResponse(result)

@csrf_exempt
def get_fba(request):
    if request.method=="POST":
        r=request.POST.get('r')
        result=fba(r)
        return HttpResponse(result)

def download_gene(request):
    file=open('data//gibbs.csv','rb')  
    response =HttpResponse(file)
    response['Content-Type']='application/octet-stream'  
    response['Content-Disposition']='attachment;filename="gibbs.csv"'  
    return response