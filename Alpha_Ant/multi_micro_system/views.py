from django.shortcuts import render_to_response
from multi_micro_system.forms import InputForm
from multi_micro_system.smiles import get_similar_compound

#from software.compound_path import main

def multi_system(request):
    form=InputForm()
    #form1=ScoreForm()
    return render_to_response('multi_micro_system.html',{'form':form})

def multi_sysytem_results(request):
    if request.method=='GET':
        form =InputForm(request.GET)
        
        if form.is_valid():
            data=form.cleaned_data
            result=get_similar_compound(data)
            return render_to_response('multi_results.html',{'result':result})
    else:
        return render_to_response('main.html')