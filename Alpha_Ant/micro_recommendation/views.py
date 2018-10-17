from django.shortcuts import render_to_response
from micro_recommendation.forms import recommendation_Form
from micro_recommendation.compound_path import main



def recommendation(request):
    form=recommendation_Form
    return render_to_response('micro_recommendation.html',{'form':form})

def recommendation_results(request):
    if request.method=='GET':
        form =recommendation_Form(request.GET)
        if form.is_valid():
            data=form.cleaned_data
            result=main(data)
            return render_to_response('recom_results.html',{'result':result})