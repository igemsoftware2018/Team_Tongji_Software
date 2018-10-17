from django.shortcuts import render_to_response
from micro_recommendation.forms import recommendation_Form
from micro_recommendation.compound_path import main
import pickle

with open('data//organ_herf.pkl','rb') as f:
    organ_herf=pickle.load(f)


def recommendation(request):
    form=recommendation_Form
    return render_to_response('micro_recommendation.html',{'form':form})

def recommendation_results(request):
    if request.method=='GET':
        form =recommendation_Form(request.GET)
        if form.is_valid():
            data=form.cleaned_data
            result=main(data)
            for i in range(len(result)):
                organ=result[i][0]
                if organ in organ_herf:
                    result[i]+=(organ_herf[organ],)
                elif organ.split(' ')[0] in organ_herf:
                    result[i]+=(organ_herf[organ.split(' ')[0]],)
            return render_to_response('recom_results.html',{'result':result})