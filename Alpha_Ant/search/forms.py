# -*- coding: utf-8 -*-


from django import forms

class search_form(forms.Form):
    Microorganism = forms.CharField(initial='ECO')   
    Input = forms.CharField(required=False)
    Output = forms.CharField(required=False)
    MaxLength = forms.IntegerField()
    result_conservation = forms.IntegerField(initial=50)
    
    Gibbs=forms.FloatField(initial=1)
    Toxicity=forms.FloatField(initial=0.1)
    Frequency=forms.FloatField(initial=0.3)
    requrired = forms.CharField(required=False)
    not_requrired = forms.CharField(required=False)