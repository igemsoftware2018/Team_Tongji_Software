# -*- coding: utf-8 -*-


from django import forms

class recommendation_Form(forms.Form):   
    Input = forms.CharField(required=False)
    Output = forms.CharField(required=False)
    MaxLength = forms.IntegerField()
    result_conservation = forms.IntegerField(initial=50)
    
