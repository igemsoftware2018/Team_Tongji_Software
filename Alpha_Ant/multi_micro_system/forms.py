# -*- coding: utf-8 -*-


from django import forms

class InputForm(forms.Form):
    smiles = forms.CharField(initial='CCCO')   
    MaxLength = forms.IntegerField(initial=10)