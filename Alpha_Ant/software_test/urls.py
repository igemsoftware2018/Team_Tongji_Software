"""software_test URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/2.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path
from django.conf.urls import url
from search.views import search,main,input_ajax,organ_ajax
from results.views import results,ajax_test,graph_data,get_compound_infor
from path.views import path_information,reaction_information,get_fba,download_gene
from micro_recommendation.views import recommendation,recommendation_results
from multi_micro_system.views import multi_system,multi_sysytem_results
#from FDA.views import FDA_result

urlpatterns = [
    path('admin/', admin.site.urls),
    url(r'^search/$',search),
    url(r'^results/$',results),
    url(r'^recommendation/$',recommendation),
    url(r'^recom_results/$',recommendation_results),
    url(r'^multi_sysytem/$',multi_system),
    url(r'^multi_sysytem_results/$',multi_sysytem_results),
    url(r'^Alpha ant/$',main),
    url(r'^test_ajax/$',ajax_test),
    url(r'^graph_data/$',graph_data),
    url(r'^input_ajax/$',input_ajax),
    url(r'^compound_name_ajax/',get_compound_infor),
    url(r'^path_information/(?P<c_list>.*?)/',path_information),
    url(r'^reaction_ajax/$',reaction_information),
    url(r'^fba_ajax/$',get_fba),
    url(r'^organ_ajax/$',organ_ajax),
    url(r'^download/(?P<gene>.*?)/',download_gene),
    #url(r'^FDA/(?P<reaction>.*?)/',FDA_result),
    url(r'^admin/', admin.site.urls,name='admin'),
]
