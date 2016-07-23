import moldb.views

"""molinf URL Configuration

The `urlpatterns` list_molecules routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.9/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.conf.urls import url, include
    2. Add a URL to urlpatterns:  url(r'^blog/', include('blog.urls'))
"""
from django.conf.urls import url
from django.contrib import admin
from django.contrib.auth import login as django_login

urlpatterns = [
    # pages
    #url(r'^admin/', admin.site.urls),
    url(r'^$', moldb.views.home, name='home'),
    url(r'^add-molecules$', moldb.views.add_molecules, name='add_molecules'),
    url(r'^list-molecules$', moldb.views.list_molecules, name='list_molecules'),
    url(r'^search-molecules$', moldb.views.search_molecules, name='search_molecules'),
    url(r'^purchases$', moldb.views.purchases, name='purchases'),
    url(r'^orders$', moldb.views.orders, name='orders'),
    url(r'^login$', moldb.views.login_user, name='login_user'),
    url(r'^logout$', moldb.views.logout_user, name='logout_user'),

    url(r'^test$', moldb.views.test, name='test'),

    # API
    url(r'^api/molConverter$', moldb.views.api_molConverter, name='api_molConverter'),
    url(r'^api/addMolecule$', moldb.views.api_addMolecule, name='api_addMolecule'),
    url(r'^api/uploadMolecules$', moldb.views.api_uploadMolecules, name='api_uploadMolecules'),
    url(r'^api/uploadMolecules/status$', moldb.views.api_uploadMoleculesStatus, name='api_uploadMoleculesStatus'),
    url(r'^api/downloadMolecules$', moldb.views.api_downloadMolecules, name='api_downloadMolecules'),
    url(r'^api/searchMoleculesByStructure$', moldb.views.api_searchMoleculesByStructure, name='api_searchMoleculesByStructure'),
    url(r'^api/searchMoleculesByFilter$', moldb.views.api_searchMoleculesByFilter, name='api_searchMoleculesByFilter'),
    url(r'^api/addOrder$', moldb.views.api_addOrder, name='api_addOrder'),
    url(r'^api/addPurchase$', moldb.views.api_addPurchase, name='api_addPurchase'),
    url(r'^api/rmMolFromInProgressPurchase$', moldb.views.api_rmMolFromInProgressPurchase, name='api_rmMolFromInProgressPurchase'),
    url(r'^api/rmInProgressPurchase$', moldb.views.api_rmInProgressPurchase, name='api_rmInProgressPurchase'),
    url(r'^api/finishInProgressPurchase$', moldb.views.api_finishInProgressPurchase, name='api_finishInProgressPurchase'),
    url(r'^api/rmMolFromInProgressOrder$', moldb.views.api_rmMolFromInProgressOrder, name='api_rmMolFromInProgressOrder'),
    url(r'^api/rmInProgressOrder$', moldb.views.api_rmInProgressOrder, name='api_rmInProgressOrder'),
    url(r'^api/finishInProgressOrder$', moldb.views.api_finishInProgressOrder, name='api_finishInProgressOrder')
]