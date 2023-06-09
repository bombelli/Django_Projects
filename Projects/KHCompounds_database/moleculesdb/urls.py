from django.urls import path
from . import views
from django.contrib.auth.decorators import login_required



urlpatterns = [
    path('', views.index, name='index'),
    path('compounds/', views.CompoundListView.as_view(), name='compounds'),
    path('compound/<int:pk>', views.CompoundDetailView.as_view(), name='compound-detail'),
    path('add_compound/', views.add_compound, name='add-compound'),
    path('search_results/', views.search_results, name='search-results'),
    path('compound_csv', views.compound_csv, name='compound-csv'),
    path('search/', views.search, name='search-new'),
    #path('select_compounds/', views.search, name='select-compounds'),
    path('download_csv/', views.download_csv, name='download-csv'),
    path('download_csv_list/', views.download_csv_list, name='download-csv-list'),
    path('manage/', views.manage_compounds, name='manage-compounds'),
    path('manage/edit/<int:pk>/', views.edit_compound,name='edit-compound'),
    path('manage/delete/<int:pk>/', views.delete_compound,name='delete-compound'),
    path('manage/delete_compounds/', views.delete_compounds_list,name='delete-compounds-list'),

   
]


#path('book/create/', login_required(book_create, login_url='labmembers:login'), name='book_create'),
