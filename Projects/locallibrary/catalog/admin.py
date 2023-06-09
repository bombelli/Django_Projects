from django.contrib import admin
from .models import Author, Book, Genre, BookInstance, Language

#admin.site.register(Book)
#admin.site.register(Author)
admin.site.register(Genre)
#admin.site.register(BookInstance)


class BookInstanceInline(admin.TabularInline):
	model = BookInstance

	extra = 0

class BookInline(admin.TabularInline):
	model = Book

	extra = 0 

class BookAdmin(admin.ModelAdmin):
	list_display = ('title', 'author', 'display_genre') # Display genre is a function because Genre is a foreign key

	inlines = [BookInstanceInline]

admin.site.register(Book, BookAdmin)



class AuthorAdmin(admin.ModelAdmin):
	list_display = ('last_name', 'first_name', 'date_of_birth', 'date_of_death') # The fields/colums you want to display in your admin page 

	fields = ['first_name','last_name', ('date_of_birth','date_of_death')]  #This will be displayed on the form in order as listed, the tuple is horizontal display

	#exclude = ['last_name']
	inlines = [BookInline]

admin.site.register(Author,AuthorAdmin)




class BookInstanceAdmin(admin.ModelAdmin):

	list_display = ('book','status','due_back','id')

	list_filter = ('status', 'due_back')  # A filter session will be displayed

	fieldsets = (

		(None, {'fields': ('book','imprint','id')}),
		('Availability', {'fields':('status','due_back')}),
		)       #Adding sections, each section has a title # All these is done on the forms
admin.site.register(BookInstance,BookInstanceAdmin)

# Register your models here.
