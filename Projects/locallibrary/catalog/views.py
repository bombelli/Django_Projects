from django.shortcuts import render
from .models import Book, Author, BookInstance, Genre
from django.views import generic


'''
def index(request):
	num_books = Book.objects.all().count()

	num_instances = BookInstance.objects.all().count()

	num_instances_available = BookInstance.objects.filter(status__exact='a').count()

	genre = Genre.objects.all().count()

	# The 'all()' is implied by default.
	num_authors = Author.objects.count()


	context = {
		'num_books': num_books,
		'num_instances': num_instances,
		'num_instances_available': num_instances_available,
		'num_authors' : num_authors,
		'genre': genre
	}


	return render(request,'index.html',context=context)
'''
def index(request):
    # …
	num_books = Book.objects.all().count()

	num_instances = BookInstance.objects.all().count()

	num_instances_available = BookInstance.objects.filter(status__exact='a').count()

	genre = Genre.objects.all().count()

	num_authors = Author.objects.count()  # The 'all()' is implied by default.

	# Number of visits to this view, as counted in the session variable.
	num_visits = request.session.get('num_visits', 0)
	request.session['num_visits'] = num_visits + 1

	context = {
		'num_books': num_books,
		'num_instances': num_instances,
		'num_instances_available': num_instances_available,
		'num_authors': num_authors,
		'num_visits': num_visits,
	}

	# Render the HTML template index.html with the data in the context variable.
	return render(request, 'index.html', context=context)



class BookListView(generic.ListView):
	model = Book
	paginate_by = 2

# Create your views here.


class BookDetailView(generic.DetailView):
    model = Book


class AuthorListView(generic.ListView):
	model = Author
	paginate_by= 2


class AuthorDetailView(generic.DetailView):
	model = Author
