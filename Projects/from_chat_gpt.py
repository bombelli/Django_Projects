{% extends 'base.html' %}

{% block content %}
  <h1>Search Results</h1>
  <form method="post" action="{% url 'select_books' %}">
    {% csrf_token %}
    {% if books %}
      <p>Select the books you want to add:</p>
      <ul>
        {% for book in books %}
          <li>
            <input type="checkbox" name="books" value="{{ book.id }}">
            <a href="{% url 'book_detail' book.id %}">{{ book.title }}</a> by {{ book.author }}
          </li>
        {% endfor %}
      </ul>
      <input type="submit" value="Add Selected Books">
      <br><br>
      <p>Or, download a CSV file with the search results:</p>
      <a href="{% url 'download_csv' %}?q={{ q }}">Download CSV</a>
      <p>Preview:</p>
      <table>
        <tr>
          <th>Title</th>
          <th>Author</th>
          <th>Publisher</th>
          <th>Publication Date</th>
          <th>ISBN</th>
        </tr>
        {% for book in books %}
          <tr>
            <td>{{ book.title }}</td>
            <td>{{ book.author }}</td>
            <td>{{ book.publisher }}</td>
            <td>{{ book.publication_date }}</td>
            <td>{{ book.isbn }}</td>
          </tr>
        {% endfor %}
      </table>
    {% else %}
      <p>No books found.</p>
    {% endif %}
  </form>
{% endblock %}



'''
---------------------------------------------------------------------------------------------------------------------------------
This template is similar to the previous one, with the addition of a preview table that shows how the CSV file will look like. 
The table has columns for the book's title, author, publisher, publication date, and ISBN.

When the user selects one or more books and clicks the "Add Selected Books" button, the selected books are added to their collection.
 When the user clicks the "Download CSV" link, a CSV file containing the search results is generated and downloaded. 
The preview table shows the same information as the CSV file, to give the user an idea of how it will look like.

--------------------------------------------------------------------------------------------------------------------------
This template is similar to the previous one, with the addition of a link to download the search results as a CSV file. 
It assumes that you have a view called download_csv that generates the CSV file, and that you pass the search query (q) as a parameter.

When the user clicks the "Download CSV" link, the current search query is passed as a parameter in the URL. 
The download_csv view can use this query to generate a CSV file containing the search results.

Note that the CSV download link is outside the form, since it doesn't need to submit any form data. 
It's also a good practice to include a line break (<br><br>) between the form and the download link to make the layout more clear.
---------------------------------------------------------------------------------------------------------------------------------------
'''


# View to download 
import csv
from django.http import HttpResponse

def download_csv(request):
    # Get the selected books from request parameters
    selected_books = ...

    # Define the CSV headers and rows
    headers = ['Title', 'Author', 'Publisher', 'Publication Year']
    rows = [[book.title, book.author, book.publisher, book.publication_year] for book in selected_books]

    # Create a new CSV response object
    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="book_selections.csv"'

    # Write the headers and rows to the CSV response
    writer = csv.writer(response)
    writer.writerow(headers)
    for row in rows:
        writer.writerow(row)

    return response


'''
Note that this is just an example implementation, 
and you'll need to modify it to suit your specific use case. 
In particular, you'll need to replace the selected_books 
variable with the appropriate data source for your application.
'''