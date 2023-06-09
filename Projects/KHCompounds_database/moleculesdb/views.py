from django.shortcuts import render, redirect
from django.views import generic
from django.http import HttpResponseRedirect
from django.http import HttpResponse

import csv
from django.contrib import messages

# Create your views here.
#from .models import Compound, CompoundLibraryName, CompoundOrigin, PlateId, CompoundType
from .models import Compound
from .forms import CompoundForm, SearchForm
from django.db.models import Q 

from moleculesdb.search import smarts_search, fast_fp_search
#from .forms import SearchForm, SubmitSingle, UploadFileForm
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from moleculesdb.pagination import get_page_wo_page
#import pybel
from django.contrib.auth.decorators import login_required
from openbabel import pybel
from rdkit.Chem import Draw
from rdkit.Chem.Draw import DrawingOptions




def draw_substructure_matches(mol, smarts):
    # Convert the molecule and query to RDKit format
    mol = mol.to_rdkit()
    query = pybel.Smarts(smarts).to_rdkit()
    
    # Find the atom indices that match the query
    atom_indices = mol.GetSubstructMatch(query)
    
    # Set the drawing options to highlight the matching atoms
    options = DrawingOptions()
    options.atomHighlights = {i:(0, 1, 0) for i in atom_indices}
    
    # Draw the molecule with highlights and return the image as a PNG byte string
    img = Draw.MolToImage(mol, options=options)
    img_str = Draw._getCanvasImage(img).GetData()
    return img_str


def search(request):
    mollist = Compound.objects.all()
    #--------------------------#
    #This takes care of taking the input from the user using the form
    error = []
    form = SearchForm(request.GET)
    form.is_valid()
    # Basic form validation as defined in forms.py
    if form.is_valid() == False:
        for item in form.errors:
            error.append(form[item].errors)
        return render(request, 'search.html', {'error': error, 'form':form})
    moltext = form.cleaned_data['moltext']
    smilesstring = str(form.cleaned_data['smiles'])
    
    # Handling of mol or smi input
    if moltext:
        moltext = pybel.readstring("mol", str(moltext))
        smiles = moltext.write("smi").split("\t")[0]
        print(smiles)
    elif smilesstring:
        smiles = str(smilesstring)
        print(smiles)
    #-----------------------------------#

    #If user makes a similarity of substructure request
    if 'similarity' in request.GET or 'substructure' in request.GET:
        try:
            test = pybel.readstring("smi", smiles)    #Get the smiles
            tanimoto = float(form.cleaned_data['tanimoto'])    #Get the tanimoto value
        except:
            error.append("Molecule or tanimoto index not valid!") # Else send this error message

    if error:
        return render(request, 'search.html', {'error': error, 'form':form})

    if 'similarity' in request.GET and smiles:
        #similartiy search
        mollist = fast_fp_search(mollist, smiles, tanimoto)
        #mollist = smiles_search(mollist, smiles, tanimoto)

        #nofmols = len(mollist)
       
        return render(request, 'moleculesdb/similarity.html', {'compounds':mollist})
        

    elif 'substructure' in request.GET and smiles:
        #substructure search
        #smarts = smiles
        #mollist = smarts_search(mollist, smarts)
        smarts =pybel.Smarts(smiles)
        matches = []

        for compound in mollist:
            try:
                mol = pybel.readstring('smi', compound.smiles)
                if smarts.findall(mol):
                    matches.append(compound)

            except:
                pass

        
        #nofmols = len(mollist)
        return render(request, 'moleculesdb/substructure.html', {'compounds':matches})

    else:
        form = SearchForm()
    return render(request, 'moleculesdb/search.html', {'error': error, 'form':form})




def edit_compound(request, pk):
    compound = Compound.objects.get(id=pk)
    if request.method == 'POST':
        form = CompoundForm(request.POST, request.FILES, instance=compound)
        if form.is_valid():
            new_compound = form.save(commit=False)

            if 'compound_structure' in request.FILES:
                new_compound.compound_structure = request.FILES['compound_structure']


            new_compound.save()

            return redirect('manage-compounds')
    else:
        form = CompoundForm(instance=compound)
    return render(request, 'moleculesdb/edit_compound.html', {'form': form, 'compound': compound})


def delete_compound(request, pk):
    compound = Compound.objects.get(id=pk)
    if request.method == 'POST':
        compound.delete()
        return redirect('manage-compounds')
    return render(request, 'moleculesdb/delete_compound.html', {'compound': compound})


def delete_compounds_list(request):
    if request.method == 'POST':
        selected_compounds = request.POST.getlist('selected_compounds')
        if selected_compounds:
            compounds = Compound.objects.filter(id__in=selected_compounds)
            compounds.delete()
            messages.success(request, 'Selected compounds have been deleted.')

        else:
            messages.warning(request, 'No compounds were selected for deletion.')
            
        return redirect('manage-compounds')

    return render(request, 'moleculesdb/delete_compound_list.html', {'compounds': compounds})


@login_required(login_url='user-login')
def manage_compounds(request):
    compounds = Compound.objects.all()
    return render(request, 'moleculesdb/manage.html', {'compounds': compounds})





def search_results(request):
	if 'query' in request.GET:
		query = request.GET['query']
		#data = Compound.objects.filter(compound_id__icontains=query)
		multiple_query = Q(Q(compound_id__icontains=query) | Q(local_name__icontains=query) | Q(smiles__icontains=query) |
		 Q(compound_library_name__icontains=query) | Q(compound_origin__icontains=query) )
		compounds = Compound.objects.filter(multiple_query)

	else:
		compounds = Compound.objects.all()

	context = {
		'compounds' : compounds
	}

	return render(request,'moleculesdb/search_results.html', context)




def download_csv(request):
    if request.method == 'POST':
        selected_compounds = request.POST.getlist('selected_compounds')
        compounds = Compound.objects.filter(id__in=selected_compounds)
        response = HttpResponse(content_type='text/csv')
        response['Content-Disposition'] = 'attachment; filename="selected_compounds.csv"'
        writer = csv.writer(response)
        writer.writerow(['Compound ID', 'Local name', 'Cas #', 'Smiles', 'Date added', 'compound library name'])
        for compound in compounds:
            writer.writerow([compound.compound_id, compound.local_name, compound.cas_number, compound.smiles, compound.date_added,compound.compound_library_name])
        return response
    else:
        return render(request, 'search_results.html')


def download_csv_list(request):
    if request.method == 'POST':
        selected_compounds = request.POST.getlist('selected_compounds')
        compounds = Compound.objects.filter(id__in=selected_compounds)
        response = HttpResponse(content_type='text/csv')
        response['Content-Disposition'] = 'attachment; filename="selected_compounds.csv"'
        writer = csv.writer(response)
        writer.writerow(['Compound ID', 'Local name', 'Cas #', 'Smiles', 'Date added', 'compound library name'])
        for compound in compounds:
            writer.writerow([compound.compound_id, compound.local_name, compound.cas_number, compound.smiles, compound.date_added,compound.compound_library_name])
        return response
    else:
        return render(request, 'compound_list.html')



def index(request):
	"""View function for home page of site."""

	# Generate counts of some of the main objects
	num_compounds = Compound.objects.all().count()


	context = {
		'num_compounds': num_compounds,
	}

	# Render the HTML template index.html with the data in the context variable
	return render(request, 'index.html', context=context)
    


#@login_required(login_url='labmembers:login_user')

def add_compound(request):
    
    if request.method == "POST":
        form = CompoundForm(request.POST, request.FILES)
        if form.is_valid():
            form.save()

            compound_name = form.cleaned_data['compound_id']
            messages.success(request,f'{compound_name} has been added')
            return redirect('compounds')


    else:
        form = CompoundForm()

    return render(request,'moleculesdb/add_compound.html',{'form':form})





# Generate CSV File Venue List
@login_required(login_url='user-login')
def compound_csv(request):
	response = HttpResponse(content_type='text/csv')
	response['Content-Disposition'] = 'attachment; filename=compound_list.csv'
	
	# Create a csv writer
	writer = csv.writer(response)

	# Designate The Model
	#compounds = Compound.objects.all()

	# Add column headings to the csv file
	writer.writerow(['Compound ID', 'Local name', 'Cas #', 'Smiles', 'Date added', 'compound library name'])

	# Loop Thu and output
	#for compound in compounds:
	#	writer.writerow([Compound.compound_id, Compound.local_name, Compound.cas_number, Compound.smiles, Compound.date_added])

	for compound in Compound.objects.all().values_list('compound_id', 'local_name','cas_number','smiles','date_added','compound_library_name'):
		writer.writerow(compound)

	return response








class CompoundListView(generic.ListView):
    model = Compound

    #paginate_by = 10




class CompoundDetailView(generic.DetailView):
    model = Compound





