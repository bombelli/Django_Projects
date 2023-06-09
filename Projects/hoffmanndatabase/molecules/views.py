from django.shortcuts import render, redirect
from .models import Molecule
from .forms import MoleculeForm, SearchForm
from django.contrib import messages
from django.db.models import Q 
from django.views import generic
from django.http import HttpResponseRedirect
from django.http import HttpResponse
from django.contrib.auth.decorators import login_required
from molecules.search import smarts_search, fast_fp_search
from openbabel import pybel
import csv




def index(request):

	num_compounds = Molecule.objects.all().count()

	context = {
		'num_compounds': num_compounds,
	}

	return render(request, 'index.html', context=context)

#-----------------------------------------------------------------------------
@login_required
def add_compound(request):
    
    if request.method == "POST":
        form = MoleculeForm(request.POST, request.FILES)
        if form.is_valid():
            form.save()

            compound_name = form.cleaned_data['compound_id']
            messages.success(request,f'{compound_name} has been added')
            
    else:
        form = MoleculeForm()

    return render(request,'molecules/add_compound.html',{'form':form})



#-----------------------------------------------------------------------------------

def search_results(request):

	if 'query' in request.GET:

		query = request.GET['query']

		multiple_query = Q(Q(compound_id__icontains=query) | Q(local_name__icontains=query) | Q(smiles__icontains=query) |
		 Q(compound_library_name__icontains=query) | Q(compound_origin__icontains=query) )

		compounds = Molecule.objects.filter(multiple_query)

	else:
		compounds = Molecule.objects.all()

	context = {
		'compounds' : compounds
	}

	return render(request,'molecules/search_results.html', context)


#---------------------------------------------------------------------------------------------------------

class MoleculeListView(generic.ListView):
    model = Molecule

    #paginate_by = 10

#----------------------------------------------------------------------------------------------------------


class MoleculeDetailView(generic.DetailView):
    model = Molecule


#--------------------------------------------------------------------------------------------------------

@login_required
def compound_csv(request):
	response = HttpResponse(content_type='text/csv')
	response['Content-Disposition'] = 'attachment; filename=compound_list.csv'
	
	# Create a csv writer
	writer = csv.writer(response)

	# Add column headings to the csv file
	writer.writerow(['Compound ID', 'Local name', 'Cas #', 'Smiles', 'Date added', 'compound library name'])

	for compound in Molecule.objects.all().values_list('compound_id', 'local_name','cas_number','smiles','date_added','compound_library_name'):
		writer.writerow(compound)

	return response

#---------------------------------------------------------------------------------------------------------------------------------

def download_csv(request):
    if request.method == 'POST':
        selected_compounds = request.POST.getlist('selected_compounds')
        compounds = Molecule.objects.filter(id__in=selected_compounds)
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
        compounds = Molecule.objects.filter(id__in=selected_compounds)
        response = HttpResponse(content_type='text/csv')
        response['Content-Disposition'] = 'attachment; filename="selected_compounds.csv"'
        writer = csv.writer(response)
        writer.writerow(['Compound ID', 'Local name', 'Cas #', 'Smiles', 'Date added', 'compound library name'])
        for compound in compounds:
            writer.writerow([compound.compound_id, compound.local_name, compound.cas_number, compound.smiles, compound.date_added,compound.compound_library_name])
        return response
    else:
        return render(request, 'compound_list.html')


#-------------------------------------------------------------------------------------------------------------------------------------------

@login_required
def manage_compounds(request):
    compounds = Molecule.objects.all()
    return render(request, 'molecules/manage.html', {'compounds': compounds})

def edit_compound(request, pk):

    compound = Molecule.objects.get(id=pk)

    if request.method == 'POST':
        form = MoleculeForm(request.POST, request.FILES, instance=compound)
        if form.is_valid():
            new_compound = form.save(commit=False)

            if 'compound_structure' in request.FILES:
                new_compound.compound_structure = request.FILES['compound_structure']


            new_compound.save()

            return redirect('manage-compounds')
    else:

        form = MoleculeForm(instance=compound)

    return render(request, 'molecules/edit_compound.html', {'form': form, 'compound': compound})


def delete_compound(request, pk):
    compound = Molecule.objects.get(id=pk)
    if request.method == 'POST':
        compound.delete()
        return redirect('manage-compounds')
    return render(request, 'molecules/delete_compound.html', {'compound': compound})


def delete_compounds_list(request):
    if request.method == 'POST':
        selected_compounds = request.POST.getlist('selected_compounds')
        if selected_compounds:
            compounds = Molecule.objects.filter(id__in=selected_compounds)
            compounds.delete()
            messages.success(request, 'Selected compounds have been deleted.')

        else:
            messages.warning(request, 'No compounds were selected for deletion.')
            
        return redirect('manage-compounds')

    return render(request, 'molecules/delete_compound_list.html', {'compounds': compounds})

def manage_download_csv(request):
    if request.method == 'POST':
        selected_compounds = request.POST.getlist('selected_compounds')
        compounds = Molecule.objects.filter(id__in=selected_compounds)

        response = HttpResponse(content_type='text/csv')
        response['Content-Disposition'] = 'attachment; filename="selected_compounds.csv"'
        writer = csv.writer(response)
        writer.writerow(['Compound ID', 'Local name', 'Cas #', 'Smiles', 'Date added', 'compound library name'])

        for compound in compounds:
            writer.writerow([compound.compound_id, compound.local_name, compound.cas_number, compound.smiles, compound.date_added,compound.compound_library_name])

        return response

    return HttpResponse("Invalid request")


#-------------------------------------------------------------------------------------------------------------

@login_required
def search(request):
    mollist = Molecule.objects.all()
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
       
        return render(request, 'molecules/similarity.html', {'compounds':mollist})
        

    elif 'substructure' in request.GET and smiles:
        #substructure search
        smarts =pybel.Smarts(smiles)
        matches = []

        for compound in mollist:
            try:
                mol = pybel.readstring('smi', compound.smiles)
                if smarts.findall(mol):
                    matches.append(compound)

            except:
                pass

        return render(request, 'molecules/substructure.html', {'compounds':matches})

    else:
        form = SearchForm()
    return render(request, 'molecules/search.html', {'error': error, 'form':form})


#------------------------------------------------------------------------------------------------------------