from django import forms
from django.forms import ModelForm
#from .models import Compound, CompoundLibraryName, CompoundOrigin, PlateId, CompoundType
from .models import Compound


class CompoundForm(ModelForm):


	fieldsets = [
		('Compound Basic Annotation', {'fields': ['date_added', ('compound_id','compound_structure'),('compound_type','compound_library_name','compound_origin'), ('local_name','commercial_name'),('cas_number','product_code'),('au_contact_person_name','au_contact_person_email'),('plate_id','solubility_and_solvent_to_use'),]
		}),
		('Compound Properties', {'fields': [('molecular_formula','molecular_weight','logp'),('smiles','compound_salt_version','inchl_key'),('chembl_id','pubchem_cid'),('iupac_name','commercially_available'),]
		}),
		('Storage Annotation', {'fields': ['storage_type',('quantity_received', 'weight','volume'),('storage_format','room_temperature'),('freezer_name','freezer_temperature'),('shelf_or_drawer_number','box_number'),]
		}),
		('Activity Annotation', {'fields': [('assay_date','previous_assay_date'),('affects_schistosomula', 'schistosomula_activity'),('affects_juvenile_worms',
			'juvenile_worms_activity'),('affects_adult_worms','adult_worms_activity','adultworm_screen_conducted_by'),]
		}),
		('Biological Targets Annotation', {'fields': ['putative_human_target', 'mechanism_of_action','invivo_data_available',]
		}),
		('Comments', {'fields':['comments',]}),
		]

		
	class Meta:
		model = Compound
		fields = "__all__"
		exclude = ['fingerprint']

		widgets = {
            'date_added': forms.DateInput(attrs={'type': 'date'}),
            'solubility_and_solvent_to_use': forms.Textarea(attrs={'rows': 4}),
            'mechanism_of_action':forms.Textarea(attrs={'rows': 4, 'columns': 8}),
            'comments': forms.Textarea(attrs={'class': 'form-control','rows': 5}),

        }
	
	
		
 

class SearchForm(forms.Form):
    moltext = forms.CharField(required=False, widget=forms.Textarea(attrs={'class':'large-4', 'placeholder':'Mol data from editor', 'rows':'1'}))
    smiles = forms.CharField(required=False, widget=forms.TextInput(attrs={'placeholder':'SMILES'}))
    tanimoto = forms.FloatField(required=False, max_value=1.0, min_value=0.0, initial=0.7, error_messages={'max_value': 'Tanimoto: enter a number between 0 and 1!', 'min_value': 'Tanimoto: enter a number between 0 and 1!', 'invalid': 'Tanimoto: enter a number between 0 and 1!'}, widget=forms.TextInput(attrs={'class':'large, columns', 'value':'0.7'}))
    
    