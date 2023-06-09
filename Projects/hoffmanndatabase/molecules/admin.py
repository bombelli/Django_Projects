from django.contrib import admin
from import_export.admin import ExportActionMixin

from .models import Molecule

#admin.site.register(Molecule)


class MoleculeAdmin(ExportActionMixin, admin.ModelAdmin):

	exclude = ('fingerprint',)
	
	list_display = ('compound_id','local_name','compound_library_name', 'compound_origin','date_added')

	list_filter = ('compound_library_name', 'compound_origin', 'compound_type')

	search_fields = ['compound_id','compound_library_name', 'compound_origin']


	fieldsets = (
		('Compound Basic Annotation', {'fields': ('date_added', ('compound_id','compound_structure'), ('local_name','commercial_name'),('cas_number','product_code'),('au_contact_person_name','au_contact_person_email'),'compound_type',('compound_library_name','compound_origin'),'plate_id','solubility_and_solvent_to_use')
		}),
		('Compound Properties', {'fields': ('molecular_formula',('smiles','compound_salt_version'),'inchl_key',('molecular_weight','logp'),('chembl_id','pubchem_cid'),'iupac_name','commercially_available')
		}),
		('Storage Annotation', {'fields': ('quantity_received', ('weight','volume'),'storage_type','room_temperature','storage_format',('freezer_name','freezer_temperature'),('shelf_or_drawer_number','box_number'))
		}),
		('Activity Annotation', {'fields': (('assay_date','previous_assay_date'),('affects_schistosomula', 'schistosomula_activity'),('affects_juvenile_worms',
			'juvenile_worms_activity'),('affects_adult_worms','adult_worms_activity','adultworm_screen_conducted_by'))
		}),
		('Biological Targets Annotation', {'fields': ('putative_human_target', 'mechanism_of_action','invivo_data_available')
		}),
		('Comments', {'fields':('comments',)}),
		)

admin.site.register(Molecule, MoleculeAdmin)
