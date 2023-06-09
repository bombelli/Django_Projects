from django.db import models
from django.urls import reverse
import uuid
from django.contrib.auth.models import User
from datetime import date

from openbabel import pybel



# Create your models here.

class Compound(models.Model):
	YES_OR_NO = (
		('Yes','Yes'),
		('No', 'No'),
		)

	STORAGE_FORMAT = (
		('Ind_vial','Individual vial'),
		('96wp', '96 well plate'),
		('384wp', '384 well plate'),
		)

	STORAGE_TEMPERATURE = (
		('4','4 degrees'),
		('-20', '-20 degrees'),
		('-80', '-80 degrees'),
		)

	DRY_OR_LIQUID = (
		('dry','dry'),
		('liquid', 'liquid'),
		)

	STORAGE_TYPE = (
		('Storage Pod','Storage Pod'),
		('Non Storage Pod', 'Non Storage Pod'),
		)

	#inventory_num = models.IntegerField(blank=True, null=True)
	date_added = models.DateField(null=True, blank=True) #Fixed
	compound_structure = models.ImageField(null=True, blank=True, upload_to = "compound_images/", help_text='Upload image of the compound structure') #Image of the compound structure

	compound_id = models.CharField(max_length=30,unique=True, help_text= 'Enter compound id') # This is a unique id
	local_name = models.CharField(max_length=100,blank=True, null=True, help_text='Enter local name')
	commercial_name = models.CharField(max_length=100,blank=True, null=True, help_text='Enter compound commercial name')
	cas_number = models.CharField(max_length=100,blank=True, null=True, help_text='Enter CAS number')
	product_code = models.CharField(max_length=100,blank=True, null=True, help_text='For only purchased compounds') #fixed

	#compound_library_name = models.ForeignKey('CompoundLibraryName', on_delete= models.SET_NULL,blank=True, null=True)
	#compound_origin = models.ForeignKey('CompoundOrigin',on_delete= models.SET_NULL,blank=True, null=True) #Remove inverted commas on foreign keys

	compound_library_name = models.CharField(max_length=200,blank=True, null=True, verbose_name='Compound Library Name', help_text='Enter the compound library name')
	compound_origin = models.CharField(max_length=200,blank=True, null=True, help_text='Enter the compound origin name')

	au_contact_person_name = models.CharField(max_length = 250,verbose_name ='AU Contact Name', blank=True, null=True,help_text='Enter full name' )
	au_contact_person_email = models.EmailField(max_length = 250,verbose_name ='AU Contact Email', blank=True, null=True,help_text='Enter email address' ) 

	#plate_id = models.ManyToManyField('PlateId',blank=True, help_text='Enter Plate Id')#**
	plate_id = models.CharField(max_length=100,blank=True, null=True, help_text='Enter Plate Id') #New variable because foreign key giving me a lot of issues

	solubility_and_solvent_to_use = models.TextField(max_length=200,null=True, blank=True, help_text='Solvent and solvent to use')



	comments = models.TextField(max_length=1000, null=True, blank=True, help_text='Additional information') # Comments/ additional information
	

	#Make sure you put storage under one session
	quantity_received = models.CharField(max_length =10, choices = DRY_OR_LIQUID, blank = True, help_text = 'Quantity received')
	weight = models.FloatField(null=True, blank=True, help_text= 'Enter weight if dry (in mg)')
	volume = models.FloatField(null=True, blank=True, help_text= 'Enter volume if liquid (in μL)')
	storage_type = models.CharField(max_length =20, choices = STORAGE_TYPE, blank = True, help_text = 'Type of storage used')

	room_temperature = models.FloatField(null=True, blank=True, verbose_name = 'Room Temperature')
	freezer_temperature = models.CharField(max_length =4, choices = STORAGE_TEMPERATURE, blank = True, help_text = 'Freezer Temperature')
	storage_format = models.CharField(max_length =50, choices = STORAGE_FORMAT, blank = True, help_text = 'Storage format')
	freezer_name = models.CharField(max_length=50, blank=True, null=True, help_text='Freezer name')
	shelf_or_drawer_number = models.IntegerField(blank=True, null=True)
	box_number = models.IntegerField(blank=True, null=True)

	#compound_type = models.ForeignKey('CompoundType',blank=True, on_delete= models.SET_NULL, null=True)

	compound_type = models.CharField(max_length=200,blank=True, null=True, help_text='Enter compound type (eg. Small molecule) ')

	molecular_formula = models.CharField(max_length=30,null=True, blank=True,help_text= 'Enter compound formula')
	molecular_weight = models.FloatField(null=True, blank=True, verbose_name = 'Compound molecular weight')
	smiles = models.CharField(max_length=500,null=True, blank=True, help_text= 'Enter compound Smiles')
	inchl_key = models.CharField(max_length=200,null=True, blank=True, help_text= 'Enter compound Inchl key')
	logp = models.FloatField(null=True, blank=True, help_text= 'Enter LogP/LogD ')

	#name annotation
	pubchem_cid = models.IntegerField(blank=True, null=True)
	chembl_id = models.CharField(max_length=30,null=True, blank=True,help_text= 'Enter compound id')
	iupac_name = models.CharField(max_length=500,blank=True, null=True, help_text= 'Enter IUPAC name') #Fixed
	commercially_available = models.CharField(max_length =4, choices = YES_OR_NO, blank = True, help_text = 'Commercial Availability')
	compound_salt_version = models.CharField(max_length=500,blank=True, null=True,help_text= 'Enter compound salt version') #Fixed

	#Activity data
	affects_schistosomula = models.CharField(max_length =4, choices = YES_OR_NO, blank = True, help_text = 'Does compound affect Schistosomula')
	schistosomula_activity = models.FloatField(null=True, blank=True, verbose_name = 'Schistosomula Activity(EC50)', help_text='(in μM)')
	affects_juvenile_worms = models.CharField(max_length =4, choices = YES_OR_NO, blank = True, help_text = 'Does compound affect Juvenile worms')
	juvenile_worms_activity = models.FloatField(null=True, blank=True, verbose_name = 'Juvenile worms Activity(EC50)', help_text='(in μM)')
	affects_adult_worms = models.CharField(max_length =4, choices = YES_OR_NO, blank = True, help_text = 'Does compound affect Adultworms')
	adult_worms_activity = models.FloatField(null=True, blank=True, verbose_name = 'Adult worms Activity(EC50)', help_text='(in μM)')
	adultworm_screen_conducted_by= models.CharField(max_length=100,blank=True, null=True, help_text='Enter name of person who conducted the adultworm screening')
	assay_date = models.DateField(null=True, blank=True)
	previous_assay_date = models.DateField(null=True, blank=True)

	#Biological targets
	putative_human_target = models.CharField(max_length=100,blank=True, null=True, help_text='Enter its putative human target')
	mechanism_of_action = models.TextField(max_length=1000,blank=True, null=True, help_text='Compound mechanism of action')
	invivo_data_available = models.CharField(max_length =4, choices = YES_OR_NO, blank = True, help_text = 'Is there any invivo data on the compound')
	fingerprint = models.TextField(blank=True, verbose_name='Fingerprint')



	def save(self, *args, **kwargs):
		if not self.smiles: # Check if the 'smiles' field is empty
			return
		mol = pybel.readstring("smi",self.smiles)
		fprint = mol.calcfp()
		bitson = fprint.bits
		self.fingerprint = bitson
		super().save(*args, **kwargs)

	
	def __str__(self):

		return self.compound_id 


	class Meta:
		ordering = ['-date_added','compound_id','local_name','compound_library_name','compound_origin']


	def get_absolute_url(self):
		return reverse('compound-detail', args=[str(self.id)])



'''
class CompoundLibraryName(models.Model):

	library_name = models.CharField(max_length=200,verbose_name='Compound Library Name', help_text='Enter the compound library name')


	def __str__(self):

		return self.library_name


class CompoundOrigin(models.Model):

	compound_origin_name = models.CharField(max_length=200, help_text='Enter the compound origin name')


	def __str__(self):

		return self.compound_origin_name


class PlateId(models.Model):

	plate_id = models.CharField(max_length=200,verbose_name='Plate ID', help_text='Enter Plate ID')


	def __str__(self):

		return self.plate_id

class CompoundType(models.Model):

	compound_type = models.CharField(max_length=200, help_text='Enter compound type (eg. Small molecule) ')


	def __str__(self):
  
		return self.compound_type


'''



