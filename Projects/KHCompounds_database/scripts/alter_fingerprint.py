from moleculesdb.models import Compound, CompoundLibraryName, CompoundOrigin, PlateId, CompoundType

def run():
	compounds = Compound.objects.all()

	for compound in compounds:
		compound.save()