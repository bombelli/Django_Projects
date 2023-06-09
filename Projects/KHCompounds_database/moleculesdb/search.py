from moleculesdb.models import Compound
#import pybel
#from sets import Set
from openbabel import pybel

try:
    from sets import Set
except ImportError:
    Set = set



def fast_fp_search(mollist, smiles, tanimoto):
    #Should be faster, since it uses precalculated fingerprints
    ret = []
    tanimoto = float(tanimoto)
    query = pybel.readstring("smi", smiles)
    fpquery = query.calcfp()
    fpquery = fpquery.bits
    for mol in mollist.all():
        try:
            fp = eval(str(mol.fingerprint))
            a = Set(fpquery)
            b = Set(fp)
            un = float(len(a.union(b)))
            inx = float(len(a.intersection(b)))
            tan = inx/un
            #print tan
            if tan > tanimoto:
                ret.append(mol)
        except:
            pass
    return ret

def smarts_search(mollist, smarts):
    ret = []
    query = pybel.Smarts(smarts)
    #print(query)
    for mol in mollist.all():
        try:
            smiles = pybel.readstring("smi", str(mol.SMILES))
            if query.findall(smiles):
                ret.append(mol)
        except:
            pass
    return ret
