    try:
        mol = pybel.readstring("smi",smiles) # Read the smiles
        descs = mol.calcdesc()  # Calculate descriptors of smiles
        
        #generate 2D coordinates, needs openbabel
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("smi", "mdl")
        obmol = openbabel.OBMol()
        obConversion.ReadString(obmol, smiles)
        gen2d = openbabel.OBOp.FindType("gen2d")
        gen2d.Do(obmol)
        MDL = obConversion.WriteString(obmol)
        outMDL = MDL.replace("\n", r"\n")
        CMW=descs["MW"]
        CHN=mol.formula
        HBA=descs["HBA1"]
        HBD=descs["HBD"]
        logP=descs["logP"]
        tpsa=descs["TPSA"]
        #Get number  of rotatable bonds
        smarts = pybel.Smarts(r"[!$([NH]!@C(=O))&!D1&!$(*#*)]\&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]")
        rb = smarts.findall(mol)
        nrb = len(rb)
        #Calculate Fsp3
        sp3c = pybel.Smarts("[CX4]")
        nsp3c = sp3c.findall(mol)
        nsp3c = float(len(nsp3c))
        allc =  pybel.Smarts("[#6]")
        nallc = allc.findall(mol)
        nallc = float(len(nallc))
        if nallc > 0:
            fsp3 = nsp3c/nallc
        else:
            fsp3 = ""
        #Get fingerprint and molecular complexity
        fprint = mol.calcfp()
        bitson = fprint.bits
        nbitson = len(bitson)
        if 'hts' in molclass.lower() or 'compound' in molclass.lower():
            #print "hts"
            pains = detect_pains(mol)
        else:
            pains = 'Not checked'
        m = Molecule(name=name,SMILES=smiles, altname=altname, supplier=supplier, supplierID=supplierID, 
        CMW=descs["MW"], CHN=CHN, HBA=HBA, HBD=HBD, logP=logP, tpsa=tpsa, amount=amount, unit=unit, 
        CAS=cas, storage=storage, storageID=storageID, molfile=outMDL, nrb=nrb, fingerprint=bitson, complexity=nbitson, 
        comment=comment, molclass=molclass, fsp3=fsp3, pains=pains, platebarcode=platebarcode, 
        samplebarcode=samplebarcode, randomstring=randomstring)
        m.save()
    except:
        # OpenBabel failed, no properties, etc..
        m = Molecule(name=name,SMILES=smiles, altname=altname, supplier=supplier, supplierID=supplierID, 
        amount=amount, unit=unit, CAS=cas, storage=storage, storageID=storageID, comment=comment, molclass=molclass, 
        platebarcode=platebarcode, samplebarcode=samplebarcode, randomstring=randomstring)
        m.save()
        #Save data to database
