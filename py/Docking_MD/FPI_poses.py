def FPIMatrix(cdocking, pprotein, prFPI):

    pmatrixFPI = prFPI + "MFPI.txt"

    i = 0
    imax = len(cdocking.lposefiles)
    #imax = 3
    while i < imax:
        j = i + 1
        while j < imax:
            cprot = PDB.PDB(PDB_input=pprotein, hydrogen=1)
            pligPDBi = runExternalSoft.babelConvertSDFtoPDB(cdocking.lposefiles[i])
            pligPDBj = runExternalSoft.babelConvertSDFtoPDB(cdocking.lposefiles[j])

            print pligPDBi
            cposei = PDB.PDB(PDB_input=pligPDBi, hydrogen=1)
            cposej = PDB.PDB(PDB_input=pligPDBj, hydrogen=1)

            sFPIi = FPI.ligFPI(cPDB=cprot, ligin=cposei, prFPI=prFPI)
            sFPIj = FPI.ligFPI(cPDB=cprot, ligin=cposej, prFPI=prFPI)

            dout = sFPIi.compareFPI(sFPIj)

            print dout

            j = j + 1
        i = i + 1
