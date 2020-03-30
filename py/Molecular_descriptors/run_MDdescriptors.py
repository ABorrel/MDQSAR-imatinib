def computeMDdesc(prMD, prout, istart=0, iend=0, descLig=1, descBS=1, descFPI=1 ):

    lpMDresult = listdir(prMD)
    if iend == 0:
        iend = len(lpMDresult)
    print prMD

    i = istart
    for MDresult in lpMDresult[istart:iend]:
        jobname = MDresult.split("_")[0]
        prlig = prMD + MDresult + "/lig/"
        prBS = prMD + MDresult + "/BSs/"
        prframes = prMD + MDresult + "/framesMD/"
        prMDdesc = pathFolder.createFolder(prout + jobname + "/")


        # control run
        print jobname, i
        i +=1

        # compute different descriptors
        if not "cMD" in locals().keys():
            cMD = MDdescriptors.MDdescriptors(jobname, prlig, prBS, prframes, prMDdesc)
        else:
            cMD.jobname = jobname
            cMD.prlig = prlig
            cMD.prBSs = prBS
            cMD.prframe = prframes
            cMD.prout = prMDdesc

        if descLig == 1:
            cMD.computeLigDesc()
        if descBS == 1:
            cMD.computeBSDesc()
        if descFPI == 1:
            cMD.computeFPI()

