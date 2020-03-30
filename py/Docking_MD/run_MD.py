import MD



def computeMD(prLigand, prMD, pprotein, pranalysis, nameLig, BSCutoff, timeMD, timeframe, stepWait, stepFrame, water, nbCPU, nbGPU):
    """
    :param prLigand: folder of ligands or poses
    :param pprotein: protein in PDB already prepared using the maestro
    :param timeMD: time MD in ps
    :param timeframe: time between two frame
    :param stepWait: time to wait for multiprocess
    :param stepFrame: frame extraction
    :param water: extract water molecules
    :param nbCPU: CPU max in parrallel
    :param nbGPU: GPU max in paralele, random selection
    :return:
    """
    # 1. Merge poses and proteins
    cMDs = MD.MD(prMD, pranalysis, water, timeMD, timeframe, stepWait, nbGPU, nbCPU, stepFrame)
    cMDs.initialisation(prLigand, pprotein)
    cMDs.runMultipleMD()# run MD

    # 2. Extract frames
    # extract frame
    cMDs.centerFrame()
    cMDs.extractFrame()

    # extract BS and ligand
    cMDs.extractLigBSbyFrame(BSCutoff, nameLig, clean=0)

    # for complete run have to control all job done
    #sleep(100)
    # 3. compute RMSD
    cMDs.analyseRMSD()