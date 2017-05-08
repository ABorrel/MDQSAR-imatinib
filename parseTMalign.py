from os import path


def parseOutputTMalign(p_filin):
    # print p_filin
    if not path.exists(p_filin):
        return {}
    d_out = {}

    try:
        filin = open(p_filin, "r")
        f_read = filin.read()
        filin.close()

        d_out["RMSD"] = retrieveRMSD(f_read)
        d_out["IDseq"] = retrieveIDseq(f_read)
        l_score = retrieveTMScore(f_read)
        d_out["TMscore1"] = l_score[0]
        d_out["TMscore2"] = l_score[1]

        return d_out
    except:
        return d_out


def retrieveRMSD(filin_read):
    """
    """
    part_file = filin_read.split("RMSD=")[1]
    RMSD = part_file.split(",")[0]
    try:
        RMSD = RMSD.replace(" ", "")
    except:
        pass

    return RMSD


def retrieveIDseq(filin_read):
    """
    Retrieve RMSD in TMalign out file
    args: -> path file RMSD
    return: value of RMSD
    """

    IDseq = filin_read.split("Seq_ID=n_identical/n_aligned=")[1].split("\n")[0]
    try:
        IDseq = IDseq.replace(" ", "")
    except:
        pass

    return IDseq


def retrieveTMScore(filin_read):
    TMScore1 = filin_read.split("TM-score= ")[1].split(" (")[0]
    TMScore2 = filin_read.split("TM-score= ")[2].split(" (")[0]

    return [TMScore1, TMScore2]


# print parseOutputTMalign("/home/borrel/Yue_project/alignment/AMP/1Z8D__3E3O_D/RMSD")