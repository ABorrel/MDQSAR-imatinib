from os import path, system

def generateMovie(lfilePDB, pfiloutMovie):

    nbframes = len(lfilePDB)
    pbasedir = path.dirname(lfilePDB[0]) + "/" + lfilePDB[0].split("/")[-1].split("_")[0]

    # change pymol script
    filescript = open("pymolCMD.py", "r")
    llinesscript = filescript.readlines()
    print llinesscript
    filescript.close()

    filescript = open("pymolCMD.py", "w")
    filescript.write("".join(llinesscript[:-1]))
    filescript.write("saveMultiplepng(\"" + str(pbasedir) + "\"," + str(nbframes) + ")")
    filescript.close()

    system("pymol pymolCMD.py")

    # command for convert png to mpeg
    cmdMPEG = "mencoder " + str(pbasedir) + "*png -mf type=png:fps=1 -opt 200 -ovc lavc -o " + str(pfiloutMovie)
    system(cmdMPEG)
    print cmdMPEG
