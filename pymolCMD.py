from pymol import cmd, finish_launching


def saveMultiplepng(pnamein, nbstates):

    finish_launching()
    for state in range(0,nbstates):
        print(state, nbstates)
        finish_launching()

        pfilin = pnamein + "_state" + str(state) + ".pdb"
        print (pfilin)
        cmd.load(pfilin, "state" + str(state))
        cmd.hide("everything", "all")



    for state in range(0, nbstates):
        cmd.show("surface", selection="state" + str(state))
        cmd.set("transparency", value=0.4)
        cmd.show("cartoon", selection="state" + str(state))
        cmd.spectrum("b", "chocolate_firebrick_red_tv_red_deepsalmon_darksalmon_salmon_lightpink_white",
                     selection="all", minimum=0, maximum=100)
        cmd.select("lig" + str(state), "resn UNK and state" + str(state))
        cmd.show("stick", selection="lig" + str(state))
        cmd.bg_color("white")

        cmd.ray(2000, 2000)
        cmd.save(pnamein + "_state" + str(state) + ".png", selection="state" + str(state))
        cmd.hide("everything", "all")

    cmd.quit()



saveMultiplepng("/home/aborrel/imitanib/results/analysis/MD_FPI/tempFPI/CHEMBL207986_2hyy_MD/FPItanimoto/protBfact",10)