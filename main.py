import tableParse

















##########
#  MAIN  #
##########


pCHEMBL = "/home/aborrel/imitanib/CHEMBL/bioactivity-TK-ABL_CHEMBL1862.txt"


table = tableParse.CHEMBL(pCHEMBL)
table.parseCHEMBLFile()
print len(table.table)

table.getOnlyExactConstant()
print len(table.table)

table.delIdenticCHEMBLID()
print len(table.table)

table.writeTable("/home/aborrel/imitanib/CHEMBL/bioactivity-TK-ABL_CHEMBL1862_filtered.txt")