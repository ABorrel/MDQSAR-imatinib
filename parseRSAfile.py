from re import search



class RSA:
    def __init__(self, pRSA):
        self.pRSA = pRSA
        self.buildClass()

    def buildClass(self):

        filin = open(self.pRSA, "r")
        llinesRSA = filin.readlines()
        filin.close()

        self.lres = []
        for lineRSA in llinesRSA:
            if search("^RES", lineRSA):
                self.lres.append(RESrsa(lineRSA))


class RESrsa:
    def __init__(self, lineRSA):

        if lineRSA != "":
            self.recoder = lineRSA[0:3].replace(" ", "")
            self.resName = lineRSA[4:7].replace(" ", "")
            self.chainID = lineRSA[8]
            self.resSeq = lineRSA[10:15].replace(" ", "")
            self.ABSall = lineRSA[16:22].replace(" ", "")
            self.RELall = lineRSA[23:28].replace(" ", "")
            self.ABSside = lineRSA[29:35].replace(" ", "")
            self.RELside = lineRSA[36:41].replace(" ", "")
            self.ABSmain = lineRSA[42:48].replace(" ", "")
            self.RELmain = lineRSA[49:54].replace(" ", "")
            self.ABSnonpolar = lineRSA[55:61].replace(" ", "")
            self.RELnonpolar = lineRSA[62:67].replace(" ", "")
            self.ABSpolar = lineRSA[68:74].replace(" ", "")
            self.RELpolar = lineRSA[75:80].replace(" ", "")

        print self.recoder, self.resName, self.chainID, self.resSeq, self.ABSall, self.RELall, self.ABSside, self.ABSmain, self.RELmain
        print self.ABSnonpolar, self.RELnonpolar, self.ABSpolar, self.RELpolar
