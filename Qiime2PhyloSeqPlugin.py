import PyPluMA

def isNumeric(s):
    for i in range(len(s)):
        if ((not s[i].isdigit()) and (s[i] != '.')):
            return False
    return True

class Qiime2PhyloSeqPlugin:
    def input(self, filename):
       infile = open(filename, 'r')
       self.parameters = dict()
       for line in infile:
          contents = line.strip().split('\t')
          self.parameters[contents[0]] = contents[1]
       self.otu_table = []
       otufile = open(PyPluMA.prefix()+"/"+self.parameters['otufile'], 'r')
       for line in otufile:
          self.otu_table.append(line.strip().split('\t'))
       self.mapping_table = []
       mapfile = open(PyPluMA.prefix()+"/"+self.parameters['mapping'], 'r')
       for line in mapfile:
           self.mapping_table.append(line.strip().split('\t'))

    def run(self):
       self.tax_table = []
       self.tax_table.append(["","Rank1","Rank2","Rank3","Rank4","Rank5","Rank6","Rank7"])
       for i in range(1, len(self.otu_table)):
          taxon = self.otu_table[i][0]  # First entry
          level = "Unassigned"
          #print(taxon)
          if (taxon.startswith("Unassigned")):
              taxarow = ["OTU"+str(i),"Unassigned","Unassigned_unclassified","Unassigned_unclassified","Unassigned_unclassified","Unassigned_unclassified","Unassigned_unclassified","Unassigned_unclassified"]
          else:
             taxarow = ["OTU"+str(i)]
             # KINGDOM
             if (taxon.find("k__") != -1):
                 kingdom = taxon[taxon.index("k__")+3:taxon.index(";")]
                 if (len(kingdom) != 0):
                    taxarow.append(kingdom)
                    level = kingdom
                 else:
                    taxarow.append(level+"_unclassified")
             else:
                 taxarow.append(level+"_unclassified")
             taxon = taxon[taxon.index(";")+1:]
             
             # PHYLUM
             if (taxon.find("p__") != -1):
                 phylum = taxon[taxon.index("p__")+3:taxon.index(";")]
                 if (len(phylum) != 0):
                    taxarow.append(phylum)
                    level = phylum
                 else:
                    taxarow.append(level+"_unclassified")
             else:
                 taxarow.append(level+"_unclassified")
             taxon = taxon[taxon.index(";")+1:]
             
             # CLASS
             if (taxon.find("c__") != -1):
                 # class is a python keyword :\
                 itsclass = taxon[taxon.index("c__")+3:taxon.index(";")]
                 if (len(itsclass) != 0):
                    taxarow.append(itsclass)
                    level = itsclass
                 else:
                    taxarow.append(level+"_unclassified")
             else:
                 taxarow.append(level+"_unclassified")
             taxon = taxon[taxon.index(";")+1:]

             # ORDER
             if (taxon.find("o__") != -1):
                 order = taxon[taxon.index("o__")+3:taxon.index(";")]
                 if (len(order) != 0):
                    taxarow.append(order)
                    level = order
                 else:
                    taxarow.append(level+"_unclassified")
             else:
                 taxarow.append(level+"_unclassified")
             taxon = taxon[taxon.index(";")+1:]

             # FAMILY
             if (taxon.find("f__") != -1):
                 family = taxon[taxon.index("f__")+3:taxon.index(";")]
                 if (len(family) != 0):
                     taxarow.append(family)
                     level = family
                 else:
                     taxarow.append(level+"_unclassified")
             else:
                 taxarow.append(level+"_unclassified")
             taxon = taxon[taxon.index(";")+1:]

             # GENUS
             if (taxon.find("g__") != -1):
                 genus = taxon[taxon.index("g__")+3:taxon.index(";")]
                 if (len(genus) != 0):
                     taxarow.append(genus)
                     level = genus
                 else:
                     taxarow.append(level+"_unclassified")
             else:
                 taxarow.append(level+"_unclassified")
             taxon = taxon[taxon.index(";")+1:]

             # SPECIES
             if (taxon.find("s__") != -1):
                 species = taxon[taxon.index("s__")+3:]
                 if (len(species) != 0):
                     taxarow.append(species)
                 else:
                     taxarow.append(level+"_unclassified")
             else:
                 taxarow.append(level+"_unclassified")

          for i in range(0, len(taxarow)):
              taxarow[i] = taxarow[i].replace('[', '')
              taxarow[i] = taxarow[i].replace(']', '')
          self.tax_table.append(taxarow)
       #print(self.otu_table[0][0])
       #self.otu_table[0][0] = ''
       self.otu_table[0] = ['',] + self.otu_table[0]
       for i in range(1, len(self.otu_table)):
          self.otu_table[i][0] = "OTU"+str(i)



    def output(self, filename):
        otufile = open(filename+".otu.csv", 'w')
        taxfile = open(filename+".tax.csv", 'w')
        mapfile = open(filename+".data.csv", 'w')
        header = self.otu_table[0]
        for i in range(len(header)):
            otufile.write('\"'+header[i]+'\"')
            if (i != len(header)-1):
                otufile.write(',')
            else:
                otufile.write('\n')

        for j in range(1, len(self.otu_table)):
            abundances = self.otu_table[j]
            otufile.write('\"'+abundances[0]+'\",')
            for i in range(1, len(abundances)):
                otufile.write(abundances[i])
                if (i != len(abundances)-1):
                    otufile.write(',')
                else:
                    otufile.write('\n')
        for taxon in self.tax_table:
            for i in range(len(taxon)):
                taxfile.write('\"'+taxon[i]+'\"')
                if (i != len(taxon)-1):
                    taxfile.write(',')
                else:
                    taxfile.write('\n')
        header = self.mapping_table[0]
        mapfile.write('\"Sample\",')
        for i in range(1, len(header)):
            mapfile.write('\"'+header[i]+'\"')
            if (i != len(header)-1):
                mapfile.write(',')
            else:
                mapfile.write('\n')
        for i in range(1, len(self.mapping_table)):
            myrow = self.mapping_table[i]
            for j in range(len(myrow)):
                if (isNumeric(myrow[j])):
                    mapfile.write(myrow[j])
                else:
                    mapfile.write('\"'+myrow[j]+'\"')
                if (j != len(myrow)-1):
                    mapfile.write(',')
                else:
                    mapfile.write('\n')

