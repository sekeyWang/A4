class Protein:
    def __init__(self, name, peptide):
        self.name = name
        self.peptide = peptide
    def __repr__(self):
        return "Name:%s\nPeptide:%s" % (self.name, self.peptide)

class Database:
    def __init__(self, file_name):
        self.database = self.read_fasta(file_name)
    def __getitem__(self, item):
        return self.database[item]

    def read_fasta(self, file_name):
        fr = open(file_name, 'r')
        l = fr.readline()
        name, peptide = "", ""
        database = []

        while(l):
            if (l.startswith('>')):
                if name != "":
                    database.append(Protein(name, peptide))
                name = l[:min(11,len(l))]
                peptide = ""
            else:
                peptide += l[:-1]
            l = fr.readline()
        database.append(Protein(name, peptide))
        return database
if __name__ == '__main__':
    database = Database("ups.fasta")
    print(database[50])
