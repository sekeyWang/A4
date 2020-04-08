class Protein:
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence
    def __repr__(self):
        return "Name:%s\nPeptide:%s" % (self.name, self.sequence)


class Database:
    def __init__(self, file_name):
        self.database = self.read_fasta(file_name)

    def __getitem__(self, item):
        return self.database[item]

    def read_fasta(self, file_name):
        fr = open(file_name, 'r')
        l = fr.readline()
        name, protein = "", ""
        database = []

        while(l):
            if (l.startswith('>')):
                if name != "":
                    peptide_list = self.trypsin_cut(protein)
                    for peptide in peptide_list:
                        database.append(Protein(name, peptide))
                name = l[:min(10,len(l))]
                protein = ""
            else:
                protein += l[:-1]
            l = fr.readline()
        peptide_list = self.trypsin_cut(protein)
        for peptide in peptide_list:
            database.append(Protein(name, peptide))
        return database

    def trypsin_cut(self, protein):
        peptide_list = []
        peptide = ""
        for idx, aa in enumerate(protein):
            peptide += aa
            if aa == 'R' or aa == 'K':
                if idx + 1 < len(protein) and protein[idx + 1] != 'P':
                    peptide_list.append(peptide)
                    peptide = ""
        return peptide_list

if __name__ == '__main__':
    database = Database("ups.fasta")
    for peptide in database:
        print(peptide)
