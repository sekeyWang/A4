import config


class Protein:
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

    def __repr__(self):
        return "Name:%s\nPeptide:%s" % (self.name, self.sequence)


class Fasta:
    def __init__(self, file_name):
        self.fasta = self.read_file(file_name)

    def __getitem__(self, item):
        return self.fasta[item]

    def read_file(self, file_name):
        fr = open(file_name, 'r')
        name, protein = "", ""
        fasta = []
        for l in fr.readlines():
            if l.startswith('>'):
                if name != "":
                    peptide_list = self.error_cut(protein)
                    for peptide in peptide_list:
                        fasta.append(Protein(name, peptide))
                name = l[:-1]
                protein = ""
            else:
                protein += l[:-1]
        peptide_list = self.error_cut(protein)
        for peptide in peptide_list:
            fasta.append(Protein(name, peptide))
        return fasta

    def error_cut(self, protein):
        peptide_list = []
        st = 0
        for idx, aa in enumerate(protein):
            if aa not in config.mass_AA:
                peptide_list.append(protein[st:idx])
                st = idx + 1
        if st < len(protein):
            peptide_list.append(protein[st:])
        return peptide_list


class Database:
    def __init__(self, file_name):
        self.database = self.read_fasta(Fasta(file_name))

    def __getitem__(self, item):
        return self.database[item]

    def read_fasta(self, fasta:Fasta):
        database = []
        for protein in fasta:
            name = protein.name[:10]
            peptide_list = self.trypsin_cut(protein.sequence)
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
    database = Database("../data/ups_decoy.fasta")
    for peptide in database:
        print(peptide)
