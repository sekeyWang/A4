from database import Database
from spectrum import Spectra, Spectrum
from candidate import find_candidate
from score import match_score

if __name__ == '__main__':
    database = Database("ups.fasta")
    spectra = Spectra("test.mgf")
    output_file = "out.csv"
    with open(output_file, 'w') as writer:
        writer.write("Id m/z z score peptide protein\n")
        for id, spectrum in enumerate(spectra):
            print(spectrum.title)
            peptide_list = find_candidate(database, spectrum)
            if len(peptide_list) == 0:
                continue
            max_score = 0
            best_peptide = peptide_list[0]
            for peptide in peptide_list:
                score = match_score(peptide, spectrum)
                if score > max_score:
                    max_score = score
                    best_peptide = peptide
            print(peptide.sequence, match_score(peptide, spectrum))
            writer.write("%d %f %d %f %s %s\n"%(id, spectrum.pepmass, spectrum.charge, max_score, peptide.sequence, peptide.name))