from database import Database
from spectrum import Spectra, Spectrum
from candidate import find_candidate


def match_score(peptide, spectrum:Spectrum):
    print(peptide, spectrum.mz_list)


if __name__ == '__main__':
    database = Database("ups.fasta")
    spectra = Spectra("test.mgf")
    test_spectrum = spectra[61]
    peptide_list = find_candidate(database, test_spectrum)
    for peptide in peptide_list:
        print(match_score(peptide, test_spectrum))
