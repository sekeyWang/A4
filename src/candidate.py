from database import Database
from spectrum import Spectra, Spectrum
import config


def find_candidate(database:Database, spectrum:Spectrum):
    theo_mass = spectrum.pepmass * spectrum.charge - config.mass_z * spectrum.charge - config.mass_H2O
    mass_error = config.precursor_mass_tol
    candidate = []
    for peptide in database:
#        candidate = candidate + find_peptide(protein.peptide, theo_mass, mass_error)
        if config.cal_mass_tol(config.mass_peptide(peptide.sequence), theo_mass) < mass_error:
            candidate.append(peptide)
#    for seq in candidate:
#        print(seq, round(mass_peptide(seq) - theo_mass, 3))
    return candidate


if __name__ == '__main__':
    database = Database("../data/ups.fasta")
    spectra = Spectra("../data/test.mgf")

#    print(find_candidate(database, spectra[10]))
    match = 0
    for spectrum in spectra:
        if len(find_candidate(database, spectrum)) > 0:
            match += len(find_candidate(database, spectrum))
    print(match)