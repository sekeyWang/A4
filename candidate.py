from database import Database
from spectrum import Spectra, Spectrum
import config

def find_peptide(sequence, theo_mass, mass_tol):#
    L = len(sequence)
    en = 0
    current_mass = config.mass_AA[sequence[0]]
    peptide_list = []
    for st in range(0, L):
        while(en < L - 1 and current_mass < theo_mass):
            en += 1
            current_mass += config.mass_AA[sequence[en]]
        if config.cal_mass_tol(current_mass - config.mass_AA[sequence[en]], theo_mass) < mass_tol:
            peptide_list.append(sequence[st : en])
        if config.cal_mass_tol(current_mass, theo_mass) < mass_tol:
            peptide_list.append(sequence[st : en+1])
        current_mass -= config.mass_AA[sequence[st]]
#    print(sequence, peptide_list)
    return peptide_list


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
    database = Database("ups.fasta")
    spectra = Spectra("test.mgf")

#    print(find_candidate(database, spectra[10]))
    match = 0
    for spectrum in spectra:
        if len(find_candidate(database, spectrum)) > 0:
            match += len(find_candidate(database, spectrum))
    print(match)