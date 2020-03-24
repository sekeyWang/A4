from database import Database
from spectrum import Spectra, Spectrum
mass_AA = {'A': 71.03711,  # 0
           'R': 156.10111,  # 1
           'N': 114.04293,  # 2
           'D': 115.02694,  # 3
           'C': 160.03065,  # C(+57.02)
           'E': 129.04259,  # 5
           'Q': 128.05858,  # 6
           'G': 57.02146,  # 7
           'H': 137.05891,  # 8
           'I': 113.08406,  # 9
           'L': 113.08406,  # 10
           'K': 128.09496,  # 11
           'M': 131.04049,  # 12
           'F': 147.06841,  # 13
           'P': 97.05276,  # 14
           'S': 87.03203,  # 15
           'T': 101.04768,  # 16
           'W': 186.07931,  # 17
           'Y': 163.06333,  # 18
           'V': 99.06841,  # 19
           }

def mass_peptide(peptide):
    mass = 0
    for AA in peptide:
        mass += mass_AA[AA]
    return mass

def cal_ppm(mass1, mass2):
    return abs(mass1 - mass2)*1e6 / mass2

def cal_mass_tol(mass1, mass2):
    return abs(mass1 - mass2)

def find_peptide(sequence, theo_mass, mass_tol):
    L = len(sequence)
    en = 0
    current_mass = mass_AA[sequence[0]]
    peptide_list = []
    for st in range(0, L):
        while(en < L - 1 and current_mass < theo_mass):
            en += 1
            current_mass += mass_AA[sequence[en]]
        if cal_mass_tol(current_mass - mass_AA[sequence[en]], theo_mass) < mass_tol:
            peptide_list.append(sequence[st : en])
        if cal_mass_tol(current_mass, theo_mass) < mass_tol:
            peptide_list.append(sequence[st : en+1])
        current_mass -= mass_AA[sequence[st]]
#    print(sequence, peptide_list)
    return peptide_list


def find_candidate(database:Database, spectrum:Spectrum):
    theo_mass = spectrum.pepmass * spectrum.charge - 1.0073 * spectrum.charge - 18.0105
    mass_error = 0.1
    candidate = []
    for protein in database:
        candidate = candidate + find_peptide(protein.peptide, theo_mass, mass_error)
#    for seq in candidate:
#        print(seq, round(mass_peptide(seq) - theo_mass, 3))
    return candidate

if __name__ == '__main__':
    database = Database("ups.fasta")
    spectra = Spectra("test.mgf")
    find_candidate(database, spectra[61])