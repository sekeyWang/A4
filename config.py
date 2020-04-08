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
mass_z = 1.0073
mass_H2O = 18.0105
precursor_mass_tol = 0.1
fragment_mass_tol = 0.5


def mass_peptide(peptide):
    mass = 0
    for AA in peptide:
        mass += mass_AA[AA]
    return mass


def cal_ppm(mass1, mass2):
    return abs(mass1 - mass2)*1e6 / mass2


def cal_mass_tol(mass1, mass2):
    return abs(mass1 - mass2)