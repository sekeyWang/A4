from math import log10 as log
from database import Database
from spectrum import Spectra, Spectrum
import config
import statistics

def theo_ion_mass(peptide, ion, z):
    l = len(peptide)
    mass_list = []
    if ion == 'y':
        fragment_mass = config.mass_H2O
        for idx in range(l - 1):
            AA_mass = config.mass_AA[peptide[l - idx - 1]]
            fragment_mass += AA_mass
            mass_list.append((fragment_mass + config.mass_z) / z)

    if ion == 'b':
        fragment_mass = 0
        for idx in range(l - 1):
            AA_mass = config.mass_AA[peptide[idx]]
            fragment_mass += AA_mass
            mass_list.append((fragment_mass + config.mass_z) / z)

    return mass_list


def ion_match_score(mass_list, spectrum:Spectrum):
    mz_list = spectrum.mz_list
    it_list = spectrum.intensity_list
    idx = 0
    l = len(mz_list)
    score_list = []
#    filter = sorted(spectrum.intensity_list, reverse=True)[int(0.5*len(spectrum.intensity_list))]
#    filter = min(filter, statistics.mean(spectrum.intensity_list) + statistics.stdev(spectrum.intensity_list) * 2)
    filter = 0.01 * max(spectrum.intensity_list)
    for ion_mass in mass_list:
        while idx < l - 1 and mz_list[idx] < ion_mass: idx += 1
        peak_idx = idx
        if idx > 0 and mz_list[idx] - ion_mass > ion_mass - mz_list[idx - 1]:
            peak_idx = idx - 1

        peak_mz = mz_list[peak_idx]
        score = 0
        if abs(peak_mz - ion_mass) < config.fragment_mass_tol:
            peak_it = it_list[peak_idx]
            if peak_it > filter:
                score = log(peak_it * 100)
        score_list.append(score)
    return score_list


def llr_score(score_list, p, q):
    score = 0
    for x in score_list:
        if x > 0:
            score += log(p / q) * x
        else:
            score += log((1 - p) / (1 - q))
    return score


def match_score(peptide, spectrum:Spectrum):
    b1 = theo_ion_mass(peptide.sequence, 'b', 1)
#    b2 = theo_ion_mass(peptide.sequence, 'b', 2)
    y1 = theo_ion_mass(peptide.sequence, 'y', 1)
#    y2 = theo_ion_mass(peptide.sequence, 'y', 2)
    score_b1 = ion_match_score(b1, spectrum)
#   score_b2 = ion_match_score(b2, spectrum)
    score_y1 = ion_match_score(y1, spectrum)
#    score_y2 = ion_match_score(y2, spectrum)

#    score = llr_score(score_y1, config.yp, config.yq) + llr_score(score_b1, config.bp, config.bq)
    score = sum(score_y1)
#    print(spectrum.title[6:], peptide, score_y1, score)
    return score



if __name__ == '__main__':
    database = Database("../data/ups.fasta")
    spectra = Spectra("../data/test.mgf")
    for test_spectrum in spectra:
        print(test_spectrum.title)
        peptide_list = test_spectrum.find_candidate(database)
        for peptide in peptide_list:
            print(peptide.sequence, match_score(peptide, test_spectrum))
