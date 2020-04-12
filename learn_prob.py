from database import Database
from spectrum import Spectra, Spectrum
from score import theo_ion_mass, ion_match_score


def non_zero(arr):
    num = 0
    for x in arr:
        if x > 0:
            num += 1
    return num


def match_score(peptide, spectrum:Spectrum, ion):
    y1 = theo_ion_mass(peptide.sequence, ion, 1)
    score_y1 = ion_match_score(y1, spectrum)
    return score_y1


def learn_prob(database_name, spectra_name, ion):
    database = Database(database_name)
    spectra = Spectra(spectra_name)
    target_match, target_total = 0, 0
    decoy_match, decoy_total = 0, 0
    for id, spectrum in enumerate(spectra):
        peptide_list = spectrum.find_candidate(database)
        for peptide in peptide_list:
            score_list = match_score(peptide, spectrum, ion)
            match = non_zero(score_list)
            total = len(score_list)
            if peptide.name.startswith('>de'):
                decoy_match += match
                decoy_total += total
            else:
                target_match += match
                target_total += total
    print(ion + 'p=%.5f' % (target_match / target_total))
    print(ion + 'q=%.5f' % (decoy_match / decoy_total))


if __name__ == '__main__':
    learn_prob("data/ups_decoy.fasta", "data/test.mgf", 'b')