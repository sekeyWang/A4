from database import Database
from spectrum import Spectra
from score import match_score
import sys


def search(database_name, spectra_name):
    database = Database(database_name)
    spectra = Spectra(spectra_name)
    result = []
    for id, spectrum in enumerate(spectra):
#        print(spectrum.title)
        peptide_list = spectrum.find_candidate(database)
        if len(peptide_list) == 0:
            continue
        max_score = -1000000
        best_peptide = peptide_list[0]
        for peptide in peptide_list:
            score = match_score(peptide, spectrum)
            if score > max_score:
                max_score = score
                best_peptide = peptide
#        print(best_peptide.sequence, match_score(best_peptide, spectrum))
        result.append({
            "Id": spectrum.title[6:],
            "m/z": spectrum.pepmass,
            "z": spectrum.charge,
            "score": max_score,
            "peptide": best_peptide.sequence,
            "protein": best_peptide.name
        })
    return result


def write_result(result, output_file = "out.csv"):
    with open(output_file, 'w') as writer:
        writer.write("Id m/z z score peptide protein\n")
        for r in result:
            row = ""
            for item in r:
                row += str(r[item]) + ' '
            row += '\n'
            writer.write(row)


def main(spectra_name, database_name, output_file):
    result = search(database_name, spectra_name)
    result = sorted(result, key=lambda t:t["score"], reverse=True)
    write_result(result, output_file)


if __name__ == '__main__':
    spectra_name, database_name, output_file = sys.argv[1], sys.argv[2], sys.argv[3]
    main(spectra_name, database_name, output_file)