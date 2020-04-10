from database import Database
import config


class Spectrum:
    title: str
    pepmass: float
    charge: int
    scans: str
    rt: float
    mz_list: list
    intensity_list: list

    def __init__(self, title, pepmass, charge, scans, rt, mz_list, intesity_list):
        self.title = title
        self.pepmass = pepmass
        self.charge = charge
        self.scans = scans
        self.rt = rt
        self.mz_list = mz_list
        self.intensity_list = intesity_list

    def __repr__(self):
        return "Title:%s\nPepmass:%s\nCharge:%s\nScans:%s\nRT:%s" % (self.title, self.pepmass, self.charge, self.scans, self.rt)

    def find_candidate(self, database: Database):
        theo_mass = self.pepmass * self.charge - config.mass_z * self.charge - config.mass_H2O
        mass_error = config.precursor_mass_tol
        candidate = []
        for peptide in database:
            if config.cal_mass_tol(config.mass_peptide(peptide.sequence), theo_mass) < mass_error * config.mass_z:
                candidate.append(peptide)
        return candidate


class Spectra:
    def __init__(self, file_name):
        self.spectra = self.read_mgf(file_name)

    def __getitem__(self, item):
        return self.spectra[item]

    def read_mgf(self, file_name):
        fr = open(file_name, 'r')
        spectra = []
        flag = False
        title, pepmass, charge, scans, rt, mz_list, intensity_list = "", 0, 0, 0, 0, [], []
        for l in fr.readlines():
            if l.startswith("BEGIN IONS"):
                title, pepmass, charge, scans, rt, mz_list, intensity_list = "", 0, 0, 0, 0, [], []
            elif l.startswith("TITLE="):
                title = l[6:-1]
            elif l.startswith("PEPMASS="):
                pepmass = float(l[8:-1])
            elif l.startswith("CHARGE="):
                charge = int(l[7:-2])
            elif l.startswith("SCANS="):
                scans = l[6:-1]
            elif l.startswith("RTINSECONDS="):
                rt = l[11:-1]
                flag = True
            elif l.startswith("END IONS"):
                spectra.append(Spectrum(title, pepmass, charge, scans, rt, mz_list, intensity_list))
                flag = False
            elif flag:
                array = l[:-1].split(' ')
                mass, intensity = float(array[0]), float(array[1])
                mz_list.append(mass)
                intensity_list.append(intensity)
        return spectra

if __name__ == '__main__':
    spectra = Spectra("../data/test.mgf")
    s = spectra[-1]
    print(s)
    print(len(s.mz_list), s.mz_list)
    print(len(s.intensity_list), s.intensity_list)