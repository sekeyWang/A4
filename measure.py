def measure_result(file_name):
    with open(file_name, 'r') as fr:
        decoy, target, max_n = 0, 0, 0
        for idx, l in enumerate(fr.readlines()):
            row = l.split(' ')
            if row[5].startswith('>de'):
                decoy += 1
            else:
                target += 1
            if decoy / target < 0.01:
                max_n = idx + 1
    return max_n


if __name__ == '__main__':
    print(measure_result("out.csv"))