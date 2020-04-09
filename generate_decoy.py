import random
from database import Fasta


def generate_decoy(input_name, output_name):
    fasta = Fasta(input_name)
    with open(output_name, 'w') as writer:
        for protein in fasta:
            sequence = protein.sequence
            name = protein.name
            writer.write(name + '\n')
            writer.write(sequence + '\n')

            sequence = ''.join(random.sample(sequence, len(sequence)))
            name = ">de" + protein.name
            writer.write(name + '\n')
            writer.write(sequence + '\n')


if __name__ == '__main__':
    generate_decoy("ups.fasta", "ups_decoy.fasta")