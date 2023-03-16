"""
coding=utf-8
author = Alejandro González

"""


def is_seq_valid(seq):
    valid_bases = ['A', 'T', 'G', 'C', 'R', 'Y', 'K', 'M', 'B', 'V', 'D', 'H', 'N', 'S', 'W']
    for base in seq:
        if base not in valid_bases:
            return False
    return True


def reverse(seq):
    return seq[::-1]


def complement(seq):
    complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'R': 'Y', 'Y': 'R', 'K': 'M', 'M': 'K',
                       'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N', 'S': 'S', 'W': 'W'}
    seq_list = list(seq)
    seq_list = [complement_dict[base] for base in seq_list]

    return ''.join(seq_list)


def reverse_complement(seq):
    seq = reverse(seq)
    seq = complement(seq)
    return seq


def main(sec):
    # Prompt user here
    sequence = sec.upper()

    # Loop until invalid char
    result = is_seq_valid(sequence)
    if not result:
        print(f'{sequence} no es una secuencia valida')
        return None

    return reverse_complement(sequence)


def main_v2(sec):
    print('---------------------- Reverse Complement Program ----------------------')
    # Prompt user here
    sequence = sec.upper()

    # Loop until invalid char
    result = is_seq_valid(sequence)
    if not result:
        print(f'{sequence} no es una secuencia valida')
        return None

    print('Secuencia original: ' + sequence)
    print('Reverso Complementario: ' + reverse_complement(sequence))

    print('------------------------------------------------------------------------')
