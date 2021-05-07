# This is a backup file for unused algorithms


def parse_to_allele(sequence, window_size):
    allele = []
    for i, sequence in enumerate(sequence):
        num = len(sequence)//window_size
        allele.append([]*num)
        for index in range(num):
            allele[i].append(sequence[index*window_size:(index + 1)*window_size])
        allele[i].append(sequence[num*window_size:])
    return allele
