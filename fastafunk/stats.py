def get_sequence(sequence, unaligned=True):
    if not isinstance(sequence, str):
        sequence = str(sequence.seq)
    if unaligned:
        sequence = sequence.replace("-","")
    return sequence

def remove_terminal_gaps(sequence):
    i = 0
    while sequence[i] in ["N","-"] and i < len(sequence):
        i += 1
    j = len(sequence)
    while j > 0 and sequence[j-1] in ["N","-"]:
        j -= 1
    return sequence[i:j]

def get_length(sequence, unaligned=True, exclude_terminal_gaps=True):
    sequence = get_sequence(sequence, unaligned)
    if exclude_terminal_gaps:
        sequence = remove_terminal_gaps(sequence)
    return len(sequence)

def get_number_missing_bases(sequence, unaligned=True, exclude_terminal_gaps=True):
    sequence = get_sequence(sequence, unaligned)
    if exclude_terminal_gaps:
        sequence = remove_terminal_gaps(sequence)
    return sequence.count("N")

def get_proportion_gaps(sequence, unaligned=True, exclude_terminal_gaps=True):
    length = get_length(sequence, unaligned, exclude_terminal_gaps)
    if length == 0:
        return 0.0
    gaps = get_number_missing_bases(sequence, unaligned, exclude_terminal_gaps)
    return float(gaps)/length

def get_stat(stat, sequence, unaligned=True, exclude_terminal_gaps=True):
    switcher = {
        "length": get_length(sequence, unaligned, exclude_terminal_gaps),
        "missing": get_number_missing_bases(sequence, unaligned, exclude_terminal_gaps),
        "gaps": get_proportion_gaps(sequence, unaligned, exclude_terminal_gaps)
    }
    return switcher.get(stat, "Invalid statistic")
