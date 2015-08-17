HASH_SIZE = 150
MIN_OVERLAP = 300
MAX_ERRORS = 100
ERROR_PENALTY = 10

class Alignment:
    def __init__(self, score, shift, errors, strand):
        self.score = score
        self.shift = shift
        self.errors = errors
        self.strand = strand

def compare_sequences_impl(seq1, seq2, max_errors):
    best_score = -100000
    shift = 0
    best_errors = 100000
    max_shift = len(seq1) - MIN_OVERLAP
    for i in xrange(max_shift):
        errors = 0
        max_score = min(len(seq1) - i, len(seq2))
        for j in xrange(min(len(seq1) - i, len(seq2))):
            if seq1[i + j] != seq2[j]:
                errors += 1
            if max_score - errors * ERROR_PENALTY < best_score:
                break
            if errors > max_errors:
                break

        if errors <= max_errors:
            score = max_score - errors * ERROR_PENALTY
            if score > best_score:
                shift = i
                best_score = score
                best_errors = errors

    return (best_score, shift, best_errors)

def compare_sequences(seq1, seq2, max_errors):
    score1, shift1, errors1 = compare_sequences_impl(seq1, seq2, max_errors)
    score2, shift2, errors2 = compare_sequences_impl(seq2, seq1, max_errors)
    if score1 >= score2:
        return (score1, shift1, errors1)
    else:
        return (score2, -shift2, errors2)

def align_sequences(seq1, seq2, max_errors = MAX_ERRORS):
    if len(seq1) < MIN_OVERLAP or len(seq2) < MIN_OVERLAP:
        return Alignment(-100000, 0, 100000, True)
    score1, shift1, errors1 = compare_sequences(str(seq1), str(seq2), max_errors)
    score2, shift2, errors2 = compare_sequences(str(seq1), str(seq2.reverse_complement()), max_errors)
    if score1 >= score2:
        return Alignment(score1, shift1, errors1, True)
    else:
        return Alignment(score2, shift2, errors2, False)
