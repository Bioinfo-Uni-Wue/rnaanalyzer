import re

def correctgensq(seq):
    """
    Takes a nucleotide sequence and checks it for incorrect bases.
    Returns the corrected sequence.
    """
    # Remove non-alphabet characters
    seq = re.sub(r'[^a-zA-Z]', '', seq)
    # Convert to lowercase
    seq = seq.lower()
    # Replace non-standard nucleotides with 'n'
    seq = re.sub(r'[^agtuc]', 'n', seq)
    # Replace 't' with 'u'
    seq = seq.replace('t', 'u')
    return seq

# Example usage
sequence = "aBcDXYZ123!@#agtTTT"
corrected_sequence = correctgensq(sequence)
print(f"Corrected Sequence: {corrected_sequence}")
