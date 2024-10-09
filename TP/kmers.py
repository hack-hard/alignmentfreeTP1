
def kmer2str(val, k):
    """ Transform a kmer integer into a its string representation
    :param int val: An integer representation of a kmer
    :param int k: The number of nucleotides involved into the kmer.
    :return str: The kmer string formatted
    """
    letters = ['A', 'C', 'T', 'G']
    str_val = []
    for _ in range(k):
        str_val.append(letters[val & 0b11])
        val >>= 2

    str_val.reverse()
    return "".join(str_val)

def get_id_nucleotide(s:str):
    return (ord(s) >> 1) & 0b11

def get_id_rev(s:str):
    return get_id_nucleotide(s) ^0b10

def stream_kmers(text:str, k:int):
    assert k <= len(text) 
    res = 0
    rev_res = 0
    print(f"streaning with k = {k} and {len(text)} nucleotides")
    for i in range(k-1):
        res |= get_id_nucleotide(text[i])
        rev_res |= get_id_rev(text[i]) << 2*(k-1)
        res <<=2
        rev_res >>=2
    for i in range(k,len(text)):
        res |= get_id_nucleotide(text[i])
        rev_res |= get_id_rev(text[i]) <<  2*(k-1)
        # print(f"string {text[i-k+1:i+1]}, got {kmer2str(res,k)} and {kmer2str(rev_res,k)}")
        yield min(res,rev_res)
        res <<=2
        rev_res >>=2
