import heapq
import math
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

def id_rev(id:int):
    return id^0b10
def stream_kmers(seq:str, k:int):
    # Initialize the kmer and its reverse complement
    kmer = 0
    rkmer = 0

    # Add the first k-1 nucleotides to the first kmer and its reverse complement
    for i in range(k-1):
        nucl = get_id_nucleotide(seq[i])
        rnucl = id_rev(nucl)
        kmer |= nucl << (2*(k-2-i))
        rkmer |= rnucl << (2*(i+1))

    mask = (1 << (2*(k-1))) - 1

    # Yield the kmers
    for i in range(k-1, len(seq)):
        nucl = get_id_nucleotide(seq[i])
        rnucl = id_rev(nucl)
 
        # Remove the leftmost nucleotide from the kmer 
        kmer &= mask
        # Shift the kmer to make space for the new nucleotide
        kmer <<= 2
        # Add the new nucleotide to the kmer
        kmer |= nucl
        # Make space for the new nucleotide in the reverse kmer (remove the rightmost nucleotide by side effect)
        rkmer >>= 2
        # Add the new nucleotide to the reverse kmer
        rkmer |= rnucl << (2*(k-1))

        yield min(kmer, rkmer)
def argmax(s:list):
    i = 0
    maximum = s[0]
    for k,m in enumerate(s):
        if m > maximum:
            maximum = m
            i = k
    return i
def filter_smallest(stream,s:int):
	heap =[-math.inf] * s
	for kmer in stream:
		if heap[0]< -kmer:
		    heapq.heappushpop(heap,-kmer)
	return [-h for h in heap]


def xorshift64(x:int):
    x ^= (x << 13) & 0xFFFFFFFFFFFFFFFF
    x ^= (x >> 7 ) & 0xFFFFFFFFFFFFFFFF
    x ^= (x << 17) & 0xFFFFFFFFFFFFFFFF
    return x
    
def min_hash_sketch(seq:str,k:int,s:int,hash=xorshift64)->list:
    hash = filter_smallest(map(hash,stream_kmers(seq,k)),s)
    return hash
