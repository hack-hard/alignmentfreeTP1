import heapq
import math
from itertools import chain
import numpy as np
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
import numpy as np

def pushpop(heap:np.ndarray, item:int):
    """
    Pushes item onto the heap, then pops and returns the smallest item.
    
    Args:
    heap (np.array): NumPy array representing a min heap
    item: The item to push onto the heap
    
    Returns:
    The smallest item in the resulting heap
    """
    if heap.size == 0:
        return item
    
    if item < heap[0]:
        item, heap[0] = heap[0], item
        _siftdown(heap, 0)
    
    return item

def _siftdown(heap:np.ndarray, pos):
    endpos = len(heap)
    startpos = pos
    newitem = heap[pos]
    childpos = 2*pos + 1
    while childpos < endpos:
        rightpos = childpos + 1
        if rightpos < endpos and not heap[childpos] < heap[rightpos]:
            childpos = rightpos
        heap[pos] = heap[childpos]
        pos = childpos
        childpos = 2*pos + 1
    heap[pos] = newitem
    _siftup(heap, startpos, pos)

def _siftup(heap, startpos, pos):
    newitem = heap[pos]
    while pos > startpos:
        parentpos = (pos - 1) >> 1
        parent = heap[parentpos]
        if newitem < parent:
            heap[pos] = parent
            pos = parentpos
            continue
        break
    heap[pos] = newitem

def filter_smallest(stream,s:int,dtype=np.int64):
    heap = np.full(s,np.iinfo(dtype).max,dtype = dtype)
    for kmer in stream:
        pushpop(heap,kmer)
    return heap

def filter_smallest_list(stream,s:int):
    heap = [-np.inf]*s
    for kmer in stream:
        heapq.heappushpop(heap,-kmer)
    return [- h for h in heap]


def xorshift64(x:int):
    x ^= (x << 13) & 0xFFFFFFFFFFFFFFFF
    x ^= (x >> 7 ) & 0xFFFFFFFFFFFFFFFF
    x ^= (x << 17) & 0xFFFFFFFFFFFFFFFF
    return x
    
def min_hash_sketch(seq:list,k:int,s:int,hash=xorshift64)->list:
    hash = filter_smallest_list(map(hash,chain.from_iterable((stream_kmers(s,k)) for s in seq)),s)
    return hash

def make_smaller(x:int):
    #assume that x is a 64 bit integer that is "small".
    dominant_bit = 63
    while x & ( 1 << dominant_bit) == 0 and dominant_bit > 0:
        dominant_bit -=1
    return (dominant_bit << 10) | (x & (1 <<10 -1))       

def compressed_min_hash_sketch(seq:list,k:int,s:int,hash=xorshift64) -> np.ndarray:
    hash = filter_smallest_list(map(lambda x: make_smaller(hash(x)),chain.from_iterable((stream_kmers(s,k)) for s in seq)),s)
    return hash
