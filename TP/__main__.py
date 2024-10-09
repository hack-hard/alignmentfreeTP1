from TP.loading import load_directory
from TP.kmers import stream_kmers, kmer2str
from itertools import product,chain
from statistics import mean

def jaccard(seq_a:list,seq_b:list):
    kmers_a = list(chain.from_iterable(stream_kmers(a,k) for a in seq_a))
    kmers_b = list(chain.from_iterable(stream_kmers(b,k) for b in seq_b))
    kmers_a.sort()
    kmers_b.sort()
    i ,j = 0,0
    n_unions =0
    while i < len(kmers_a) or j < len(kmers_b):
        if kmers_a[i] == kmmers_b[j]:
            n_union +=1
            i+=1
            j+=1
        elif kmers_a[i] < kmers_b[j]:
            i+=1
        else:
            j+=1
    return n_unions /(len(kmers_a) + len(kmers_b) - n_unions)

if __name__ == "__main__":
    print("Computation of Jaccard similarity between files")

    # Load all the files in a dictionary
    files = load_directory("data")
    k = 21
    
    print("Computing Jaccard similarity for all pairs of samples")
    filenames = list(files.keys())
    for i in range(len(files)):
        for j in range(i+1, len(files)):
            
            # --- Complete here ---

            j = jaccard(files[filenames[i]], files[filenames[j]], k)
            print(filenames[i], filenames[j], j)
