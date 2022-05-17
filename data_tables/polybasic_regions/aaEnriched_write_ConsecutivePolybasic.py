from pyfasta import Fasta
import numpy as np
import pandas as pd

#fasta = '/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/FastaSubsets/Scer_protein_RiboSeq.fa'
fasta = '/Users/KevinStein/Desktop/Lab/Bioinformatics/ReferenceData/FastaSubsets/Cele_protein_RiboSeq.fa'
#fasta = '/Users/KevinStein/Desktop/stallNOTkenyon_aggregates.fa'
outFN = '/Users/KevinStein/Desktop/DE6of6.csv'
window = 6
cutoff = 6

f = Fasta(fasta)
kern = np.ones((window,))

enriched = 0
sites = 0
data = []
names = []

for k in sorted(f.keys()):
    seq = f[k]
    seq_array = np.array(f[k], dtype='c')     
    a = (seq_array == 'D') | (seq_array == 'E') # get a boolean array with True where AA is one of 'K', 'R'
    mw = np.convolve(a, kern, mode = 'same')    # take a moving window of 'window' aa
    if mw.max() >= cutoff:     # if it doesn't have any windows with >= 'cutoff' aa, skip it.
        indices = np.where(mw >= cutoff)[0] 
        # for 3 consecutive
        #indices = [i for i in indices if (i-1 not in indices and i-2 not in indices and i-3 not in indices and i-4 not in indices and i-5 not in indices)]
        # for 4 or 5 consecutive
        #indices = [i-1 for i in indices if (i-1 not in indices and i-2 not in indices and i-3 not in indices and i-4 not in indices and i-5 not in indices)]
        # for 6 consecutive
        indices = [i-2 for i in indices if (i-1 not in indices and i-2 not in indices and i-3 not in indices and i-4 not in indices and i-5 not in indices)]
        sites += len(indices)
        enriched += 1
        for index in indices:
            data.append(index)
            names.append(k.split()[0])

print "\n%s sites in %i proteins are enriched out of %g sequences" % (sites, enriched, len(f.keys()))
df = pd.DataFrame(data=data, index=names)

df.to_csv(outFN)