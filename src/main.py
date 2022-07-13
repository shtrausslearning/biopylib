from biopylib.read_sequence import readSeq
from biopylib.motif import dMotif
from biopylib.sequence import SQ
import termcolor
from termcolor import colored

# https://dev.to/zirkelc/extract-highlighted-text-from-a-book-using-python-e15

# # read_seq = readSeq('../sequences/NC_000932.faa')   # multiple
# read_seq = readSeq('../sequences/AAH12844.1.faa') # just the one
# sequence = read_seq.store()
# print(type(sequence))
# print(sequence.seq)

# read_seq = readSeq('../sequences/NC_000932.gb')
# sequence = read_seq.store()
# sequence.info()

# print(sequence[0].seq)

# seq_str = 'TTCCACGTCGGGTTCCCACATTCCTG'
# seq_SQ = SQ(seq_str,
#             seq_type='dna',
#             id='pd|1',
#             origin='/name.fna',
#             description='sequence description')
# seq_SQ.info()

# # Split the sequence at AA|CC
# cuts = seq_SQ.cut.cut_pattern(cut_pattern='GG|GT')
# print(cuts[0].info())
# # print(cuts[1].info())

# text = seq_SQ.seq
# l1=['CC']

# formattedText = []
# for t in text:
#     if t in l1:
#         print('pass')
#         formattedText.append(colored(t,'white','on_red'))
#     else: 
#         formattedText.append(t)

# print("".join(formattedText))

# seq_SQ.highlight('TTCC')
# sequence.features['1_gene']

# Read Case
    
# def test1():  
#     sm = dMotif(size=3,method='heuristic')
#     sm.read_aln("example.txt","dna")
#     result = sm.get()
#     print(result)
    
# test1()

# Defined

# Test Motif
# A[TAG]AGCTGA 1 
# ACG[TAG]ATGA 3
# AAGA[TAG]GGG 4

def test2():
    seq1 = SQ("ATAGAGCTGA","dna")
    seq2 = SQ("ACGTAGATGA","dna")
    seq3 = SQ("AAGATAGGGG","dna")
    mf = dMotif(size=3,
                method='exhaustive',
                seqs=[seq1,seq2,seq3])
    result = mf.get()
    print(result)
    print(mf.motif_size)
    
test2()


