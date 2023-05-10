from biopylib.read_sequence import readSeq
from biopylib.sequence import SQ
# import termcolor
# from termcolor import colored

# https://dev.to/zirkelc/extract-highlighted-text-from-a-book-using-python-e15

# # read_seq = readSeq('../sequences/NC_000932.faa')   # multiple
# read_seq = readSeq('../sequences/AAH12844.1.faa') # just the one
# sequence = read_seq.store()
# print(type(sequence))
# print(sequence.seq)

#read_seq = readSeq('../sequences/NC_000932.gb')
#sequence = read_seq.store()
# sequence.info()

# print(sequence[0].seq)

seq_str = 'TTCCACGTCGGGTTCCCACATTCCTG'
seq_SQ = SQ(seq_str,
             seq_type='dna',
             id='pd|1',
             origin='/name.fna',
             description='sequence description')

seq_SQ.info()

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


