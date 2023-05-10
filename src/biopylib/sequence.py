import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from re import finditer,sub,search
from colorama import Back, Style, Fore
from termcolor import colored
import random
import copy

'''

#####################################################

# General SQ Methods

# info - show sequence information
# view - visualise the sequence
# validate - check if the sequence contains no errors
# reverse_comp - reverse complement strand of DNA 
# add_annotation - add annotation to sequence

#####################################################

'''

class SQ: 
    
    # Constructor
    def __init__(self, seq=None,            # Sequence String
                        seq_type = "dna",   # Sequence Type 
                        id=None,            # Sequence Identifier (accessions)
                        origin=None,        # Sequence origin
                        description=None,   # Sequence Description
                        valid = None,        # Sequence validity status 
                        
                        locus=None,         # locus tag
                        accession=None,     # accession identifier
                        version=None,       # version of accession ID
                        dblink=None,        # database link
                        keywords=None,      # keywords
                        source=None,        # source
                        reference=None,     # references (dict)
                        comment=None,       # comments
                        features=None,      # stored features
                        ): 
        
        # core SQ
         
        self.seq = seq.upper()    # sequence in string format
        self.seq_type = seq_type  # sequence type
        
        # If more detail is provided
        
        self.id = id  # sequence identifier (eg. database id)
        self.description = description # give the sequence some description
        self.annotations = {} # Annotate some parts of the seqence
        
        # genbank specific
        
        self.locus = locus
        self.accession = accession
        self.version = version
        self.dblink = dblink
        self.keywords = keywords
        self.source = source
        self.reference = reference
        self.comment = comment
        self.features = features
        
        # if some fields are not filled
        if(self.id is None):
            self.id = f'sequence_{random.randint(1,1000)}'
        if(self.description is None):
            self.description = 'not given'
            
        self.valid = valid          # store sequence validity status 
        self.origin = origin       # sequence origin
        
        # If sequence was defined but origin wasn't provided (from read file)
        if(self.seq is not None and self.origin is None):
            self.origin = 'input'  # file origin, read/input/other sequence
        
        # Instantiate 
        self.decode = Decode(self.seq,self.seq_type) 
        self.count = Count(self.seq,self.seq_type,self.id)
        self.cut = Cut(copy.deepcopy(self))
        self.pattern = Pattern(self.seq,self.seq_type)
    
    # class instance operations
    
    def __len__(self):
        return len(self.seq)
    def __getitem__(self, n):
        return self.seq[n]
    def __getslice__(self, i, j):
        return self.seq[i:j]
    def __str__(self):
        return self.seq
    def __add__(self,other):
        if(self.seq_type == other.seq_type):
            return SQ(self.seq + other.seq,seq_type=self.seq_type)
        else:
            print('[note] sequences must of be same type')
            
    @staticmethod
    def colored(lseq):
        
        bcolors = {
            # standard color-codes
            'A': Back.GREEN,
            'a': Back.GREEN,
            'C': Back.YELLOW,
            'c': Back.YELLOW,
            'G': Back.RED,
            'g': Back.RED,
            'T': Back.BLUE,
            't': Back.BLUE,
            'reset': Style.RESET_ALL
        }

        tmpStr = ""
        for nuc in lseq:
            if nuc in bcolors:
                tmpStr += bcolors[nuc] + nuc
            else:
                tmpStr += bcolors['reset'] + nuc
        return tmpStr + '\033[0;0m'
    
    @staticmethod
    def split_nucleotide_seq(genome,cut_id):
        genes = []
        for ix, char in enumerate(genome):
            if ix != 0 and ix%cut_id == 0:
                genes.append(' ')
            genes.append(char)
        return ''.join(genes)
    
    def __seq_repr(self,genome_str, strand ='dna',split_id=10):
        
        nu_clr_switcher = {
            # standard color-codes
            'A': Back.GREEN,
            'a': Back.GREEN,
            'C': Back.YELLOW,
            'c': Back.YELLOW,
            'G': Back.RED,
            'g': Back.RED,
            'T': Back.BLUE,
            't': Back.BLUE,
            ' ': Style.RESET_ALL
        }
        
        if(strand == 'dna'):
            genome_str = self.split_nucleotide_seq(genome=genome_str,
                                                   cut_id=split_id)
            line_break_cntr = 0
            for i in range(len(genome_str)):
                if genome_str[i] == ' ':
                    line_break_cntr += 1
                    if line_break_cntr>0 and line_break_cntr%6==0:
                        text = "\n"
                    else:
                        text = nu_clr_switcher[genome_str[i]] + genome_str[i]
                else:
                    text = nu_clr_switcher[genome_str[i]] + genome_str[i]
                print(text, end="")
            Style.RESET_ALL
            
    # Print Sequence Information
    
    def info(self):
        
        print('[sequence information]')
        if(len(self.seq) > 10):
            
            tseq = self.seq[0:10]
            if(self.seq_type == 'dna' or self.seq_type == 'rna'):
                tseq = self.colored(tseq)
                
            tups =  (f"seq: {tseq} ...\n"
                     f"type: {self.seq_type}\n"
                     f"id: {self.id}\n"
                     f"valid: {self.valid}\n"
                     f"origin: {self.origin}\n"
                     f"description: {self.description}")
            
        else:
            
            # colourcode dna
            if(self.seq_type == 'dna' or self.seq_type == 'rna'):
                tseq = self.colored(self.seq)
            else:
                tseq = self.seq
            
            tups =  (f"seq: {self.colored(self.seq)}\n"
                     f"type: {self.seq_type}\n"
                     f"id: {self.id}\n"
                     f"valid: {self.valid}\n"
                     f"origin: {self.origin}\n"
                     f"description: {self.description}")
            
        print(tups);pass

    # convert sequence to pandas series

    def to_series(self):
        return pd.Series({'seq':self.seq,'type':self.seq_type,
                          'id':self.id,'valid':self.valid,
                          'origin':self.origin,'description':self.description},
                         name=self.id)
        
    # Highlight parts of the sequence

    def highlight(self,pattern=None):
        
        if(pattern != None):
            print(pattern in self.seq)
            hl = self.seq.replace(pattern,
                                  colored(pattern,
                                         'white','on_red'))
            print(hl)
        else:
            print('[note] enter pattern to highlight')
    
    # Visualise the nucleotide sequence
            
    def view(self,split_id=10):
        self.__seq_repr(self.seq,self.seq_type,split_id)
        
    # Check Validity
    
    def validate(self,verbose=False):
        alp = self.count.abc()
        res = True; i = 0
        while (res and i < len(self.seq)):
            if self.seq[i] not in alp: 
                res = False
            else: i += 1
        if(res):
            self.valid = True
            if(verbose):
                print(f'{self.seq_type} is valid')
            return res
        else:
            self.valid = False
            if(verbose):
                print(f'{self.seq_type} is invalid')
            return res
        
    # Reverse Complement of DNA molecule
    
    def reverse_comp(self):
        
        if (self.seq_type != "dna"): 
            print('[note] input only DNA sequence')
            return None
        
        lst_seq = ['A','T','G','C']
        lst_comp = ['T','A','C','G']
        
        comp = ''
        for char in self.seq:
            ii=-1
            for c in lst_seq:
                ii+=1
                if(char == c ):
                    comp = lst_comp[ii] + comp
                    
        return SQ(comp,"dna")
    
    #####################################################
            
    # Annotating Sequence, Notes
    # add_annotation - add_annotation
            
    #####################################################
    
    # method to add annotation to the current sequence
    def annotate(self,pattern,description=None):
        
        # input sequence (string/SQ)
        lseq = pattern      
        
        # If input pattern is a String
        
        if(type(pattern) == str):
            
            print('[note] creating bare SQ' )
            if(description is not None):
                lseq = SQ(pattern,
                          description=description)
            else:
                lseq = SQ(pattern)
                
            lid = lseq.id              # randomly created
            ldesc = lseq.description   # no description
                
        # Input pattern is a SQ sequence
                
        else:
            
            # just get the attributes from the class instance
            lid = lseq.id              # input sequence identifier
            ldesc = lseq.description   # input sequence description
                
        # Find the subsequence index in the main sequnce 
        # always returns a list, even if one found
        idx = self.pattern.find(pattern=str(lseq),
                                find_id='all',
                                search_id='standard',
                                verbose=False)
        
        # Add found input sequence pattern in main sequence
        
        if(idx != None):
            for i in idx:
                self.annotations[f'({i}...{i+len(lseq)})'] = f"{lid}_{ldesc}"
        else:
            print('[note] pattern not found in sequence')
            
    # Remove existing annotation
            
    def remove_annotation(self,annotation=None):
        
        if(annotation is not None):
            self.annotations.pop(annotation)
        else:
            print('[note] enter (annotation) to be removed')
            
'''

########################################################

Cutting the Sequence
cut_pattern - cut sequence based on particular pattern
dnaseq_features - cut the sequence 

########################################################

'''

class Cut:
    
    # Constructor
    def __init__(self,sq_id):
        self.__sq_id = sq_id # SQ class instance
            
    # converts IUB ambiguity code into RE
    # returns cut position of a restriction enzyme 
    # (in IUB code) in a sequence
            
    @staticmethod
    def divide_loc(enzyme, sequence):
        
        def iubrex(IUB):   
            
            # main 4 bases
            # purine, pyrimidine, amino
            # keto, strong, weak
            # not A, not C, not G, not T
            #  not T, any 
            
            dic = {"A":"A", "C":"C", "G":"G", "T":"T", 
                    # Additional Cases
                    "R":"[GA]", "Y":"[CT]", "M":"[AC]", 
                    "K":"[GT]", "S":"[GC]", "W": "[AT]",
                    "B":"[CGT]", "D":"[AGT]", "H":"[ACT]",
                    "V":"[ACG]", "N":"[ACGT]"}
            
            site = IUB.replace("|","")
            rex = ""
            
            for c in site:
                rex += dic[c]
            return rex
        
        regexp = iubrex(enzyme) # convert pattern to IUB format
        matches = finditer(regexp, sequence)
        locs = []
        for match in matches:
            locs.append(match.start() + enzyme.find("|"))
            
        return locs # indicies of cuts
    
    # determines subsequences resulting from a sequence 
    # cut in a list of positions
    def cut_pattern(self,cut_pattern=None):
        
        if(cut_pattern is not None):
        
            res = []
            positions = self.divide_loc(cut_pattern,self.__sq_id.seq)
            positions.insert(0,0)
            positions.append(len(self.__sq_id.seq))
            for i in range(len(positions)-1):
                lseq = self.__sq_id.seq[positions[i]:positions[i+1]]
                cut = copy.deepcopy(self.__sq_id)
                cut.seq = lseq
                cut.id = cut.id + f'_cut({positions[i]})-({positions[i+1]})'
                res.append(cut)
            return res
        
        else:
            print('[note] enter cut_pattern')
    
    # Function for when you want to prepare DNA sequence 
    # feature for ML applications
    def dnaseq_features(self,start=0,n_segs=101,seq_name=None):
        
        print(f'[note] cutting @{start} w/ {n_segs} segments')
        
        print(f"Input Sequence Length: {len(self.__sq_id.seq)}")
        remaind = len(self.__sq_id.seq)%n_segs
        if(remaind is not 0):
            last_id = len(self.__sq_id.seq) - remaind
        print(f"# Bases cut-off: {int(remaind)}")
        
        upd_seq = self.__sq_id.seq[start:last_id]
        
        print(f"Updated sequence length: {len(upd_seq)}")
        print(f"# Segments: {int(len(upd_seq)/n_segs)} created")
        if(seq_name is None):
            seq_name = 'seq'
            
        # store sequence subsets in a dictionary
        dic_seq = {}
        for i in range(0,3):
            a = int(i*n_segs) ; b = int(i*n_segs)+n_segs 
            identifier = f"{seq_name}_{a}:{b}"
            dic_seq[identifier] = upd_seq[a:b]
            
        lst_seq = dic_seq.values()
        index = list(dic_seq.keys())
        
        # One hot encode
        
        for ii,data in enumerate(lst_seq):
            
            abc = 'acgt'.upper()
            char_to_int = dict((c, i) for i, c in enumerate(abc))
            int_enc = [char_to_int[char] for char in data]
            
            ohe = []
            for value in int_enc:
                base = [0 for _ in range(len(abc))]
                base[value] = 1
                ohe.append(base)
            np_mat = np.array(ohe)
            np_mat = np.expand_dims(np_mat,axis=0)
            
            if(ii is not 0):
                matrix = np.concatenate([np_mat,matrix],axis=0)
            else:
                matrix = np_mat
                
        return matrix,index
    
    
'''

#######################################################

# Attributes:
- seq - sequence in string format (passed from SQ)
- seq_type - sequence type (passed from SQ)

# Methods:
- transcription() - Change DNA to RNA
- protein(mins=0) - Decode all amino acid chains for 
                    all 6 reading frames

# Decoding of instructions for making proteins from DNA
# transcription (DNA -> RNA)
# get_protein (RNA -> AA chains containing proteins)

#######################################################

'''

class Decode(SQ):
    
    # constructor
    def __init__(self,seq,seq_type):
        self.seq = seq
        self.seq_type = seq_type
        
    @staticmethod
    def dictmap(tid=None):

        tc = {
            "GCT":"A", "GCC":"A", "GCA":"A","GCG":"A",
            "TGT":"C", "TGC":"C","GAT":"D","GAC":"D",   
            "GAA":"E", "GAG":"E","TTT":"F","TTC":"F",   
            "GGT":"G", "GGC":"G","GGA":"G","GGG":"G",
            "CAT":"H", "CAC":"H","ATA":"I","ATT":"I", "ATC":"I",  
            "AAA":"K", "AAG":"K","TTA":"L","TTG":"L", "CTT":"L",  
            "CTC":"L", "CTA":"L","CTG":"L",
            "ATG":"M", # starting codon
            "AAT":"N", "AAC":"N","CCT":"P","CCC":"P", "CCA":"P", "CCG":"P",
            "CAA":"Q", "CAG":"Q","CGT":"R","CGC":"R", "CGA":"R",
            "CGG":"R", "AGA":"R","AGG":"R","TCT":"S", "TCC":"S", "TCA":"S",
            "TCG":"S", "AGT":"S","AGC":"S","ACT":"T", "ACC":"T", "ACA":"T", 
            "ACG":"T","GTT":"V", "GTC":"V","GTA":"V","GTG":"V","TGG":"W",
            "TAT":"Y", "TAC":"Y",
            "TAA":"_","TAG":"_","TGA":"_" # ending codon
            }

        if tid in tc: 
          return tc[tid]
        else: 
          return None
        
    ''' DNA Transcription '''
    # Convert DNA -> RNA chain
    
    def transcription(self):
        if (self.seq_type == "dna"):
            return SQ(self.seq.replace("T","U"), 
                      seq_type="rna")
        else:
            print('[note] seq_type != dna')
            return None
        
    # Translate sequence
    def trans(self,seq,p0=0):
        seq_aa = ""
        for pos in range(p0,len(seq)-2,3):
            cod = seq[pos:pos+3]
            seq_aa += self.dictmap(tid=cod)
        return seq_aa
    
    '''Get All Possible open reading frames (ORF)'''
    # store all possible collections of aa
    # groups in all 6 frames
    
    def frames(self):
        res = []
        for i in range(0,3):
            res.append(self.trans(self.seq,i))
        rc = self.reverse_comp()
        for i in range(0,3):
            res.append(self.trans(rc,i)) 
        return res
    
    ''' Computes all possible proteins in an aa sequence in RF '''
    # using the knowledge that it starts with M and ends with _, 
    # filter out rule breaking ORFs
    # aa_seq -> full converted amino acid sequence
    
    @staticmethod
    def all_rf(aa_seq):
        
        current_prot = []; proteins = []
        for aa in aa_seq:
            if(aa == "_"): # stopping gap
                if(current_prot):
                    for p in current_prot:
                        proteins.append(p)
                    current_prot = []
            else: # not stopping gap
                if(aa == "M"):      # starting amino acid            
                    current_prot.append("")
                for i in range(len(current_prot)):
                    current_prot[i] += aa
                    
        return proteins
    
    ''' Computes all possible putative proteins for all ORF '''
    # and sort them based on size
    
    def protein(self,mins=0):
        
        # order 
        def insert_prot_ord(prot, lst_prot):
            i = 0
            while(i < len(lst_prot) and len(prot)<len(lst_prot[i])):
                i += 1
            lst_prot.insert(i, prot)
            
        # get all ORF conversions
        rfs = self.frames();res = []
        for rf in rfs:
            # return only protein cases
            prots = self.all_rf(rf) 
            # additionally sort based on protein size
            for p in prots: 
                if(len(p) > mins): 
                    insert_prot_ord(p,res)
        return res
    
    
'''

#####################################################

# Counting
# freq - count bases in sequence
# count_purines - count purines & pyrimidines
# groupfreq - count grouped bases

#####################################################

'''

class Count:
    
    def __init__(self,seq,seq_type,id):
        self.seq = seq
        self.seq_type = seq_type
        self.id = id
        
    @staticmethod
    def dict_sum(dictlist):
        outdic = {}
        for d in dictlist:
            for k in d.keys():
                outdic[k] = 0
        for d in dictlist:
            for k in d.keys():
                outdic[k]+=d[k]
        return outdic
        
    # Get ABC
    def abc(self):
        if(self.seq_type=="dna"): 
          return "ACGT"
        elif(self.seq_type=="rna"):
          return "ACGU"
        elif (self.seq_type=="aa"): 
          return "ACDEFGHIKLMNPQRSTVWY"
        else: 
          return None
    
    # Mapping Dictionary
    @staticmethod
    def dictmap(map_id='iupac_amino',tid=None):

        # IUPAC Amino Acids
        if(map_id is 'iupac_amino'):
            tc   = {'A':'Alanine','C':'Cysteine','D':'Aspartic Acid',
                    'E':'Glutamic Acid','F':'Phenylalanine','G':'Glycine',
                    'H':'Histidine','I':'Isoleucine','L':'Lysine',
                    'M':'Methionine','N':'Asparagine','P':'Proline',
                    'Q':'Glutamine','R':'Arginine','S':'Serine',
                    'T':'Threonine','V':'Valine','W':'Tryptophan',
                    'Y':'Tryosine','_':'Gap'}

        # IUPAC nuceotides
        elif(map_id is 'iupac_nucleotide'):
            tc  = {'A':'Adenine','C':'Cytosine','G':'Guanine','T':'Thymine',
                   'U':'Uracil'}

        if tid in tc: 
          return tc[tid]
        else: 
          return None
        
    ''' Frequency of Alphabet '''
    # nucleotide or amino acid sequences
    
    def freq(self,compare:list=None,       # list of comparing sequences
                  show_id:str='perc',     # perc/count
                  fsize:list=[None,None],  # figure size
                  title:str=None,           # title
                  barwidth:float=0.25,            # bar width
                  barmode:str='group',        # bar display type
                  bargap:float=0.4,
                  bargroupgap:float=0.3):     # bar settings
        
        c1 = dict(Counter(self.seq))  # abc counter for s1
        abc = list(self.abc())
        count = Counter(abc)
        abc_c = dict(Counter({x:0 for x in count}))
        c_all1 = self.dict_sum([c1,abc_c])
            
        lst = []
        for i in c_all1.keys():
           if(self.seq_type == 'dna' or self.seq_type == 'rna'):
               lst.append(self.dictmap('iupac_nucleotide',i))
           elif(self.seq_type == 'aa'):
               lst.append(self.dictmap('iupac_amino',i))
                    
        perc = [round(x/len(self.seq),3)*100 for x in [*c_all1.values()]]
        if(show_id is 'perc'):
            show1 = lst; show2 = perc;
            xaxis_id = 'Character (%)'
        elif(show_id is 'count'):
            show1 = lst; show2 = [*c_all1.values()]
            xaxis_id = 'Count'
            
        if(self.id is not None):
            lname = self.id
        else:
            lname = 'main seq'
            
        fig = go.Figure(go.Bar(y=show1,x=show2,
                               orientation='h',name=lname))
        
        # Compare other sequences
        if(compare is not None):
            
            # Cycle through all SQ in list
            
            for ii,lseq in enumerate(compare):
                
                ii+=1
                c2 = dict(Counter(lseq.seq))  # abc counter for s2
                c_all2 = self.dict_sum([c2,abc_c])    
            
                lst2 = []
                for i in c_all2.keys():
                   if(self.seq_type == 'dna' or self.seq_type == 'rna'):
                       lst2.append(self.dictmap('iupac_nucleotide',i))
                   elif(self.seq_type == 'aa'):
                       lst2.append(self.dictmap('iupac_amino',i))
                
                perc = [round(x/len(lseq.seq),3)*100 for x in [*c_all2.values()]]
                if(show_id is 'perc'):
                    show1 = lst2;show2 = perc
                elif(show_id is 'count'):
                    show1 = lst2; show2 = [*c_all2.values()]
                
                if(lseq.id is not None):
                    lname = lseq.id
                else:
                    lname = f'seq{ii}'
                    
                fig.add_trace(go.Bar(y=show1,x=show2,
                                     orientation='h',
                                     name=lname))
                
        # Set title
        
        if(title is None):
            title = f'{self.seq_type} sequence content'
        else:
            title = title
            
        fig.update_layout(template='plotly_white',
                          barmode=barmode, # stack,group,overlay,relative
                          height=fsize[0],width=fsize[1],
                          bargroupgap=bargroupgap, 
                          bargap=bargap,
                          font={'size':11},
                          title=title)
                      
        fig.update_traces(width=barwidth,
#                           marker_color='rgb(158,202,225)', 
                          marker_line_color='#212121',
                          marker_line_width=1.0, opacity=0.6)
        fig.update_xaxes(nticks=20,title_text=xaxis_id);fig.show()
        
    # Count frequency of grouped nucleotides
    
    def gfreq(self,compare:list=None,
                   count_id='di',
                   fsize=[None,None],
                   barmode:str='group',
                   bargap:float=0.4,
                   bargroupgap:float = 0.3,
                   barwidth:float = 0.25,
                   xlim:int = 20):
        
        if(count_id is 'di'):
            lst_count_id = ['AA','AC','AG','AT',
                            'CA','CC','CG','CT',
                            'GA','GC','GG','GT',
                            'TA','TC','TG','TT']
            
        elif(count_id is 'tri'):
            lst_count_id = ['AAA','AAC','AAG','AAT','ACA','ACC','ACG',
                            'ACT','AGA','AGC','AGG','AGT','ATA','ATC',
                            'ATG','ATT','CAA','CAC','CAG','CAT','CCA',
                            'CCC','CCG','CCT','CGA','CGC','CGG','CGT',
                            'CTA','CTC','CTG','CTT','GAA','GAC','GAG',
                            'GAT','GCA','GCC','GCG','GCT','GGA','GGC',
                            'GGG','GGT','GTA','GTC','GTG','GTT','TAA',
                            'TAC','TAG','TAT','TCA','TCC','TCG','TCT',
                            'TGA','TGC','TGG','TGT','TTA','TTC','TTG',
                            'TTT']
            
        if(self.seq_type is 'dna'):
            
            lst_c = []
            for i in lst_count_id:
                lst_c.append(self.seq.count(i))
                
            df = pd.DataFrame(data=lst_c,
                              index=lst_count_id).T
            df.index = [self.id]
            
            if(compare is not None):
  
                ii=-1
                for seq in compare: # cycle through all SQ
                
                    ii+=1;lst_c = []
                    for jj in lst_count_id:
                        lst_c.append(compare[ii].seq.count(jj))
                        
                    ldf = pd.DataFrame(data=lst_c,
                                       index=lst_count_id).T
                    ldf.index = [seq.id]
                    df = pd.concat([df,ldf],axis=0)
            
            if(count_id is 'tri'):
                xaxis=dict(rangeslider=dict(visible=True))
            else:
                xaxis = None
                
            # Plot
            fig = px.bar(df.T,
                         x=df.T.index,
                         y=df.T.columns)
                      
            fig.update_layout(template='plotly_white',
                              barmode=barmode,
                              title=f'{count_id} nucleotide | count',
                              font={'size':11},
#                               margin=dict(l=20, r=20, t=80, b=20),
                              xaxis=xaxis,
                              bargroupgap=bargroupgap, 
                              bargap=bargap,
                              height=fsize[0],width=fsize[1]); 
                          
            fig.update_traces(width=barwidth,
    #                           marker_color='rgb(158,202,225)', 
                              marker_line_color='#212121',
                              marker_line_width=1.0, opacity=0.6)
            
            if(count_id is 'tri'):
                fig['layout']['xaxis'].update(title='', 
                                              range=[-1,xlim], 
                                              autorange=False) # range slider
            
            fig.show()
        else:
            print('[note] input must be dna type')
            
    # Return purines & pyrimidines (%)
    def purines(self):
        
        if (self.seq_type == "dna" or self.seq_type == "rna"):        
        
            dic_count = {}
            purines1 = self.seq.count("A") + self.seq.count("G")
            pyrimidines1 = self.seq.count("C") + self.seq.count("T")
            purines1 = purines1/len(self.seq)
            pyrimidines1 = pyrimidines1/len(self.seq)
            dic_count['purines'] = round(purines1,4)
            dic_count['pyrimidine'] = round(pyrimidines1,4)
            return dic_count
    
    # Return GC nucleotide (%)
    def gc(self):
        
        if (self.seq_type == "dna" or self.seq_type == "rna"):
            ii = 0
            for s in self.seq:
                if(s in "GC"):
                    ii += 1

            val = round(ii/len(self.seq),4)
            return val
        
''' 

#####################################################
    
# Finding Patterns in Sequence
# find - find index(ies) of particular pattern 
    
#####################################################

'''

class Pattern:
    
    def __init__(self,seq,seq_type):
        self.seq = seq
        self.seq_type = seq_type
    
    # Prosite Pattern Lines
    # - Standard IUPAC amino acid used to as bases in pattern, separated by -
    # - x -> any amino acid acceptable
    # - [] -> ambiguity represented by list, any aa in that list acceptable
    # - {} -> ambiguity represented by list, any aa other than in {} accepted
    # - repetition of pattern element shown below:
    #   x(3) -> x-x-x, x(2,4) -> to x-x or x-x-x or x-x-x-x
    
    @staticmethod
    def prosite_process(rex):
        # adjust prosite to RE format
        rex = rex.replace("(","{")
        rex = rex.replace(")","}")
        rex = rex.replace("x",".")
        rex = rex.replace("-","")
        return rex
    
    ''' Find Substring Pattern in Sequence '''
    
    def find(self,pattern,          # subsequence pattern
                  find_id='first',  # first,all,overlap
                  search_id=None,   # search pattern option (prosite)
                  verbose=True):    # Verbal Output
        
        ''' [1] Find First Match Only '''
        
        if(find_id is 'first'):
            
            if(search_id is 'prosite'):
                pattern = self.prosite_process(pattern)
                
            # General search as well
            re_search = search(pattern,self.seq)
            if (re_search != None):
                if(verbose):
                    print(f"showing first for {pattern}")
                result = re_search.span()[0]
                return result
            else:
                if(verbose):
                    print(f'no matches for {pattern} found')
                
        elif(find_id == 'all'):
            
            if(search_id is 'prosite'):
                pattern = self.prosite_process(pattern)
                
            re_search = finditer(pattern,self.seq)
            result = []
            for x in re_search:
                result.append(x.span()[0])
                
            if(len(result) is not 0):
                if(verbose):
                    print(f"found {len(result)} matches")
                return result
            else:
                if(verbose):
                    print(f'no matches for {pattern} found')
                
        elif(find_id == 'overlap'):
            
            if(search_id is 'prosite'):
                pattern = self.prosite_process(pattern)
            mos = finditer("(?="+pattern+")",self.seq)
            result = []
            for x in mos:
                result.append(x.span()[0])
                
            if(len(result) is not 0):
                if(verbose):
                    print(f"found {len(result)} matches")
                return result
            else:
                if(verbose):
                    print(f'no matches for {pattern} found')
        else:
            print('first,all,overlap options')
            
            
# s1 = SQ("ATAGCATTACAG",id='new')
# s2 = SQ("AACC",id='news')
# s3 = SQ("ATGACTAAGAC",id='seq3')
# s4 = SQ("ATGAAGACACACAC",id='seq4')

# # # Visualise sequence
# # view_example = SQ("ATAGCATTACAGATAGCATTATAGCATTACAGATGCATTG",id='new')
# # view_example.view(split_id=5)

# # Defining a detailed sequence
# s1 = SQ(seq="ATAGCA",
#         seq_type='dna',
#         id='db|12321',
#         name='short sequence',
#         description='we will use this sequence for analysis')

# # print sequence info
# s1.info()

# s1.validate()

# # convert sequence to pandas series
# series = s1.to_series()
# display(series)