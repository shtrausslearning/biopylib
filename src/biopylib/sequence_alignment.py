import numpy as np
import pandas as pd
import panel as pn
import panel.widgets as pnw
pn.extension()
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Plot, Grid, Range1d
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot
from bokeh.transform import dodge
import bokeh
import itertools
from itertools import islice
from biopylib.sequence import SQ
from Bio import SeqIO, SearchIO
from Bio.Blast import NCBIWWW 
from Bio.Seq import Seq

'''

####################################################################

# Substitution Matrix Class
# Class for storing the substitution matrix, reading & making

# Attributes:
- abc
- subm 

# Methods:
- __getitem__ - obtain the score from the substitution matrix
- score_pair - obtain the score from the substitution matrix
- head - show the first few entries in the substitution matrix 
- read - read a substitution matrix file
- make - create a simple hit/miss substitution matrix

####################################################################

'''

class SUBM:
    
    # Constructor & class operation
    def __init__(self):
        self.abc = ""        # characters
        self.subm = {}       # Substitution matrix mapping
        
    # get score pair class operation
    def __getitem__(self, ij):
        i,j = ij
        return self.score_pair(i, j)
    
    # get specific score of a pair
    def score_pair(self,c1,c2):
        if(c1 not in self.abc or c2 not in self.abc):
            return None
        return self.subm[c1+c2]
    
    # show first n number of pairs 
    def head(self,n=5):
        return list(islice(self.subm.items(),n))
    
    # Read Substitution Matrix from file (.mat)
    # Format includes separation by tabs
    
    def read(self,name):
        
        # Read the code
        file = open(name, "r")
        abc = file.readline().split('\t')
        
        self.abc = ""
        for i in range(0, len(abc)): 
            self.abc += abc[i][0]
        for i in range(0,len(abc)):
            line = file.readline();
            abc = line.split('\t');
            for j in range(0, len(abc)):
                k = self.abc[i]+self.abc[j]
                self.subm[k] = int(abc[j])
        file.close()
        
    # Create Substitution Matrix
    # requires match & mismatch score & alphabet
    # A constant match, mismatch score is used
    
    def make(self,code,match,missmatch):
        self.abc = code
        for values in itertools.product(code,code):
            c1 = values[0]; c2 = values[1]
            if (c1 == c2):
                self.subm[c1+c2] = match
            else:
                self.subm[c1+c2] = missmatch
                

####################################################################

# Alignment Visualisation using Bokeh
# @https://www.kaggle.com/nwheeler443/understanding-phylogenetic-trees

####################################################################
class aliView:
    
    def __init__(self,aln=None,colcod=None):
        if(aln is not None):
            self.aln = aln
        else:
            print('alignment not defined')
        if(colcod is not None):
            self.colcod = colcod
        else:
            self.colcod = 'general'
                
    @staticmethod
    def get_colors(seqs,ids):
        """make colors for bases in sequence"""
        text = [i for s in list(seqs) for i in s]

        # DNA
        # clrs =  {'A':'red','T':'green','G':'orange','C':'blue','-':'white'}
        
        if(ids is 'general'):
            # General colourcoding (n/aa chain)
            clrs   = {'A':'#2C679E','C':'#2C679E','D':'#B842B2','E':'#B842B2',
                    'F':'#2C679E','G':'#FF5733','H':'#37ADBB','I':'#2C679E',
                    'L':'#2C679E','M':'#2C679E','N':'#24CE5D','P':'#E3E710',
                    'Q':'#24CE5D','R':'#D3385E','S':'#24CE5D','T':'#24CE5D',
                    'V':'#2C679E','W':'#2C679E','Y':'#37ADBB','_':'white',
                    'K':'#D3385E','-':'white'}
        
        elif(ids is 'charge_aa'):
            # +/-ve charged colourcoding for aa chain only
            clrs   = {'A':'#8081f7','C':'#8081f7','D':'#fe8081','E':'#fe8081',
                    'F':'#8081f7','G':'#8081f7','H':'#b3ff80','I':'#8081f7',
                    'L':'#8081f7','M':'#8081f7','N':'white','P':'#8081f7',
                    'Q':'white','R':'#b3ff80','S':'white','T':'#white',
                    'V':'#8081f7','W':'#8081f7','Y':'white','_':'white',
                    'K':'#b3ff80','-':'white'}

        colors = [clrs[i] for i in text]
        return colors
            
    # View BLAST alignment using Bokeh
    def view(self, fontsize="9pt", plot_width=800):

        #make sequence and id lists from the aln object
        seqs = [rec.seq for rec in (self.aln)] # list of sequences
        ids = [rec.id for rec in self.aln]    
        text = [i for s in list(seqs) for i in s]
        colors = self.get_colors(seqs,self.colcod)    
        N = len(seqs[0]); S = len(seqs)    

        x = np.arange(0.5,N+0.5); y = np.arange(0,S,1)
        #creates a 2D grid of coords from the 1D arrays
        xx, yy = np.meshgrid(x, y)
        #flattens the arrays
        gx = xx.ravel(); gy = yy.flatten()
        #use recty for rect coords with an offset
        recty = gy+0.5; h= 1/S
        #now we can create the ColumnDataSource with all the arrays
        source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, text=text, colors=colors))
        plot_height = len(seqs)*15+25
        x_range = Range1d(0,N+1, bounds='auto')
        if(N>100):
            viewlen=50
        else:
            viewlen=N
        #view_range is for the close up view
        view_range = (0,viewlen)
        tools="xpan, xwheel_zoom, reset, save" 

        #entire sequence view (no text, with zoom)
        p = figure(title=None, plot_width= plot_width, plot_height=50,
                   x_range=x_range, y_range=(0,S), tools=tools,
                   min_border=0, toolbar_location='below')
        rects = Rect(x="x", y="recty",  width=1, height=0.5, fill_color="colors",
                     line_color=None, fill_alpha=0.7)
        p.add_glyph(source, rects)
        p.yaxis.visible = False
        p.grid.visible = False  

        #sequence text view with ability to scroll along x axis
        p1 = figure(title=None, plot_width=plot_width, plot_height=plot_height,
                    x_range=view_range, y_range=ids, tools="xpan,reset",
                    min_border=0, toolbar_location='below')#, lod_factor=1)          
        glyph = Text(x="x", y="y", text="text", text_align='center',text_color="black",
                    text_font_size=fontsize,text_font_style='bold')
        rects = Rect(x="x", y="recty",  width=1.0, height=1, fill_color="colors",
                    line_color=None, fill_alpha=0.7)
    
        p1.add_glyph(source, rects)
        p1.add_glyph(source, glyph)

        p1.grid.visible = False
        p1.xaxis.major_label_text_font_style = "bold"
        p1.yaxis.minor_tick_line_width = 0
        p1.yaxis.major_tick_line_width = 0
        
        return gridplot([[p],[p1]],toolbar_location='below') 
    
####################################################################
    
''' Alignment Class ''' 
# Class Used for Alignment Storage & MSA consensus
    
####################################################################

class ALI(aliView):
    
    ''' Constructor '''
    def __init__(self, lst_seqs, al_type = "aa"):
        self.lst_seqs = lst_seqs  # stored alignment sequences
        self.al_type = al_type    # alignment type
        self.nseq = len(self.lst_seqs) # number of alignment
        
        it = iter(self.lst_seqs)
        the_len = len(next(it))
        if not all(len(l) == the_len for l in it):
             raise ValueError('Sequences must all be the same length')
            
        # Keep SQ storage
        self.lst_SQ = []; ii=-1
        for seq in lst_seqs:
            ii+=1; self.lst_SQ.append(SQ(seq=seq,id=f'sequence{ii}'))
            
    ''' 2. Print Alignment '''
    def __str__(self):
        pprint = ""
        print(f'Alignment with {len(self.lst_seqs)} rows and {len(self)} columns.')
        for seq in self.lst_seqs:
            pprint += seq + "\n"
        return pprint
    
    # len(class); number of characters
    def __len__(self): 
        return len(self.lst_seqs[0]) # all should be identical
    
    ''' 3. Get nth row in alignment '''
    def __getitem__(self, n):
        if( type(n) is tuple and len(n) == 2):
            i,j = n
            return self.lst_seqs[i][j]
        elif(type(n) is int): 
            return self.lst_seqs[n]
        return None
    
    ''' 4. Get nth column in alignment '''
    
    def col(self,n):
        res = []
        for k in range(len(self.lst_seqs)):
            res.append(self.lst_seqs[k][n])
        return res
    
    # Find a consensus between multiple sequences
    def consensus(self):
        
        cons = ""
        for i in range(len(self)):
            cont = {}
            for k in range(self.nseq):
                c = self.lst_seqs[k][i]
                if(c in cont):
                    cont[c] = cont[c] + 1
                else: 
                    cont[c] = 1
                    
            maxs = 0; cmax = None
            for ii in cont.keys():
                if(ii != "-" and cont[ii] > maxs): 
                    maxs = cont[ii]
                    cmax = ii
            cons+= cmax
        return cons
    
    # View alignment (seconday function)
    def view(self):
        plot = aliView(aln=self.lst_SQ,
                       colcod=None).view()
        display(pn.pane.Bokeh(plot))
        
####################################################################
        
''' Pairwise Sequence Alignment Custom Class '''
# Needleman-Wunsch & Smith Waterman DP approaches
        
####################################################################
        
class PSA(SQ,aliView):
    
        # PSA Constructor
        def __init__(self,seqs=None, # List of Sequences
                                            subm=None, # Substitution Matrix Class
                                            g=None,  # Constand Gap Penalty
                                            colcod=None):
            
                # If Sequences are initially given
                if(seqs is not None):
                        self.seq1 = seqs[0]   # sequence class SQ
                        self.seq2 = seqs[1]   # sequence class SQ
                        self.l1 = len(self.seq1)  # length of seq1
                        self.l2 = len(self.seq2)  # length of seq2
                        if(len(seqs) != 2):
                                print('Require 2 sequences')
                            
                # Gap Penalty & Substitution Matrix
                self.g = g         # constant gap penalty model
                self.subm = subm   # substitution matrix class
                self.score = 0    # Either Max Score @sW or final value @nW
                if(colcod is not None):
                        self.colcod = colcod
                else:
                        self.colcod = 'general'
                    
        @staticmethod
        def maxidx(v1,v2,v3):
                if(v1 > v2):
                        if(v1 > v3): return 1
                        else: return 3
                else:
                        if(v2 > v3): return 2
                        else: return 3
                    
        def score_pos(self,c1,c2):
                if(c1 != '-' or c2 != '-'):
                        return self.subm[c1,c2]
                elif(c1 == "-" or c2=="-"):
                        return self.g
            
        ''' Initialisation of Scoring & Traceback Matrices '''
        # Fill Top left corner values
        # Traceback Movement interpretation/mapping
        # 0 - done, 1 - diagonal, 2 - up, 3 - left
            
        def initial_condition(self,ids='nW'):
            
                # Initialisation of matrices; top left -> 0
                self.SM = [[0]]; self.TM = [[0]]
            
                # Fill Matrices SM & TM
                if(ids is 'nW'):
                        # Fill Rows of first column
                        for i in range(1,self.l1+1):
                                self.SM.append([self.g * i]); 
                                self.TM.append([2])
                        # Fill the First Row
                        for j in range(1,self.l2+1):
                                self.SM[0].append(self.g * j); 
                                self.TM[0].append(3)
                elif(ids is 'sW'):
                        for i in range(1,self.l1+1):
                                self.SM.append([0]); 
                                self.TM.append([0])
                        for j in range(1,self.l2+1):
                                self.SM[0].append(0); 
                                self.TM[0].append(0)
                            
        ''' (1) Sequence Alignment '''
                            
        # Needleman-Wunsch Global Alignment
        def nW(self):
            
                # Initialise Matrices (fill first row/column)
                self.initial_condition(ids='nW')
            
                # Fill in matrix
                for values in itertools.product(range(0,self.l1),range(self.l2)):
                        i = values[0]; j = values[1]
                        get_score = self.score_pos(self.seq1[i],self.seq2[j])
                        s1 = self.SM[i][j] + get_score
                        s2 = self.SM[i][j+1] + self.g
                        s3 = self.SM[i+1][j] + self.g
                        # Find Maximum
                        self.SM[i+1].append(max(s1,s2,s3))
                        self.TM[i+1].append(self.maxidx(s1,s2,s3))
                    
                # Return alignment score
                self.score = self.SM[self.l1][self.l2]
            
        # Smith-Waterman Local Alignment (get max score)
        def sW(self):
            
                # Initialise Matrices (fill first row/column)
                self.initial_condition(ids='sW')
            
                # Fill in Score Matrix & Traceback Matrix
                max_score = 0
                for values in itertools.product(range(self.l1),range(self.l2)):
                        i = values[0]; j = values[1]
                        get_score = self.score_pos(self.seq1[i],self.seq2[j]) 
                    
                        s1 = self.SM[i][j] + get_score
                        s2 = self.SM[i][j+1] + self.g
                        s3 = self.SM[i+1][j] + self.g
                        b = max(s1, s2, s3)
                        if b <= 0:
                                self.SM[i+1].append(0)
                                self.TM[i+1].append(0)
                        else:
                                self.SM[i+1].append(b)
                                self.TM[i+1].append(self.maxidx(s1,s2,s3))
                                if(b > max_score): 
                                        max_score = b
                                    
                # Store Final Value
                self.score = max_score
            
        ''' (2) Recover Alignment '''
            
        # Global Alignment (nW)
        def realign(self,ids='nW'):
            
                ali_seq = {'sq1':'','sq2':''}
            
                # Get Max Score Case
                if(ids is 'nW'):
                        i = self.l1; j = self.l2
                elif(ids is 'sW'):
                    
                        max_score = 0
                        lst_max = [0,0]
                        for i in range(1,len(self.SM)):
                                for j in range(1, len(self.SM[i])):
                                        # get maximum score location
                                        if(self.SM[i][j] > max_score):
                                                max_score = self.SM[i][j]
                                                lst_max[0] = i; lst_max[1] = j
                                            
                        i = lst_max[0]  # maximum in row 
                        j = lst_max[1]  # maximum in column
                    
                while(i>0 or j>0):       # for all except done
                        if(self.TM[i][j]==1): # diagonal
                                ali_seq['sq1'] = self.seq1[i-1] + ali_seq['sq1']
                                ali_seq['sq2'] = self.seq2[j-1] + ali_seq['sq2']
                                i -= 1; j -= 1
                        elif(self.TM[i][j] == 2): # up
                                ali_seq['sq1'] = self.seq1[i-1] + ali_seq['sq1']
                                ali_seq['sq2'] = '-' + ali_seq['sq2']
                                i -= 1
                        elif(self.TM[i][j] == 3): # left
                                ali_seq['sq1'] = '-' + ali_seq['sq1']
                                ali_seq['sq2'] = self.seq2[j-1] + ali_seq['sq2']
                                j -= 1
                        else:
                                break
            
                self.pali_str = [ali_seq['sq1'], # alignment sequences
                                                    ali_seq['sq2']] # string format
                self.pali = ALI(self.pali_str,self.seq1.seq_type)
            
                # return alignment sequence (twin sequence class)
                return self.pali
    
        # Visualise Alignment; only after realignment
    
        def view(self):
            
                # if id not given, use default names
                self.lst_SQ = []
                if(self.seq1.id is None):
                        self.seq1.id = 'Main Sequence'
                if(self.seq2.id is None):
                        self.seq2.id = 'Comparing Sequence'
                    
                self.lst_SQ.append(SQ(seq=self.pali_str[0],id=self.seq1.id))
                self.lst_SQ.append(SQ(seq=self.pali_str[1],id=self.seq2.id))
                plot = aliView(aln=self.lst_SQ,
                                                colcod=self.colcod).view()
                display(pn.pane.Bokeh(plot))
                
'''

####################################################################
            
Progressive Multiple Sequence Alignmnent
# Using consensus between two sequences we can iteratively add 
# and align a list of sequences 
            
#####################################################################

'''
        
class MSA(PSA):
    
    def __init__(self, seqs, subm=None,g=None,colcod=None):
        
        tlst = []
        self.seqs = seqs # list of sequences SQ objects
        self.subm = subm # substitution matrix
        self.g = g       # gap penalty
        
        if(colcod is not None):
            self.colcod = colcod
        else:
            self.colcod = 'general'
        
    # total number of sequences in input
    def num_seqs(self):
        return len(self.seqs)
    
    # helper function for align_consensus
    def add_ali(self, alignment, seq):
        
        res = []; 
        upd_len = len(alignment.lst_seqs)+1
        for i in range(upd_len):
            res.append("")
            
        # Using consensus of first two, create a sequence
        cons = SQ(alignment.consensus(),alignment.al_type)
        # PSA global alignment of consensus seq & next seq in list
        lpsa = PSA(seqs=[cons,seq],subm=self.subm,g=self.g)
        lpsa.nW()
        align2 = lpsa.realign(ids='nW')
        
        orig = 0
        for i in range(len(align2)):
            if align2[0,i]== '-':
                for k in range(len(alignment.lst_seqs)):
                    res[k] += "-"
            else:
                for k in range(len(alignment.lst_seqs)):
                    res[k] += alignment[k,orig]
                orig+=1
        res[len(alignment.lst_seqs)] = align2.lst_seqs[1]
        return ALI(res, alignment.al_type)
    
    # Find iterative concensus for all sequences
    def align_consensus(self):
        
        # First PSA global alignment
        cons = PSA(seqs=[self.seqs[0], self.seqs[1]],
                   subm=self.subm,g=self.g)
        cons.nW()
        ali = cons.realign(ids='nW')
        
        # Subsequent sequences are added the first two
        for i in range(2, len(self.seqs)):
            ali = self.add_ali(ali, self.seqs[i])
            
        self.lst_SQ = []; ii=-1
        for i in range(len(self.seqs)):
            ii+=1; self.lst_SQ.append(SQ(seq=ali[ii],id=self.seqs[ii].id))
            
        return ali # return final alignmnent
    
    # Visualise alignmnent using Bokeh
    def view(self):
        
        plot = aliView(aln=self.lst_SQ,
                        colcod=self.colcod).view()
        display(pn.pane.Bokeh(plot))