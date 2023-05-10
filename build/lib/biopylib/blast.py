from Bio import SeqIO, SearchIO
from Bio.Blast import NCBIWWW 
from Bio.Blast import NCBIXML 
from Bio.Seq import Seq
from biopylib.sequence import SQ
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
import re
import numpy as np

'''

Class to access BioPython's blastwww & plot alignments 

'''

class BLAST:
    
    # Constructor
    def __init__(self,seq=None,        # input sequence (Seq/string)
                      program='blastp', # BLAST program (blastn/blastp)
                      database='pdb',  # BLAST query database
                      verbose=False,   # write outputs
                      show_top=10,     # maximum number of hits
                      hit_id=0,        # show query match id #
                      read_xml=False,   # read XML BLAST query over new query
                      name=None):      # query identifier
        
        self.seq = seq 
        self.database = database
        self.verbose = verbose
        self.show_top = show_top
        self.hit_id = hit_id
        self.read_xml = read_xml # if present -> read xml
        self.name = name  # search identifier (used for save)
        self.program = program # blast search program
        self.aquery_df = None   # dataframe that stored query data
        self.gquery_df = None  # Query Output using 
        
        # Colourcoding dictionary
        self.dict_cc = {'A':'#3386FF','C':'#3386FF','D':'#B842B2','E':'#B842B2',
                        'F':'#3386FF','G':'#FF5733','H':'#37ADBB','I':'#3386FF',
                        'L':'#3386FF','M':'#3386FF','N':'#24CE5D','P':'#E3E710',
                        'Q':'#24CE5D','R':'#D3385E','S':'#24CE5D','T':'#24CE5D',
                        'V':'#3386FF','W':'#3386FF','Y':'#37ADBB','_':'white',
                        'K':'#D3385E'}
        
    ''' Save BLAST search in CSV format '''
    # read xml & convert to pandas dataframe
    
    def save(self):
        
        # save search xml
        sf = open(f"/kaggle/working/blast_{self.name}.xml", "w")
        sf.write(self.__result_handle.read()) 
        self.__result_handle.close()
        sf.close() 

        # save metadata
        self.df.to_csv(f'/kaggle/working/csv_{self.name}.csv')
        
    ''' Search for aminoacid chain in database '''
    # Find amino acid sequence in NCBI databases
    # by default, program -> blastp (protein/aa chain search)
    # To prevent multiple fetches, read xml is active & reads xml if fetched
    # more than once
    
    def fetch(self):

        # Input in BioPython Sequence Format
        if(type(self.seq) is str):
            self.seq = Seq(self.seq)
        
        if(self.read_xml == True):
            self.__result_handle = open(f"/kaggle/working/blast_{self.name}.xml")
            blast_qresult = SearchIO.read(self.__result_handle,"blast-xml")
        else:
            self.__result_handle = NCBIWWW.qblast(self.program,
                                                  self.database,
                                                  self.seq)
            
            blast_qresult = SearchIO.read(self.__result_handle,"blast-xml")
            self.read_xml = True # blast query has been saved
 
        ii=-1
        lst_id = []; lst_descr = []; lst_bitscore = []
        lst_ali = []; lst_eval = []
        
        for f in blast_qresult: 
            ii+=1
            seqid = blast_qresult[ii]
            details = seqid[0]

            lst_id.append(seqid.id)
            lst_descr.append(seqid.description)
            lst_eval.append(details.evalue)
            lst_bitscore.append(details.bitscore)
            lst_ali.append(details.aln)

        dic_data = {'id':lst_id,'description':lst_descr,
                    'evalue':lst_eval,'bitscore':lst_bitscore,
                    'alignment':lst_ali}

        self.__query_df = pd.DataFrame(dic_data)

        if(self.verbose):

            print(blast_qresult[0:self.show_top])
            #fetch the id, description, evalue, bitscore & alignment
            print(f'\nShowing BLAST query result #{self.hit_id}')
            print('------------------------------------------------------------')
            seqid = blast_qresult[self.hit_id]
            details = seqid[self.hit_id]

            print(f'Sequence ID: {seqid.id}')
            print(f'Description: {seqid.description}')
            print(f'e-value: {details.evalue}')
            print(f'Bit Score: {details.bitscore}')

            print('\nShowing Alignment:')
            print('------------------------------------------------------------')
            print(f"Alignment:\n{details.aln}")
            
    def __get_colors(self,seqs):
        """make colors for bases in sequence"""
        text = [i for s in list(seqs) for i in s]
        colors = [self.dict_cc[i] for i in text]
        return colors
    
    ''' Visualise BLAST query '''
    # After fetch, view result dataframe
    
    def view_query(self):
        
        if(self.__query_df is not None):
            return self.__query_df
            
    ''' View Alignment '''
    # View BLAST alignment using Bokeh
    
    def view_alignment(self,aln_id,
                  fontsize="9pt",
                  plot_width=800):
        
        # Choose Alignment 
        if(self.__query_df is not None):
            aln_id = self.__query_df.loc[aln_id,'alignment']
            
        #make sequence and id lists from the aln object
        seqs = [rec.seq for rec in (aln_id)]
        ids = [rec.id for rec in aln_id]    
        text = [i for s in list(seqs) for i in s]
        colors = self.__get_colors(seqs)    
        N = len(seqs[0])
        S = len(seqs)    

        x = np.arange(0.5,N+0.5)
        y = np.arange(0,S,1)
        #creates a 2D grid of coords from the 1D arrays
        xx, yy = np.meshgrid(x, y)
        #flattens the arrays
        gx = xx.ravel()
        gy = yy.flatten()
        #use recty for rect coords with an offset
        recty = gy+0.5
        h= 1/S
        #now we can create the ColumnDataSource with all the arrays
        source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, text=text, colors=colors))
        plot_height = len(seqs)*15+25
        x_range = Range1d(0,N+1, bounds='auto')
        if N>100:
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
        rects = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors",
                     line_color=None, fill_alpha=0.6)
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
                    line_color=None, fill_alpha=0.4)
        p1.add_glyph(source, glyph)
        p1.add_glyph(source, rects)

        p1.grid.visible = False
        p1.xaxis.major_label_text_font_style = "bold"
        p1.yaxis.minor_tick_line_width = 0
        p1.yaxis.major_tick_line_width = 0
        
        return gridplot([[p],[p1]],toolbar_location='below') 
    
    
    # Search NCBI BLAST Query
    def gfetch(self,verbose=False,n_top=5):
        
        # if input is a string
        if(type(seq) is str):
            seq = Seq(seq)
        
        res_handle = NCBIWWW.qblast(self.program,self.database,self.seq)
        sf = open("/kaggle/working/temp_blast.xml", "w")
        sf.write(res_handle.read()) 
        sf.close() 
        res_handle.close()
        res_handle = open("/kaggle/working/temp_blast.xml")
        record = NCBIXML.read(res_handle)

        print (f"SEARCHING DATABASE: {record.database}")
        print (f"SUBSTITUTION MATRIX USED in PSA: {record.matrix}")
        print (f'GAP PENALTY IN PSA: {record.gap_penalties}')

        print (f"TOT # OF HITS: {len(record.alignments)}")
        
        if(len(record.alignments) != 0):
            
            first_alignment = record.alignments[0]
            hsp = first_alignment.hsps[0]

            if(verbose):
                print(f"\nFIRST HIT #1:")
                print('--------------------------------------------------------')
                print (f"ID: {first_alignment.hit_id}")
                print (f"ACCESSION: {first_alignment.accession}")
                print (f"DEFINITION: {first_alignment.hit_def}")
                print (f"ALI LENGTH: {first_alignment.length}")
                print (f"HSP NUMBER: {len(first_alignment.hsps)}")

                print (f"e-VALUE: ", hsp.expect)
                print (f"SCORE: {hsp.score}")
                print (f"LENGTH: {hsp.align_length}")
                print (f"IDENTITIES: {hsp.identities}")
                print ("ALI OF HSPS:")
                print (f'QUERY: {hsp.query}')
                print (f'MATCH: {hsp.match}')
                print (f'SUBJE: {hsp.sbjct}')
                print('--------------------------------------------------------\n')
                print("TOP 5 ALIGNMENTS:")
                print('--------------------------------------------------------\n')

            lst_hit = []; lst_access = []; lst_alilen = []
            lst_hsps = []; lst_eval = []; lst_score = []
            lst_len = []; lst_ident = []; lst_def = []
            lst_quer = []; lst_mat = []; lst_sub = []

            # store data for alignment
            for i in range(n_top):

                loc_al = record.alignments[i]
                hsp = loc_al.hsps[0] # store only first (can be multi)

                lst_hit.append(loc_al.hit_id)
                lst_access.append(loc_al.accession)
                lst_def.append(loc_al.hit_def)
                lst_alilen.append(loc_al.length)
                lst_hsps.append(len(loc_al.hsps))

                lst_eval.append(hsp.expect)
                lst_score.append(hsp.score)
                lst_len.append(hsp.align_length)
                lst_ident.append(hsp.identities)

                lst_quer.append(hsp.query)
                lst_mat.append(hsp.match)
                lst_sub.append(hsp.sbjct)

                if(verbose):
                    print (f"ACCESSION: {loc_al.accession}")
                    print (f"DEFINITION: {loc_al.hit_def}")
                    for hsp in first_alignment.hsps:
                        print (f"e-VALUE: {hsp.expect}")

            dic_data = \
            {"ID":lst_hit,"ACCESSION":lst_access,"DEFINITION":lst_def,
            "ALI LENGTH":lst_alilen,"HSP NUMER":lst_hsps,"e-VALUE":lst_eval,
            "SCORE":lst_score,"LENGTH":lst_len,"IDENTITIES":lst_ident,
            "SEQ-QUERY":lst_quer,"SEQ-MATCH":lst_mat,"SEQ-SUBJECT":lst_sub,
            }
            
            # Store General Query DataFrame
            self.gquery_df = pd.DataFrame(dic_data)

            lst_org = []
            for i in range(n_top):
                ali = record.alignments[i]
                defin = ali.hit_def
                x = re.search("\[(.*?)\]", defin).group(1)
                lst_org.append(x)

                if(verbose):
                    print ("ORGANISMS:")
                    print('--------------------------------------------------------\n')
                    for s in lst_org: 
                        print(s)

           