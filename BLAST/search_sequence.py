# Class to access BioPython's blastwww & plot alignments 
class BLASTwww:
    
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
        self.query_df = None   # dataframe that stored query data
        
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
        
    ''' Search for sequence in database '''
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
            print(f'\nShowing BLAST query result #{hit_id}')
            print('------------------------------------------------------------')
            seqid = blast_qresult[hit_id]
            details = seqid[hit_id]

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
    # After fetch() method, we can visualise the results in df 
    
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
    

''' Example Usage '''

# We can utlise the SQ class to get amino acid chain
str_seq = 'ATTAAAGGTTTATACCTTC'      # Sequence in string format
seq = SQ(str_seq)                    # define SQ class sequence
lst_aa = seq.get_protein(min_size=1) # extract amino acid chains from via translation
lst_aa                               # visualise all created amino acid chains

blast_query = BLASTwww(lst_aa[0])    # Instantiate class
blast_query.fetch()                  # BLAST query
blast_query.view_query()             # view BLAST query results (results stored in dataframe)\

# Visualise quary
plot = blast_query.view_alignment(0) # view first alignment from results dataframe
pn.pane.Bokeh(plot)
