# Class to access BioPython's blastwww & plot alignments 
class BLASTwww:
    
    def __init__(self,seq=None,database='pdb',verbose=False,
                 show_top=10,hit_id=0,read_xml=None,name=None):
        self.seq = seq 
        self.database = database
        self.verbose = verbose
        self.show_top = show_top
        self.hit_id = hit_id
        self.read_xml = read_xml # if present -> read xml
        self.name = name  # search identifier (used for save)
        
    def save(self):
        # save search xml
        sf = open(f"/kaggle/working/blast_{self.name}.xml", "w")
        sf.write(self.result_handle.read()) 
        self.result_handle.close()
        sf.close() 

        # save metadata
        self.df.to_csv(f'/kaggle/working/csv_{self.name}.csv')
        
    # Find amino acid sequence in NCBI databases
    def blast_aa(self):

        # if not seq format
        if(type(self.seq) is str):
            self.seq = Seq(self.seq)
        
        if(self.read_xml is not None):
            self.result_handle = open(f"/kaggle/working/blast_{self.name}.xml")
            blast_qresult = SearchIO.read(self.result_handle,"blast-xml")
        else:
            self.result_handle = NCBIWWW.qblast("blastp",self.database,self.seq)
            blast_qresult = SearchIO.read(self.result_handle,"blast-xml")
 
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

        self.df = pd.DataFrame(dic_data)

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
            
    @staticmethod
    def get_colors(seqs):
        """make colors for bases in sequence"""
        text = [i for s in list(seqs) for i in s]

        # DNA
    #     clrs =  {'A':'red','T':'green','G':'orange','C':'blue','-':'white'}
        # IUPAC aa
        clrs   = {'A':'#3386FF','C':'#3386FF','D':'#B842B2','E':'#B842B2',
                'F':'#3386FF','G':'#FF5733','H':'#37ADBB','I':'#3386FF',
                'L':'#3386FF','M':'#3386FF','N':'#24CE5D','P':'#E3E710',
                'Q':'#24CE5D','R':'#D3385E','S':'#24CE5D','T':'#24CE5D',
                'V':'#3386FF','W':'#3386FF','Y':'#37ADBB','_':'white',
                'K':'#D3385E'}

        colors = [clrs[i] for i in text]
        return colors
            
    # View BLAST alignment using Bokeh
    def view(self,aln, fontsize="9pt", plot_width=800):

        #make sequence and id lists from the aln object
        seqs = [rec.seq for rec in (aln)]
        ids = [rec.id for rec in aln]    
        text = [i for s in list(seqs) for i in s]
        colors = self.get_colors(seqs)    
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
