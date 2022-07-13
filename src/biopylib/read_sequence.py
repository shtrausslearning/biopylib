from biopylib.sequence import SQ
from re import sub, search
import numpy as np

# NCBI identifiers
identifiers_dic = {'lcl':'local(nodb)','bbs':'GenInfo backbone seqid',
                   'bbm':'GenInfo backbone moltype','gim':'GenInfo import ID',
                   'gb':'GenBank','emb':'EMBL','pir':'PIR','sp':'SWISS-PROT',
                   'pat':'patent','pgp':'pre-grant patent','ref':'RefSeq',
                   'gnl':'general database reference','prf':'PRF','pdb':'PDB',
                   'gi':'GenInfo integrated database','dbj':'DDBJ'}

# FASTA format names
dict_FASTA = {'fa':'generic','fasta':'generic','fna':'nucleic acid',
              'ffn':'nucleotide of gene regions','faa':'amino acid',
              'frn':'non-coding RNA'}

# GENBANK format names
dict_GB = {'gb':'genbank'}

# Class to read different files and store info only

class readSeq(SQ):
    
    # Constructor
    
    def __init__(self,name):
        
        self.name = name
        self.format = name.rsplit('.',1)[1]
        
        # Read FASTA or GENKBANK formats

        if(self.format in dict_FASTA):
            self.__readFASTA(self.name)
        elif(self.format in dict_GB):
            self.__readGB(self.name)
            
    ''' Read GENKBANK format '''
    
    # store data in attributes dict_data & features in dict_features
            
    def __readGB(self,filename):
        
        self.filename = filename
        section_id = []  # index at which section starts
        with open(self.name,"r") as fi:
            
            # select main rows
            for line,ln in enumerate(fi):
                if(ln[0].isalpha()):
                    section_id.append(line)
                    
                    # Read all lines in file
        with open(self.name,"r") as f:
            lines = f.readlines()  # read lines
        max_rows = len(lines)      # max number of rows in file
        
        # Add final row
        section_id.append(max_rows)
        # number of rows need to be read for each section
        read_rows = np.diff(section_id).tolist()
        section_id.pop() # return to original list¥
        
        # Data will be stored here 
        self.dict_gb = {}       # Dictionary containing GB data
        self.dict_features = {} # Dictionary to store features (part of dict_gb)
        self.dict_source = {}   # Dictionary to store sources (part of dict_gb)
        self.dict_ref = {}      # Dictionary to store references (part of dict_gb)
        
        def get_part(ids):
            
            reads_id = read_rows[ids] # number of rows to be read
            start_id = section_id[ids] # first row
            
            # LINE READ LOOP
            
            ii=0; line = ''; 
            lst_source = []; lst_features = []; lst_ref = [] # temp lists

            for i in range(0,reads_id): 
                
                ii+=1; 
                lline = lines[start_id+ii-1] 
                
                if(tline == 'SOURCE'):
                    lst_source.append(lline)
                    
                elif(tline == 'REFERENCE'):
                    lst_ref.append(lline)
                    
                elif(tline != 'FEATURES'):
                    lline = lline.replace(tline,'')
                    lline = lline.lstrip(); lline = lline.rstrip()
                    line+=lline
                else:
                    lst_features.append(lline)
                    
            max_features = len(lst_features)
            max_source = len(lst_source)
            max_ref = len(lst_ref)
            
            ''' READ SOURCE DATA IN GB FILE '''
            # return dictionary with all data relevant to source
            
            def parse_source(lst_source):
                
                lst_source_start = []  # list to store feature starting index
                lst_source_id = []     # list to store feature name
                
                # read these features only
                source_read = ['  ORGANISM  ']
                
                # [1] ADD FIRST LINE DATA TO DICTIONARY
                self.dict_source['main'] = lst_source[0].replace('SOURCE      ','').strip()
                
                # [2] ADD ALL SEGMENT DATA TO DICTIONARY 
                
                for num,line in enumerate(lst_source):
                    
                    # add source index
                    for source in source_read:
                        if(source in line):
                            lst_source_id.append(source.strip())
                            lst_source_start.append(num) 
                            
                    lst_source_start.append(max_source)                  # start of each feature
                    source_reads_id = np.diff(lst_source_start).tolist() # number of lines to be read
                    lst_source_start.pop()                               # return to original list
                    
                    # CYCLE THOUGH ALL SECTIONS
                    
                for jj,source in enumerate(lst_source_id):
                    
                    temp = [lst_source[lst_source_start[jj]+i] 
                        for i in range(0,source_reads_id[jj])]
                    
                    lst_temp = temp.copy(); temp = []
                    for ii,source_content in enumerate(lst_temp):
                        
                        # if subsection is organism
                        
                        if('  ORGANISM  ' in source_content):
                            start_id = ii;  trans_code = ""
                            for last in lst_temp[start_id:]:
                                
                                # remove name
                                if('  ORGANISM  ' in last):
                                    last = last.replace('  ORGANISM  ','')
                                    
                                    # adjust string; sum
                                last = last.strip(); last = last.lstrip()
                                last = last.rstrip()
                                trans_code += last
                                
                            self.dict_source['organism'] = trans_code
                            
                            ''' PARSE REFERENCES '''  
                            # return a dictionary with all data relevant to references
                            
            def parse_references(lst_ref):
                
                lst_ref_start = []  # list to store feature starting index
                lst_ref_id = []     # list to store feature name
                
                # read these features only
                ref_read = ['  AUTHORS   ',
                            '  TITLE     ',
                            '  JOURNAL   ',
                            '   PUBMED   ']
                
                # [1] ADD FIRST LINE DATA TO DICTIONARY
                self.dict_ref['main'] = lst_ref[0].replace('REFERENCE   ','').strip()
                
                # [2] ADD ALL SEGMENT DATA TO DICTIONARY 
                
                for num,line in enumerate(lst_ref):
                    
                    # add source index
                    for ref in ref_read:
                        if(ref in line):
                            lst_ref_id.append(ref.strip())
                            lst_ref_start.append(num) 
                            
                    lst_ref_start.append(max_ref)                  # start of each feature
                    ref_reads_id = np.diff(lst_ref_start).tolist() # number of lines to be read
                    lst_ref_start.pop()                               # return to original list
                    
                    # CYCLE THOUGH ALL SECTIONS
                    
                for jj,ref in enumerate(lst_ref_id):
                    
                    temp = [lst_ref[lst_ref_start[jj]+i] 
                        for i in range(0,ref_reads_id[jj])]
                    
                    lst_temp = temp.copy(); temp = []
                    for ii,ref_content in enumerate(lst_temp):
                        
                        # add sections to dictionary
                        
                        def add_to_dict(subsection,name):
                            
                            if(subsection in ref_content):
                                start_id = ii;  trans_code = ""
                                for last in lst_temp[start_id:]:
                                    
                                    # remove name
                                    if(subsection in last):
                                        last = last.replace(subsection,'')
                                        
                                        # adjust string; sum
                                    last = last.strip(); last = last.lstrip()
                                    last = last.rstrip()
                                    trans_code += last
                                    
                                self.dict_ref[name] = trans_code
                                
                        add_to_dict('  AUTHORS   ','authors')
                        add_to_dict('  TITLE     ','title')
                        add_to_dict('  JOURNAL   ','journal')
                        add_to_dict('   PUBMED   ','pubmed')
                        
                        ''' READ FEATURES DATA IN GB FILE '''
                        # return dictionary with all data relevant to features
                        
            def parse_features(lst_features):
                
                lst_feat_start = []   # list to store feature starting index
                lst_feat_id = []      # list to store feature name
                
                # read these features only
                features_read = ['    gene            ',
                                '    source          ',
                                '    tRNA            ',
                                '    CDS             ']
                
                # [2] CYCLE THROUGH ALL FEATURE LINES ONLY
                
                # line -> individual line in feature list
                
                for num,line in enumerate(lst_features):
                    
                    # add feature index
                    for feature in features_read:
                        if(feature in line):
                            lst_feat_id.append(feature.strip())
                            lst_feat_start.append(num) 
                            
                    lst_feat_start.append(max_features) # start of each feature
                    feat_reads_id = np.diff(lst_feat_start).tolist() # number of lines to be read
                    lst_feat_start.pop() # return to original list
                    
                    # [!] Cycle through all features in file
                    
                for jj,feature in enumerate(lst_feat_id):
                    
                    # [!] list of feature data [all]
                    
                    temp = [lst_features[lst_feat_start[jj]+i] 
                        for i in range(0,feat_reads_id[jj])]
                    
                    # feature name for dictionay entry 
                    feat_name = f'{jj}_{feature}'
                    
                    # add feature index
                    
                    for feature in features_read:
                        if(feature in temp[0]):
                            temp[0] = temp[0].replace(feature,'')
                            
                    lst_temp = temp.copy(); temp = []; feat_dict = {}
                    for ii,feat_content in enumerate(lst_temp):
                        
                        tline = feat_content.strip()
                        tline = tline.lstrip()
                        tline = tline.rstrip()
                        
                        # translation is always last, store translation code
                        
                        if('/translation="' in tline):
                            start_id = ii;  trans_code = ""
                            for last in lst_temp[start_id:]:
                                last = last.strip(); last = last.lstrip()
                                last = last.rstrip()
                                trans_code += last
                            trans_code = trans_code.replace('/translation="','')
                            trans_code = trans_code.replace('"','')
                            feat_dict['translation'] = trans_code
                            
                            # Other Entries
                            
                        if(ii>0 and '=' in tline and '/translation="' not in tline):
                            segment_name = tline.split('=')[0]
                            segment_name = segment_name.replace('/','')
                            segment_val = tline.split('=')[1]
                            feat_dict[segment_name] = segment_val 
                        elif(ii==0):
                            feat_dict['location'] = tline.strip()
                        elif(ii>0 and '=' not in tline and '/translation="' not in tline):
                            feat_dict['option'] = tline.strip()
                            
                            # store feature dictionary inside dict for each feature
                    self.dict_features[feat_name] = feat_dict
                    
            if(tline == 'COMMENT'):
                line = line.strip()               # strip whitespace on left
                line = "".join(line.splitlines()) # merge multiple lines
                
            if(tline == 'SOURCE'):
                parse_source(lst_source)
                line = self.dict_source
                
            if(tline == 'REFERENCE'):
                parse_references(lst_ref)
                line = self.dict_ref
                
            if(tline == 'FEATURES'):
                parse_features(lst_features)
                line = self.dict_features    # return dictionary
                
            if(tline == 'ORIGIN'):
                line = line.lstrip()
                nonnumeric = ''.join([i for i in line if not i.isdigit()])
                nonnumeric = nonnumeric.lstrip()
                line = nonnumeric.replace(' ','')
                line = line.replace('//','')  # remove last line
                
            return line # dictionary/string
        
        # MAIN LOOP, READ ALL SECTIONS
        
        lst_refs = []
        for ii,i in enumerate(section_id):
            
            tline = lines[i].split('  ')[0] # section name
            res = get_part(ii)              # get the relevant data
            
            if('REFERENCE' in tline):
                lst_refs.append(res.copy())
                
        for ii, i in enumerate(section_id):
            tline = lines[i].split('  ')[0] # section name
            res = get_part(ii)              # get the relevant data
            
            if('REFERENCE' not in tline):
                self.dict_gb[tline] = res
            elif('REFERENCE' in tline):
                self.dict_gb['REFERENCE'] = lst_refs
                
    ''' Read FASTA format '''
    # store read data into:
    # lst_header -> contains header data, non parsed
    # lst_seq -> contains sequence data
        
    def __readFASTA(self,filename):
        
        tseq = None
        self.filename = filename                    # store filename
        self.lst_seq = []                           # list of sequences
        self.lst_header = []                        # list of sequence id
        format = filename.rsplit('.',1)[1]          # file format
        ff = dict_FASTA[format]                     # file format interpretation

        # Store data into header & seq lists
        with open(self.filename) as file:
            maxline =  sum(1 for line in file)

        with open(self.filename) as file:

            for ii,line in enumerate(file):
                
                if(line.startswith('>')):
                    
                    self.lst_header.append(line)
                    tseq = ""
                    
                else:

                    line = line.lstrip()
                    line = line.rstrip()
                    
                    if(line.isalnum()):
                        tseq += sub("\s","",line)

                        # last line in file is reached (no gap)
                        if(ii == maxline - 1):
                            self.lst_seq.append(tseq)
                    else:
                        # if gap 
                        self.lst_seq.append(tseq)
            
        print(f'[note] read -> FASTA [{ff}] | #seq: {len(self.lst_seq)}')
        
    ''' Store Sequences in Sequence Class Attributes '''
    
    # Once the files are read, they are stored in temporary attributes 
    # in class readSeq, method store stores them in class SQ attributes
        
    def store(self):
        
        # Create SQ instances from read FASTA data
        
        if(self.format in dict_FASTA):
            
            # if only one sequence in file
        
            if(len(self.lst_seq) == 1):
                
                lst_types = ['dna','rna','aa']
                for check in lst_types:
                    
                    # if sequence is valid
                    if(SQ(self.lst_seq[0],check).validate()): 

                        # if header doesn't contain | splits (just itself)
                        if(len(self.lst_header[0].split('|')) == 1):
                            ids = self.lst_header[0].split(' ')[0].replace('>','').lstrip()
                            desc = self.lst_header[0].split(' ',1)[1].lstrip()
                        else:
                            ids = self.lst_header[0].split('|')[-2].lstrip()
                            desc = self.lst_header[0].split('|')[-1].lstrip()

                        return SQ(seq=self.lst_seq[0],
                                  seq_type=check,
                                  origin=self.filename,
                                  id=ids,
                                  description=desc,
                                  valid=True)
                       
            # if more than one sequence in file
                       
            elif(len(self.lst_seq) > 1):
                
                lst_out = []     
                for i in range(0,len(self.lst_seq)):
                    lst_types = ['dna','rna','aa']
                    for check in lst_types:
                        
                        # if valid sequence
                        if(SQ(seq=self.lst_seq[i],seq_type=check).validate()):
                            
                            # if FASTA format
                            if(self.format in dict_FASTA):

                                ids = self.lst_header[i].split('|')[-2].lstrip().rstrip()
                                desc = self.lst_header[i].split('|')[-1].lstrip().rstrip()

                                lst_out.append(SQ(self.lst_seq[i],
                                                  seq_type=check,
                                                  origin=self.filename,
                                                  id=ids,
                                                  description=desc,
                                                  valid=True,
                                                ))
                                
                            elif(self.format in dict_GB):
                                print('[note] GB currently reads only one')
                                
                return lst_out  # return list of SQ files
            
        # Create SQ instances from read GENBANK data
                        
        elif(self.format in dict_GB):
            
            lst_types = ['dna','rna','aa']
            for check in lst_types:
                
                if(SQ(self.dict_gb['ORIGIN'],check).validate()):    # if valid sq
                    return SQ(seq=self.dict_gb['ORIGIN'],           # sequence string
                            seq_type=check,                         # sequence type
                            origin=self.filename,                   # SQ file generation origin
                            id = self.dict_gb['ACCESSION'],         # use accesion in identifier
                            description=self.dict_gb['DEFINITION'], # description
                            valid=True,                             # validated sequence
                            
                            locus = self.dict_gb['LOCUS'],          # locus tag
                            accession = self.dict_gb['ACCESSION'],  # accession identifier
                            version = self.dict_gb['VERSION'],      # version 
                            dblink = self.dict_gb['DBLINK'],        # dblink
                            keywords = self.dict_gb['KEYWORDS'],    # keyword
                            source = self.dict_gb['SOURCE'],        # source 
                            reference = self.dict_gb['REFERENCE'],  # references
                            comment = self.dict_gb['COMMENT'],      # comment
                            features = self.dict_gb['FEATURES'])    # sequence features