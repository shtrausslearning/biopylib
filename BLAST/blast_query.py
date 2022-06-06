from Bio.Blast import NCBIXML 
from Bio.Blast import NCBIWWW 
from Bio.Seq import Seq
from Bio import SeqIO
import re
import pandas as pd

# Search NCBI BLAST Query
def blast_query(seq,program="blastp",database='blastp',
               verbose=False,n_top=5):
    
    # if input is a string
    if(type(seq) is str):
        seq = Seq(seq)
    
    res_handle = NCBIWWW.qblast(program,database,seq)
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

        lst_hit = []; lst_access = []; lst_alilen = [];
        lst_hsps = []; lst_eval = []; lst_score = [];
        lst_len = []; lst_ident = []; lst_def = [];
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
                for hsp in alignment.hsps:
                    print (f"e-VALUE: {hsp.expect}")

        dic_data = \
        {"ID":lst_hit,"ACCESSION":lst_access,"DEFINITION":lst_def,
         "ALI LENGTH":lst_alilen,"HSP NUMER":lst_hsps,"e-VALUE":lst_eval,
         "SCORE":lst_score,"LENGTH":lst_len,"IDENTITIES":lst_ident,
         "SEQ-QUERY":lst_quer,"SEQ-MATCH":lst_mat,"SEQ-SUBJECT":lst_sub,
        };df_data = pd.DataFrame(dic_data)

        lst_org = []
        for i in range(n_top):
            ali = record.alignments[i]
            defin = ali.hit_def
            x = re.search("\[(.*?)\]", defin).group(1)
            lst_org.append(x)

            if(verbose):
                print ("ORGANISMS:")
                print('--------------------------------------------------------\n')
                for s in lst_orgs: 
                    print(s)

        return df_data
