## Searching for similar sequence in a Database

- If we have a <code>biological sequence</code> (eg. ATTAAAGGTTTATACCTTCCCAGG)


- Below is a function & class that differ slightly, both use the **NCBI API**
- Both by default use BLAST program <code>blastp</code> & seach in database <code>pdb</code>

<code>blast_query.py</code>
- contains a small function <code>blast_query</code>
- uses <code>NCBIXML.read</code> in Biopython
- class inputs require sequences in biopython <code>Seq</code> or <code>string</code> formats

```python

''' Example Usage '''

str_seq = 'ATTAAAGGTTTATACCTTCCCAGG'
seq = SQ(str_seq)
lst_aa = seq.get_protein(min_size=1)

# Output DataFrame w/ BLAST results
results_df = blast_query(lst_aa[0])

```

           ID ACCESSION                                         DEFINITION  \
0         ---       ---                                                ---   
0  pdb|6LO8|E    6LO8_E  Cryo-EM structure of the TIM22 complex from ye...   
1  pdb|3DXR|A    3DXR_A  Crystal structure of the yeast inter-membrane ...   
2  pdb|2BSK|A    2BSK_A  Crystal structure of the TIM9 Tim10 hexameric ...   
3  pdb|7CGP|J    7CGP_J  Cryo-EM structure of the human mitochondrial t...   
4  pdb|3DXR|B    3DXR_B  Crystal structure of the yeast inter-membrane ...   

  ALI LENGTH HSP NUMER   e-VALUE  SCORE LENGTH IDENTITIES  \
0        ---       ---       ---    ---    ---        ---   
0         87         1       0.0  131.0     57         24   
1         89         1       0.0  130.0     57         24   
2         89         1       0.0  127.0     57         19   
3        103         1  0.000006   95.0     50         16   
4         95         1  0.000784   80.0     37         15   

                                           SEQ-QUERY  \
0                                                ---   
0  QMRDSMNTYNNMVNRCFATCIRSFQEKKVNAEEMDCTKRCVTKFVG...   
1  QMRDSMNTYNNMVNRCFATCIRSFQEKKVNAEEMDCTKRCVTKFVG...   
2  QMRDSMNTYNNMVNRCFATCIRSFQEKKVNAEEMDCTKRCVTKFVG...   
3  MRDSMNTYNNMVNRCFATCIRSFQEKKVNAEEMDCTKRCVTKFVGY...   
4              YNNMVNRCFATCIR-SFQEKKVNAEEMDCTKRCVTKF   

                                           SEQ-MATCH  \
0                                                ---   
0  QM+D M  Y+N+V RCF  C+  F   K+  +E  C  +C  KF+ ...   
1  QM+D M  Y+N+V RCF  C+  F   K+  +E  C  +C  KF+ ...   
2  Q ++ + TYN +   CF  C++ F  ++V  EE  C++ C+ K++ ...   
3  +RD +  YN M   CF  C+ S   + ++AEE  C   C  K +  ...   
4              +N +VN C+  CI  S+ E ++N  E  C  RCV K+   

                                         SEQ-SUBJECT  
0                                                ---  
0  QMKDFMRLYSNLVERCFTDCVNDFTTSKLTNKEQTCIMKCSEKFLK...  
1  QMKDFMRLYSNLVERCFTDCVNDFTTSKLTNKEQTCIMKCSEKFLK...  
2  QFKEFLGTYNKLTETCFLDCVKDFTTREVKPEETTCSEHCLQKYLK...  
3  LRDFLLVYNRMTELCFQRCVPSLHHRALDAEEEACLHSCAGKLIHS...  
4              FNKLVNNCYKKCINTSYSEGELNKNESSCLDRCVAKY  

<br>

<code>search_sequence.py</code>
- contains a small class <code>BLASTwww</code>
- uses <code>SearchIO.read</code> in Biopython
- class inputs require sequences in biopython <code>Seq</code> or <code>string</code> formats

```python

''' Example Usage '''

# We can utlise the SQ class to get amino acid chain
str_seq = 'ATTAAAGGTTTATACCTTC'      # Sequence in string format
seq = SQ(str_seq)                    # define SQ class sequence
lst_aa = seq.get_protein(min_size=1) # extract amino acid chains from via translation
lst_aa                               # visualise all created amino acid chains

blast_query = BLASTwww(lst_aa[0])    # Instantiate class
blast_query.fetch()                  # BLAST query
blast_query.view_query()             # view BLAST query results, stored in DataFrame

# Visualise query
plot = blast_query.view_alignment(0) # view first alignment from results dataframe
pn.pane.Bokeh(plot)

```
