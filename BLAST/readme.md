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

|  ID         | ACCESSION   | DEFINITION                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |   ALI LENGTH |   HSP NUMER |     e-VALUE |   SCORE |   LENGTH |   IDENTITIES | SEQ-QUERY                                                 | SEQ-MATCH                                                 | SEQ-SUBJECT                                               |
| :-----------|:------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------:|------------:|------------:|--------:|---------:|-------------:|:----------------------------------------------------------|:----------------------------------------------------------|:----------------------------------------------------------|
|  pdb|6LO8|E | 6LO8_E      | Cryo-EM structure of the TIM22 complex from yeast [Saccharomyces cerevisiae S288C] >pdb|6LO8|G Cryo-EM structure of the TIM22 complex from yeast [Saccharomyces cerevisiae S288C] >pdb|6LO8|I Cryo-EM structure of the TIM22 complex from yeast [Saccharomyces cerevisiae S288C]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |           87 |           1 | 1.50081e-11 |     131 |       57 |           24 | QMRDSMNTYNNMVNRCFATCIRSFQEKKVNAEEMDCTKRCVTKFVGYSQRVALRFAE | QM+D M  Y+N+V RCF  C+  F   K+  +E  C  +C  KF+ +S+RV  RF E | QMKDFMRLYSNLVERCFTDCVNDFTTSKLTNKEQTCIMKCSEKFLKHSERVGQRFQE |
|  pdb|3DXR|A | 3DXR_A      | Crystal structure of the yeast inter-membrane space chaperone assembly TIM9.10 [Saccharomyces cerevisiae]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |           89 |           1 | 1.6426e-11  |     130 |       57 |           24 | QMRDSMNTYNNMVNRCFATCIRSFQEKKVNAEEMDCTKRCVTKFVGYSQRVALRFAE | QM+D M  Y+N+V RCF  C+  F   K+  +E  C  +C  KF+ +S+RV  RF E | QMKDFMRLYSNLVERCFTDCVNDFTTSKLTNKEQTCIMKCSEKFLKHSERVGQRFQE |
|  pdb|2BSK|A | 2BSK_A      | Crystal structure of the TIM9 Tim10 hexameric complex [Homo sapiens] >pdb|2BSK|C Crystal structure of the TIM9 Tim10 hexameric complex [Homo sapiens] >pdb|2BSK|E Crystal structure of the TIM9 Tim10 hexameric complex [Homo sapiens] >pdb|7CGP|D Cryo-EM structure of the human mitochondrial translocase TIM22 complex at 3.7 angstrom. [Homo sapiens] >pdb|7CGP|E Cryo-EM structure of the human mitochondrial translocase TIM22 complex at 3.7 angstrom. [Homo sapiens] >pdb|7CGP|F Cryo-EM structure of the human mitochondrial translocase TIM22 complex at 3.7 angstrom. [Homo sapiens] >pdb|7CGP|K Cryo-EM structure of the human mitochondrial translocase TIM22 complex at 3.7 angstrom. [Homo sapiens] >pdb|7CGP|L Cryo-EM structure of the human mitochondrial translocase TIM22 complex at 3.7 angstrom. [Homo sapiens] |           89 |           1 | 4.83462e-11 |     127 |       57 |           19 | QMRDSMNTYNNMVNRCFATCIRSFQEKKVNAEEMDCTKRCVTKFVGYSQRVALRFAE | Q ++ + TYN +   CF  C++ F  ++V  EE  C++ C+ K++  +QR+++RF E | QFKEFLGTYNKLTETCFLDCVKDFTTREVKPEETTCSEHCLQKYLKMTQRISMRFQE |
|  pdb|7CGP|J | 7CGP_J      | Cryo-EM structure of the human mitochondrial translocase TIM22 complex at 3.7 angstrom. [Homo sapiens]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |          103 |           1 | 5.58706e-06 |      95 |       50 |           16 | MRDSMNTYNNMVNRCFATCIRSFQEKKVNAEEMDCTKRCVTKFVGYSQRV        | +RD +  YN M   CF  C+ S   + ++AEE  C   C  K +  + R+        | LRDFLLVYNRMTELCFQRCVPSLHHRALDAEEEACLHSCAGKLIHSNHRL        |
|  pdb|3DXR|B | 3DXR_B      | Crystal structure of the yeast inter-membrane space chaperone assembly TIM9.10 [Saccharomyces cerevisiae]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |           95 |           1 | 0.000783618 |      80 |       37 |           15 | YNNMVNRCFATCIR-SFQEKKVNAEEMDCTKRCVTKF                     | +N +VN C+  CI  S+ E ++N  E  C  RCV K+                     | FNKLVNNCYKKCINTSYSEGELNKNESSCLDRCVAKY                     |
add Codeadd Markdown

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
