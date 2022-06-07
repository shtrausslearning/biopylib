## Searching for similar sequence in a Database

We can utilise <code>BLAST</code> to find the <code>sequence</code> or its subset in a databse & subsequently identify it, if found.


What we need:

- If we have a nucleotide sequence (eg. ATTAAAGGTTTATACCTTCCCAGG)
- Or an amino acid chain sequence (eg. HWLQMRDSMNTYNNMVNRCFATCIRSFQEKKVNAEEMDCTKRCVTKFVGYSQRVALRFAE)

We can utilise the <code>NCBI API</code>, via <code>BioPython</code>

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

<br>

<code>search_sequence.py</code>
- contains a small class <code>BLASTwww</code>
- uses <code>SearchIO.read</code> in Biopython
- class inputs require sequences in biopython <code>Seq</code> or <code>string</code> formats
- Result dataframe viewable via <code>view_query</code> method, contains an <code>alignment</code> column
- This is useful, as we can visualise it via <code>view_alignment</code> method

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
