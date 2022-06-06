Below there are two functions that differ slightly, both use the **NCBI API**, to search databases:

<code>blast_query.py</code>
- contains a small function <code>blast_query</code>
- <code>NCBIXML.read</code> in Biopython

```python

str_seq = 'ATTAAAGGTTTATACCTTCCCAGG'
seq = SQ(str_seq)
lst_aa = seq.get_protein(min_size=1)

results_df = blast_query(lst_aa[0])

```

<code>search_sequence.py</code>
- contains a small class <code>BLASTwww</code>
- uses <code>SearchIO.read</code> in Biopython
- class inputs require sequences in biopython <code>Seq</code> or string formats

```python

''' Example Usage '''

# We can utlise the SQ class to get amino acid chain
str_seq = 'ATTAAAGGTTTATACCTTC'      # Sequence in string format
seq = SQ(str_seq)                    # define SQ class sequence
lst_aa = seq.get_protein(min_size=1) # extract amino acid chains from via translation
lst_aa                               # visualise all created amino acid chains

blast_query = BLASTwww(lst_aa[0])    # Instantiate class
blast_query.fetch()                  # BLAST query
blast_query.view_query()             # view BLAST query results (results stored in dataframe)\

# Visualise query
plot = blast_query.view_alignment(0) # view first alignment from results dataframe
pn.pane.Bokeh(plot)

```
