Below there are two functions that differ slightly, both use the **NCBI API**, to search databases:
- <code>blast_aav1</code> - Uses <code>NCBIXML.read</code> in Biopython, thus don't obtain the **alignment alignment**
- <code>search_sequence.py</code> - Uses <code>SearchIO.read</code> in Biopython, which enables us to visualise the alignment using <code>view</code>

<code>search_sequende.py</code>
- contains a small class <code>BLASTwww</code>

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

# Visualise quary
plot = blast_query.view_alignment(0) # view first alignment from results dataframe
pn.pane.Bokeh(plot)

```
