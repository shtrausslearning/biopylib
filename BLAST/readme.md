Below there are two functions that differ slightly, both use the **NCBI API**, to search databases:
- <code>blast_aav1</code> - Uses <code>NCBIXML.read</code> in Biopython, thus don't obtain the **alignment alignment**
- <code>BLASTwww</code> - Uses <code>SearchIO.read</code> in Biopython, which enables us to visualise the alignment using <code>view</code>
