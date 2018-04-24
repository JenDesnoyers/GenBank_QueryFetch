## GenBank_QueryFetch
The code in this repository will fetch six nucleotide sequence results for a GenBank query, translate the sequences for all six reading frames and then returns the longest ORF. The script also finds PDB files for each amino acid sequence and returns the results as a .txt file. 

/////////////////////////////////////////////////////////////Notes//////////////////////////////////////////////////////////////////////////////////

The script is set up so that it uses my Entrez email and the query is set to find Sonic Hedgehog proteins, but both of these can be easily changed -> see Entrez.email and query near the top of the script. The variables will still have something to do with Sonic Hedgehog, but the file names are ambiguous enough that they will work with any query and shouldn't cause any confusion (ex. "query_sequences_dna.fasta").  

Also, this script was made of an assignment for my masters degree, so thats why there are instructions in the code. I left them in because I think they help to explain what I was trying to accomplish at each part. 

///////////////////////////////////////////////////////////Dependencies//////////////////////////////////////////////////////////////////////////////
#Python 3
#BioPython 

