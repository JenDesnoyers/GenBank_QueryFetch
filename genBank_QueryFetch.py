from Bio.Seq import Seq #Allows me to use .reverse_complement(), .translate(), and create Seq objects.
from Bio import SeqIO #Allows me to use SeqIO.parse to easily read a FASTA file, record.seq and record.ID to extract the sequence and ID portion of a FASTA file, respectively.
from Bio import Entrez #Allows me to use e-utilities like elink, esearch, and efetch in order to get sequences (nucleotide and amino acid) from GenBank.
from Bio.Blast import NCBIWWW #Allows me to access BLAST. 

#1. Using your desired text query, find six protein coding nucleotide sequences in Genbank. Fetch the sequences and write them into a fasta file called 'query_sequences_dna.fasta'.

Entrez.email = "desnoyej@uoguelph.ca" #This line is required in order to use Entrez functions. 
query="sonic hedgehog[All Fields] AND complete cds[All Fields] AND (animals[filter] AND biomol_genomic[PROP])" #This line allows me to search GenBank for the complete coding sequences.
#of the gene called Sonic Hedgehog. My search only returns the DNA sequence (not mRNA) of this gene in animals.
result = Entrez.esearch(db="nuccore",term=query) #Entrez.esearch uses the above query to search GenBank and pull the ID of sequences matching my query. 
record = Entrez.read(result) #In order to acually aqcuire the IDs for sequences matching my query, the program must use Entrez.read.
IDs = record["IdList"] #This line stores the IDs that were acquired above in variable called "IDs" so that they can be used later. 

hedgehog_sequence_dict = {} #Creating an empty dictionary to store IDs as the keys, and their respective sequence as the values.

for i in range(6): #The next block of code iterates through six of the IDs acquired above and retrieves their sequence using Entrez.efetch().
    hedgehog_seqs = Entrez.efetch(db="nuccore",id=IDs[i],rettype="fasta",retmode="text") #This line actually pulls the sequences based on the IDs provided. It uses a nucleotide database
    #in order to retrieve nucleotide sequences, and returns these sequences as fasta. 
    hedgehog_record = SeqIO.parse(hedgehog_seqs,"fasta") #This line parses through the sequences 'fetched' above and stores them in a variable called "hedgehog_record" so that they can be accessed later.
    hedgehog_key = IDs[i] #This line stores the IDs in the 'IDs' variable to another variable named 'hedgehog_key'. The 'hedgehog_key' will later be used to define the keys
        #in the dictionary called 'hedgehog_sequence_dict'. 
    for j in hedgehog_record: #This line iterates through the sequences stored in the 'hedgehog_record' variable. 
        hedgehog_value = j.seq #This line stores the sequences in the 'hedgehog_record' variable to another variable named 'hedgehog_value'. The 'hedgehog_value' will later be used to define the values
        #in the dictionary called 'hedgehog_sequence_dict'. 
        hedgehog_sequence_dict[hedgehog_key] = hedgehog_value #This line of code actually assigns the keys and values defined above to the 'hedgehog_sequence_dict' so that the dictinary
        #will contain IDs as the keys, and their associated sequences as the values. 

with open("query_sequences_dna.fasta",'w') as outfile_sequences: #This line opens the file called "query_sequences_dna.fasta" and assigns it to a variable called 'outfile_sequences'.
    #The 'with' portion will make sure that the "query_sequences_dna.fasta" file is closed after the code below this line has been executed. 
    for key, value in hedgehog_sequence_dict.items(): #This line iterates through each key:value pair from the 'hedgehog_sequence_dict' dictionary one at a time.
        outfile_sequences.write('>' + str(key) + '\n' + str(value) + '\n' + '\n') #This line writes a '>' character followed by a key from the 'hedgehog_sequence_dict' dictionary on one line and then
        #the corresponding value on the next line and then a blank line to separate key:value pairs to the file opened above. 

#2. Change your python script from Assignment #1 to a function such that the input is a sequence of interest and the output is the longest translated ORF.

#Making a function out of Assignment #1
#I ended up having to change a lot about my script from Assignment #1 because it wouldn't work with these new sequences. The main problem was that I originally used .index("TAA" or "TAG" or "TGA")
#to find the position of the first stop codon, however, this code really only looks for "TAA", which apparently wasnt' a problem with Assignment #1, but for these new sequences
#python ran into an issue where it couldn't find "TAA" for some of the ORFs and would throw an error. So, I managed to solve this issue by using .translate(), which I wasn't able to do before
#because somewhere along the way my sequence was no longer a sequence object, so I changed that using Seq and then I just let .translate() find the stop codons.
        
def translate_longest_ORF(ID, seq): #This line is require din order to define a function. The name is i blue, whereas the componetnts that need to be given to the function in order for it to work
    #are defined inside the brackets.
    #1. Read the contents of the FASTA file.
    #The file doesn't need to be read by the function because the file is being parsed outside the function in order to provide the ID and sequence to the function.  

    #2. Save descriptor line to a variable.
    descriptor = ID 
    
    #3. Store the DNA sequence (without descriptor line) as a string containing all the DNA bases.
    sequence = seq
        
    #4. For each of the six reading frames, find and translate the first ORF.
    #Reading frame #1:
    sequence1List = [] #This creates an empty list called 'sequence1List' so that I can later add codons to the list in order to build a sequence separated into codons.
    for i in range(0,len(sequence)-4,3): #This iterates through the 'sequence' variable starting at 0, stopping at the end of the sequence (found using len(sequence)), and iterates in steps of 3.
        #I had to add len(sequence)-4 because otherwise I ran into an index out of range error where my senquence length wasn't divisible by 3 (so python couldn't make a triplet at the end).
        triplet = (sequence[i])+(sequence[i+1])+(sequence[i+2]) #This takes the first i, second i and third i in each interation, adds them together to create a triplet, and then stores the triplet in a
        #variable called triplet. 
        sequence1List.append(triplet) #This stores every triplet in the 'sequence1List' in order to build a sequence seperated into codons. 

    start_Codon_Position = sequence1List.index("ATG") #This uses .index() to find the position of the very first start codon ("ATG") in the reading frame, and stores that position
    #in the variable 'start_Codon_Position'.

    cutSeq1 = sequence1List[int(start_Codon_Position):] #This line cuts the sequence up to (but not including) the position of the first start codon in the reading frame. This leaves a sequence, stored
    #in 'cutSeq1', that starts at the first start codon.

    sequence1 = "" #This empty string will later allow me to turn the list of codons (cutSeq1) into a string of nucleotides so that I can turn sequence1 into a Seq object and then use .translate().
    for i in cutSeq1: #iterating through the list of codons created above. 
        sequence1 = sequence1 + i #This line concatenates each codon from cutSeq1 into a string called sequence1. 
    sequence_object = Seq(sequence1) #This line converts sequence1, a string object, into a sequence object and stores it in a variable called "sequence_object".  
    aa_Sequence1 = sequence_object.translate(table='Standard', to_stop=True) #This line translates sequence_object. to_stop=True makes translation stop at the very first stop codon that is encountered. 

    #Reading frame #2:
    #Refer to reading frame #1 for exlaination of code. 
    sequence2List = [] 
    for i in range(1,len(sequence)-3,3):  
        triplet = (sequence[i])+(sequence[i+1])+(sequence[i+2]) 
        sequence2List.append(triplet) 
    
    start_Codon_Position2 = sequence2List.index("ATG")

    cutSeq2 = sequence2List[int(start_Codon_Position2):] 

    sequence2 = ""
    for i in cutSeq2:
        sequence2 = sequence2 + i
    sequence_object2 = Seq(sequence2)
    aa_Sequence2 = sequence_object2.translate(table='Standard', to_stop=True)

    #Reading frame #3:
    #Refer to reading frame #1 for exlaination of code. 
    sequence3List = [] 
    for i in range(2,len(sequence)-2,3): 
        triplet = (sequence[i])+(sequence[i+1])+(sequence[i+2]) 
        sequence3List.append(triplet) 
    
    start_Codon_Position3 = sequence3List.index("ATG") 

    cutSeq3 = sequence3List[int(start_Codon_Position3):] 

    sequence3 = ""
    for i in cutSeq3:
        sequence3 = sequence3 + i
    sequence_object3 = Seq(sequence3)
    aa_Sequence3 = sequence_object3.translate(table='Standard', to_stop=True)

    revsequence = sequence.reverse_complement() #This line uses the reverse_complement command from the Seq library from Biopython to get the reverse complement of my sequence
    #so I can access the last 3 reading frames.

    #Reading frame #4:
    #Refer to reading frame #1 for exlaination of code. 
    sequence4List = []
    for i in range(0,len(revsequence)-4,3):
        triplet = (revsequence[i])+(revsequence[i+1])+(revsequence[i+2])
        sequence4List.append(triplet)

    start_Codon_Position4 = sequence4List.index("ATG")

    cutSeq4 = sequence4List[int(start_Codon_Position4):]

    sequence4 = ""
    for i in cutSeq4:
        sequence4 = sequence4 + i
    sequence_object4 = Seq(sequence4)
    aa_Sequence4 = sequence_object4.translate(table='Standard', to_stop=True)

    #Reading frame #5:
    #Refer to reading frame #1 for exlaination of code. 
    sequence5List = []
    for i in range(1,len(revsequence)-3,3):
        triplet = (revsequence[i])+(revsequence[i+1])+(revsequence[i+2])
        sequence5List.append(triplet)

    start_Codon_Position5 = sequence5List.index("ATG")

    cutSeq5 = sequence5List[int(start_Codon_Position5):]

    sequence5 = ""
    for i in cutSeq5:
        sequence5 = sequence5 + i
    sequence_object5 = Seq(sequence5)
    aa_Sequence5 = sequence_object5.translate(table='Standard', to_stop=True)

    #Reading frame #6:
    #Refer to reading frame #1 for exlaination of code. 
    sequence6List = []
    for i in range(2,len(revsequence)-2,3):
        triplet = (revsequence[i])+(revsequence[i+1])+(revsequence[i+2])
        sequence6List.append(triplet)

    start_Codon_Position6 = sequence6List.index("ATG")

    cutSeq6 = sequence6List[int(start_Codon_Position6):]

    sequence6 = ""
    for i in cutSeq6:
        sequence6 = sequence6 + i
    sequence_object6 = Seq(sequence6)
    aa_Sequence6 = sequence_object6.translate(table='Standard', to_stop=True)

    
    #5. Find the longest translated ORF from #4.
    proteinList = [aa_Sequence1, aa_Sequence2, aa_Sequence3, aa_Sequence4, aa_Sequence5, aa_Sequence6] #Creates a list of each protein sequence from above so that it is easier to find the longest protein sequence.
    longestAminoAcidSequence = (max(proteinList, key=len)) #This line finds the longest protein sequence stored in the list named 'proteinList', and saves the longest sequence to the
    #variable 'longestAminoAcidSequence'.
    longest_aa = str(longestAminoAcidSequence) #This line converts the 'longestAminoAcid' variable to a tring so that it can be added to the dictionary named 'aa_Dict'.
    aa_Dict = {descriptor:longest_aa} #This line creates a dictionary named 'aa_Dict', where the descriptor of the sequence is the key, and the longest amino acid from that sequence is the value. 

    return(aa_Dict) #This line makes the 'translate_longest_ORF' function return aa_Dict so that it can be accessed outside of the function. This is particularily useful when writting the descriptor and longest amino
    #acid sequence to a file in order to create a fasta file. 

#End of Assignment #1
        
#3. Write a loop that takes the nuceotide sequences from (1) as input for the function and writes all resulting amino acid sequences to a fasta file called 'query_sequences_aa.fasta'.

with open("query_sequences_aa.fasta",'w') as outfile_aa: #This line opens the 'query_sequences_aa.fasta' file as the variable 'outfile_aa' for writting, and will close the file when the code block has
    #been executed. 
    for record in SeqIO.parse("query_sequences_dna.fasta", "fasta"): #This line parses through the 'query_sequences_dna.fasta' file, and recognizes it as a fasta file. It is important to parse this file
        #with Seq so that Seq can provide the descriptor and sequence from the fasta file to the 'translate_longest_ORF' function using record.id and record.seq, respectively. 
        get_aa_Dict = translate_longest_ORF(record.id, record.seq) #This line calls the 'translate_longest_ORF' function defined above, and provides the required sequece ID and sequence
        #using record.id and record.seq. Since the function returns a dictionary called aa_Dict, which contains ID's as the key and the translated longest ORF as the value, each time
        #'translate_longest_ORF' is called, the returned dictionary key:values are stored in the 'get_aa_Dict' dictionary.
        for key, value in get_aa_Dict.items(): #This line parses through each key:value pair in the 'get_aa_Dict' dictionary, one at a time. 
            outfile_aa.write('>' + str(key) + '\n' + str(value) + '\n' + '\n') ##This line writes a '>' character followed by a key from the 'get_aa_Dict' dictionary on one line and then
        #the corresponding value on the next line and then a blank line to separate key:value pairs to the file opened above. 
        
#4. Find the related protein sequences to the nucleotide sequences from (1) using ELink and write them to a fasta file called 'query_sequences_linked_aa.fasta'.
link_protein = Entrez.elink(db='protein',dbfrom='nuccore',id=IDs) #This line uses Entrez.elink() to find the protein IDs that are linked to the sequence ID's from part one. 
record_link = Entrez.read(link_protein) #This line extracts the results that ELink has found.

hedgehog_protein_IDs = [] #An empty list to store protein IDs acquired below.

#The next few lines assign the IDs extracted above to variables, and then appends these variables one by one to the 'hedgehog_protein_IDs' list. 
proteinID1 = record_link[1]['LinkSetDb'][0]['Link'][0]['Id']
hedgehog_protein_IDs.append(int(proteinID1))
proteinID2 = record_link[2]['LinkSetDb'][0]['Link'][0]['Id']
hedgehog_protein_IDs.append(int(proteinID2))
proteinID3 = record_link[3]['LinkSetDb'][0]['Link'][0]['Id']
hedgehog_protein_IDs.append(int(proteinID3))
proteinID4 = record_link[4]['LinkSetDb'][0]['Link'][0]['Id']
hedgehog_protein_IDs.append(int(proteinID4))
proteinID5 = record_link[5]['LinkSetDb'][0]['Link'][0]['Id']
hedgehog_protein_IDs.append(int(proteinID5))
proteinID6 = record_link[6]['LinkSetDb'][0]['Link'][0]['Id']
hedgehog_protein_IDs.append(int(proteinID6))

hedgehog_protein_dict = {} #Creates an empty dictionary that will later store the linked protein IDs as the keys and the corresponding amino acid sequences as the values
#so that they can later be written to a file.

#For the next two blocks of code, refer to section #1 of the assignment for what each code line is doing. This part is essentially the same except the program is using linked protein IDs to
#acquire the corresponding amino acid sequences, rather than using ID's to acquire the corresponding nucleotide sequences of the CDS.
for i in range(6):
    hedgehog_proteins = Entrez.efetch(db="protein",id=hedgehog_protein_IDs[i],rettype="fasta",retmode="text") 
    hedgehog_protein_record = SeqIO.parse(hedgehog_proteins,"fasta")
    hedgehog_protein_key = hedgehog_protein_IDs[i]
    for j in hedgehog_protein_record:
        hedgehog_protein_value = j.seq
        hedgehog_protein_dict[hedgehog_protein_key] = hedgehog_protein_value

with open("query_sequences_linked_aa.fasta",'w') as outfile_proteins:
    for key, value in hedgehog_protein_dict.items():
        outfile_proteins.write('>' + str(key) + '\n' + str(value) + '\n' + '\n')
       
#5. Using the linked protein sequences from (4) write a loop that performs Protein Blast searches one at a time against the PDB database and write the results to the file 'PDB_hits.txt'.
infile_proteins = open("query_sequences_linked_aa.fasta",'r') #This line opens the 'query_sequences_linked_aa.fasta' file created in section #4 for reading. 

with open("PDB_hits.txt",'w') as outfile_pdb: #This line opens the "PDB_hits.txt" file as outfile_pdb for writting.
    for record in SeqIO.parse("query_sequences_linked_aa.fasta", "fasta"): #This line parses through the "query_sequences_linked_aa.fasta" file created in section #4, which has been opened for reading.
        protein_sequences = record.seq #This line acquires the amino acid sequences from 'query_sequences_linked_aa.fasta' using record.seq
        #and then stores these sequences in the variable named 'protein_sequences'.
        result_handle = NCBIWWW.qblast("blastp", "pdb", protein_sequences, hitlist_size=5) #This line uses the sequences stored in 'protein_sequences' to perform blast searches using
        #NBCIWWW.qblast("blastp", "pdb", protein_sequences, hitlist_size=10) against the PDB database. I had to limit the number of hits written to the file using hitlist_size=10 because the default
        #of 50 was too much for my computer to handle and it would get stuck. It still takes a bit of time to run but if patient then it does work. 
        outfile_pdb.write(result_handle.read()) #This line writes the results from the blast searches to the 'PDB_hits.txt' file. 
