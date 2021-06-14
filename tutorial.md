# Tutorials
All used/downloaded files in this tutorial are deposed in folder "data". Outputs are deposed in folder data/outs.
## HMM 
### Downloading and modifying HMM files
For downloading hmm profiles you can use the same csv file as for clans/families downloading. For downloading hmm profiles you can use also a file with only the first column with family names. The header is obligatory. 
```python3
import hmm_download
list = hmm_download.load_data("data/pfam-seq.csv")
families = hmm_download.get_names(list)
hmm_download.download_hmm(families, "hmm_folder")
```
You can download hmm profiles of a few families directly from a list of Pfam accession numbers.
```python3
family = ["PF00042", "PF00002"] #globin and 7tm_2 pfam accession numbers
hmm_download.download_hmm(family, "hmm_folder")
```
The hmm profile flat file downloaded from Pfam could be wrapped into an object. The method get_length() returns the number of positions in the profile. Also, the initial probabilities are counted. 
```python3
import hammer_to_object
hmmObj = hammer_to_object.file_to_object("hmm_folder/PF00042.hmm")
print(hmmObj.get_length())
```
The output of the code above.
```python3
111
```
The object could be modified by adding a new position with the following match emission probabilities and insert emission probabilities of all amino acids, and transaction states probabilities of all states represented by a dictionary. All probabilities should be represented as -log0.25 of probability.
```python3
m = {'A': 2.61238, 'C': 4.6652, 'D': 2.90035, 'E': 2.58371, 'F': 3.40111, 'G': 3.24748, 'H': 3.72403, 'I': 2.7102, 'K': 2.54322, 'L': 2.36402, 'M': 3.69264, 'N': 3.13736, 'P': 3.53515, 'Q': 3.009, 'R': 2.70022, 'S': 2.65142, 'T': 2.85255, 'V': 2.54643, 'W': 5.22806, 'Y': 3.67926}
i = {'A': 2.68618, 'C': 4.42225, 'D': 2.77519, 'E': 2.73123, 'F': 3.46354, 'G': 2.40513, 'H': 3.72494, 'I': 3.29354, 'K': 2.67741, 'L': 2.69355, 'M': 4.2469, 'N': 2.90347, 'P': 2.73739, 'Q': 3.18146, 'R': 2.89801, 'S': 2.37887, 'T': 2.77519, 'V': 2.98518, 'W': 4.58477, 'Y': 3.61503}
t = {('m', 'm'): 0.36202, ('m', 'i'): 5.18438, ('m', 'd'): 1.21023, ('i', 'm'): 0.61958, ('i', 'i'): 0.77255, ('d', 'm'): 0.0, ('d', 'd'): None}
hmmObj.add_position(m, i, t)
hmmObj.get_length()
```
The output of the code above.
```python3
112
```
You can print out the object.
```python3
print(str(hmmObj))
```
First a few lines of output of the code above.
```python3
M - match emmision I - insert emmision T - state transition length =112
0	M	{'A': 2.41012, 'C': 4.57436, 'D': 2.9257, 'E': 2.66767, 'F': 2.91137, 'G': 3.04448, 'H': 3.24838, 'I': 3.00985, 'K': 2.51089, 'L': 2.4397, 'M': 3.76934, 'N': 3.09462, 'P': 3.54502, 'Q': 3.01212, 'R': 3.01928, 'S': 2.73624, 'T': 2.91348, 'V': 2.69511, 'W': 4.52529, 'Y': 3.66582}
	I	{'A': 2.68618, 'C': 4.42225, 'D': 2.77519, 'E': 2.73123, 'F': 3.46354, 'G': 2.40513, 'H': 3.72494, 'I': 3.29354, 'K': 2.67741, 'L': 2.69355, 'M': 4.2469, 'N': 2.90347, 'P': 2.73739, 'Q': 3.18146, 'R': 2.89801, 'S': 2.37887, 'T': 2.77519, 'V': 2.98518, 'W': 4.58477, 'Y': 3.61503}
	T	{('m', 'm'): 0.00555, ('m', 'i'): 5.59236, ('m', 'd'): 6.31471, ('i', 'm'): 0.61958, ('i', 'i'): 0.77255, ('d', 'm'): 0.0, ('d', 'd'): None}
1	M	{'A': 3.33262, 'C': 5.86832, 'D': 1.56356, 'E': 1.7629, 'F': 5.20656, 'G': 3.69721, 'H': 4.26823, 'I': 4.70125, 'K': 3.04485, 'L': 4.17119, 'M': 4.91753, 'N': 3.22983, 'P': 3.94972, 'Q': 1.60297, 'R': 3.55387, 'S': 2.70879, 'T': 2.55376, 'V': 4.25659, 'W': 6.29954, 'Y': 4.87028}
	I	{'A': 2.68618, 'C': 4.42225, 'D': 2.77519, 'E': 2.73123, 'F': 3.46354, 'G': 2.40513, 'H': 3.72494, 'I': 3.29354, 'K': 2.67741, 'L': 2.69355, 'M': 4.2469, 'N': 2.90347, 'P': 2.73739, 'Q': 3.18146, 'R': 2.89801, 'S': 2.37887, 'T': 2.77519, 'V': 2.98518, 'W': 4.58477, 'Y': 3.61503}
	T	{('m', 'm'): 0.00555, ('m', 'i'): 5.59236, ('m', 'd'): 6.31471, ('i', 'm'): 0.61958, ('i', 'i'): 0.77255, ('d', 'm'): 0.48576, ('d', 'd'): 0.9551}
2	M	{'A': 2.49951, 'C': 3.69605, 'D': 3.30642, 'E': 2.56902, 'F': 4.46199, 'G': 3.49121, 'H': 4.26288, 'I': 2.40988, 'K': 1.75931, 'L': 2.90479, 'M': 4.33838, 'N': 3.6181, 'P': 4.43702, 'Q': 3.4457, 'R': 1.93857, 'S': 3.26036, 'T': 3.41783, 'V': 2.57855, 'W': 3.8604, 'Y': 4.47897}
	I	{'A': 2.68618, 'C': 4.42225, 'D': 2.77519, 'E': 2.73123, 'F': 3.46354, 'G': 2.40513, 'H': 3.72494, 'I': 3.29354, 'K': 2.67741, 'L': 2.69355, 'M': 4.2469, 'N': 2.90347, 'P': 2.73739, 'Q': 3.18146, 'R': 2.89801, 'S': 2.37887, 'T': 2.77519, 'V': 2.98518, 'W': 4.58477, 'Y': 3.61503}
	T	{('m', 'm'): 0.00555, ('m', 'i'): 5.59236, ('m', 'd'): 6.31471, ('i', 'm'): 0.61958, ('i', 'i'): 0.77255, ('d', 'm'): 0.48576, ('d', 'd'): 0.9551}
```
And a few last lines.
```python3
109	M	{'A': 3.08846, 'C': 4.77368, 'D': 5.38738, 'E': 4.77213, 'F': 2.06296, 'G': 4.5922, 'H': 4.9213, 'I': 1.79291, 'K': 4.55798, 'L': 1.08668, 'M': 2.51104, 'N': 4.75929, 'P': 4.93488, 'Q': 4.66572, 'R': 4.55528, 'S': 3.90753, 'T': 3.64121, 'V': 2.60202, 'W': 4.35426, 'Y': 3.72979}
	I	{'A': 2.68618, 'C': 4.42225, 'D': 2.77519, 'E': 2.73123, 'F': 3.46354, 'G': 2.40513, 'H': 3.72494, 'I': 3.29354, 'K': 2.67741, 'L': 2.69355, 'M': 4.2469, 'N': 2.90347, 'P': 2.73739, 'Q': 3.18146, 'R': 2.89801, 'S': 2.37887, 'T': 2.77519, 'V': 2.98518, 'W': 4.58477, 'Y': 3.61503}
	T	{('m', 'm'): 0.00555, ('m', 'i'): 5.59236, ('m', 'd'): 6.31471, ('i', 'm'): 0.61958, ('i', 'i'): 0.77255, ('d', 'm'): 0.48576, ('d', 'd'): 0.9551}
110	M	{'A': 2.45439, 'C': 4.82492, 'D': 4.32274, 'E': 2.8472, 'F': 3.43029, 'G': 3.52542, 'H': 4.21011, 'I': 2.74114, 'K': 3.41666, 'L': 1.67102, 'M': 3.21975, 'N': 3.64086, 'P': 3.67427, 'Q': 3.33641, 'R': 3.3995, 'S': 3.3561, 'T': 2.43023, 'V': 1.84571, 'W': 5.42175, 'Y': 3.61143}
	I	{'A': 2.68618, 'C': 4.42225, 'D': 2.77519, 'E': 2.73123, 'F': 3.46354, 'G': 2.40513, 'H': 3.72494, 'I': 3.29354, 'K': 2.67741, 'L': 2.69355, 'M': 4.2469, 'N': 2.90347, 'P': 2.73739, 'Q': 3.18146, 'R': 2.89801, 'S': 2.37887, 'T': 2.77519, 'V': 2.98518, 'W': 4.58477, 'Y': 3.61503}
	T	{('m', 'm'): 0.00374, ('m', 'i'): 5.59055, ('m', 'd'): None, ('i', 'm'): 0.61958, ('i', 'i'): 0.77255, ('d', 'm'): 0.0, ('d', 'd'): None}
111	M	{'A': 2.61238, 'C': 4.6652, 'D': 2.90035, 'E': 2.58371, 'F': 3.40111, 'G': 3.24748, 'H': 3.72403, 'I': 2.7102, 'K': 2.54322, 'L': 2.36402, 'M': 3.69264, 'N': 3.13736, 'P': 3.53515, 'Q': 3.009, 'R': 2.70022, 'S': 2.65142, 'T': 2.85255, 'V': 2.54643, 'W': 5.22806, 'Y': 3.67926}
	I	{'A': 2.68618, 'C': 4.42225, 'D': 2.77519, 'E': 2.73123, 'F': 3.46354, 'G': 2.40513, 'H': 3.72494, 'I': 3.29354, 'K': 2.67741, 'L': 2.69355, 'M': 4.2469, 'N': 2.90347, 'P': 2.73739, 'Q': 3.18146, 'R': 2.89801, 'S': 2.37887, 'T': 2.77519, 'V': 2.98518, 'W': 4.58477, 'Y': 3.61503}
	T	{('m', 'm'): 0.36202, ('m', 'i'): 5.18438, ('m', 'd'): 1.21023, ('i', 'm'): 0.61958, ('i', 'i'): 0.77255, ('d', 'm'): 0.0, ('d', 'd'): None}
```
The modified (or not) object could be written into a file for further analysis. 
```python3 
hmmObj.file_format("hmm_folder/test.hmm")
```
### Search by hmmsearch and hmmerscan and automatisation
For basic analysis, by hmmsearch you can use the following command. To use more parameters see the documentation.  
```python3
import hmmer_command
hmmer_command.hmmsearch(o="out/hem.hmmsearchout", hmm_file="hmm_folder/PF00042.hmm", fasta_file="fasta_folder/example_hem.fasta")
```
For basic analysis by hmmscan first prepare files by hmmpress, the use hmmscan as in the example. To use more parameters see the documentation.  
```python3
import hmmer_command
#prepare files
hmmer_command.hmmpress(path="hmmpress", hmm_file="hmm_folder/PF00042.hmm")
#perform analisys
hmmer_command.hmmscan(o="out/hem.hmmscanout", hmm_file="hmm_folder/PF00042.hmm", fasta_file="fasta_folder/example_hem.fasta")
```
We can also perform the automatic analysis: all hmm files vs all fasta files. For using mode search type the following lines. 
In the folder example_hmm are files: PF00042.hmm and PF00002.hmm, and in folder example_fasta are files: example_gprot.fasta and example_hem.fasta which you can easily find in the data folder. Outputs as always are in data/outs.  
```python3
import hmm_autosearch
hmm_autosearch.automatic_search("hmm_folder", "fasta_folder", "outsearch", "search")
```
For scan mode, you can use the following command. This function includes preparing files for the hmmscan by hmmpress.
```python3
import hmm_autosearch
hmm_autosearch.automatic_search("hmm_folder", "fasta_folder", "outscan", "scan")
```


## Download pfam info

You can choose one of two option when downloading the information about pfam family: create the class family, which download all data automatically or by the function, which will download only requested information.

### Pfam sequences

Pfam sequences will automatically download all sequences format specified by the user:
* fasta, stockholm, msf format in SeqIO object
* selex format in str object
* full or seed format of sequences
All of these sequences can be downloaded or downloaded and saved as file.

To download all full sequences in globin family in fasta format type:
```python3
globin_seqs = family_seq(‘PF00042’, form='fasta', type='full')
```
The result is fully compatible with the [SeqIO package](https://biopython.org/wiki/SeqIO)


### Pfam tree

Similar to the pfam sequences the tree of the family can be downloaded in the Phylo format from Biopython package. Users can specify if they wish to download the newick format into the file.
To download the tree of the globin family:
```python3
globin_tree = family_tree(‘PF00042’)
```

The result is fully compatible with the [Phylo package](https://biopython.org/wiki/Phylo)

### Alternative terms

Sometimes we don’t have the accession name of the family. With get_alternative function, we can obtain the accession name with the only id, or get id with the only accession.
To get the accession for globin:
```python3
globin_acc = get_alternative(‘Globin’) 
print(globin_acc)
```

Output is:
```python3
 ’PF00042’
```

### GO terms

pfam_sequences can also obtain all GO terms which users can find on pfam family website. They are represented as a set of string GO terms.
To get GO terms fot he family:
```python3
globin_go = go_terms(’PF00042’)
print(globin_go)
```

The output will be:
```python3
('GO:0051920', 'GO:0055114')
```
