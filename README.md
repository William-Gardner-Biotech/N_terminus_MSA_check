# N-termini Multiple Sequence Alignment Scoring and proofreading tool 

Using modern genome assembly tools, scientists can annotate genomes in record time. This revolution in genome assembly is not without its own problems however. Weak signals along with incomplete sequencing reads can create a problem of unwanted gene merging. This is an error in identification of genes that includes a lot of noise from another gene by accidental merging. Using conservation of genes protein sequences (BLAST) and alignment tools we can catch these outliers and hopefully remove these imperfections from genomes. 

![image](https://user-images.githubusercontent.com/99355149/199628209-a2cd8edd-34a2-4bb7-9a82-445b64bbd27f.png)

My program seeks to solve this problem by using multiple sequence alignment data from hundreds of genomes to identify misidentified gene start sites. The program is primarily written in Python but there is also a Bash script to automate the analysis over large data sets.

Orthologous groups and Multiple Sequence Alignments are pulled from:
<ul>
  https://omabrowser.org/</ul>
<ul>
  http://eggnog5.embl.de/
</ul>

## This program is still in development and not meant for commerical use
