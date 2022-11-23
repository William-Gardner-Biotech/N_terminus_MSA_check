# Multiple Sequence Alignment outlier IDentifier (MSAoutID) 

Using modern genome assembly tools, scientists can annotate genomes in record time. This revolution in genome assembly is not without its own problems however. Weak signals along with incomplete sequencing reads can create a problem of unwanted gene merging. This is an error in identification of genes that includes a lot of noise from another gene by accidental merging. Using conservation of gene protein sequences (BLAST) and alignment tools we can catch these outliers and hopefully remove these imperfections from genomes. 

![image](https://user-images.githubusercontent.com/99355149/199628209-a2cd8edd-34a2-4bb7-9a82-445b64bbd27f.png)

My program seeks to solve this problem by using multiple sequence alignment data from hundreds of genomes to identify misidentified gene start sites. The standalone program is written in Python and run using the Linux command line. There is also a Bash script version avalable in the Program_with_bash folder but be advised: this is not being updated anymore as the Python standalone has become the only continuously updated/supported program.

Orthologous group Multiple Sequence Alignments are pulled from:
<ul>
  http://eggnog5.embl.de/
</ul>

# Using MSAoutID tool

1. Clone the respository into a directory of your choosing and ensure the machine using the program has an up to date version of Python.
2. Download and unzip multiple sequence alignment files from your desired Taxa is a designated directory
3. Run command line and navigate to the repository directory
4. Now run the program, including the following arguments:

       python3 MET_prog_sa.py (MSA direcory)  (Run option)  (Z-score Threshold)
       
       An example: ~/(Your_Directory)$ python3 MET_prog_sa.py mammalia/40674/ 2 3
       
5. View the processed results in the newly created processed/ directory.

## Run options

1. Comprehensive: 

![](https://i.imgur.com/KocOh5Z.png)

2. Trimmed:

![](https://i.imgur.com/B3tx2EJ.png)

3. CSV format:

![](https://i.imgur.com/XYqgWmz.png)



*Additional Notes:
       
Line 56 of MET_prog_sa.py includes a hard coded parameter to only run the program on MSA's with 10 or more proteins in the file. To remove this just open MET_lib_sa.py in a python editor and change the value. 
(WARNING: Lowering the threshold will increase run time and make meaningful analysis more difficult as statistical significance is decreased with sample size.)

Page last updated:
23-Nov-22
