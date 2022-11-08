import MET_library
import sys
import os

# Allows ability to run in terminal and from IDE

if len(sys.argv) > 1:
    filename = sys.argv[1]
    # Marks the files
else:
    #filename = "FASTA_oma/HOG-B0561591.txt"
    filename = "mammalia/40674/3JI04.raw_alg.faa"


#show variable lets us see the sequences/you can also just write "T/F" in the functions#
show = False

seqs = MET_library.txt_to_seq(filename, show)

if len(seqs) < 10:
    #print("N < 10")
    exit()

'''sys.argv.append("Random")
sys.argv.append("Random")'''

'''if len(sys.argv) > 1:
    print("Alignment", sys.argv[2])
else:
    print("Alignment:", filename)
'''
PTs = MET_library.ProtIDs(filename)

seq_len = len(seqs[0])

sorted_positions = (MET_library.regx_firM(seqs))

triangle_groups = MET_library.build_triangles(sorted_positions, seq_len)

best_M_start = (max(triangle_groups, key=triangle_groups.get)-1)

#print("\nBEST M", best_M_start,"\n")

N_ter_only = MET_library.gaps(seqs, best_M_start)

Outlier = min(N_ter_only, key=N_ter_only.get)

Z_score = MET_library.Z_score(N_ter_only)

if type(Z_score) == str:
    print(Z_score)
    exit()

# Recursively build a list of all the bad alignments

'''print("ARG0", sys.argv[0])
print("ARG1", sys.argv[1])
print("ARG2", sys.argv[2])
print("ARG3", sys.argv[3])'''

if Z_score > 3.0:
    pt = PTs[Outlier]
    # Data Vis Portion
    MET_library.summary(seqs, pt, sorted_positions, N_ter_only, Outlier, sys.argv[3], Z_score, int(sys.argv[2]))
    # MET_library.summary(seqs, pt, sorted_positions, N_ter_only, Outlier, sys.argv[1], Z_score, sys.argv[2])
