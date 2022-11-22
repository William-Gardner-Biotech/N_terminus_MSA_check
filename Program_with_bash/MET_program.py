import MET_library

filename = "FASTA_oma/HOG-B0561591.txt"

#show variable lets us see the sequences/you can also just write "T/F" in the functions#
show = False

seqs = MET_library.txt_to_seq(filename, show)

PTs = MET_library.ProtIDs(filename)

seq_len = len(seqs[0])

sorted_positions = (MET_library.regx_firM(seqs))

triangle_groups = MET_library.build_triangles(sorted_positions, seq_len)

best_M_start = max(triangle_groups, key=triangle_groups.get)

print("\nBEST M", best_M_start,"\n")

N_ter_only = MET_library.gaps(seqs, best_M_start)

Outlier = min(N_ter_only, key=N_ter_only.get)

# Make up for indexing issues
Outlier = Outlier-1

Z_score = MET_library.second_step_vis(N_ter_only, seqs)

# Recursively build a list of all the bad alignments

print("Outlier Seq No", Outlier, PTs[Outlier])





