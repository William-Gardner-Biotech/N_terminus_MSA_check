import statistics
import re
from Bio import SeqIO

# Import our Fasta file of MSA
# Add a ">" character so that the for loop will append final sequence

if __name__ == "__main__":
    print("Running in the library")

def fasta_parser(filename):
    sequences = []
    fasta_seqs = SeqIO.parse(filename, "fasta")
    for fasta in fasta_seqs:
        sequence = str(fasta.seq)
        sequences.append(sequence)
        #print(sequence)
    return sequences

def ProtIDs(filename):
    PT_names = ["0"]
    fasta_names = SeqIO.parse(filename, "fasta")
    for fasta in fasta_names:
        name = fasta.id
        PT_names.append(name)
    return PT_names


# Moving to identify outliers
# which Methionine position is the most conserved and closer to the left hand side

# Concatenating the sequences to the beginning portion upstream the most common M
# Good index starting at 1
def gaps(sequences, M_start):
    gap_dict = {}
    seq_no = 1
    for i in sequences:
        i = i[0:M_start+1]
        gaps = 0
        for space in i:
            if space == "-":
                gaps += 1
            else: continue
        #print(seq_no, i)
        gap_dict[seq_no] = gaps
        seq_no += 1

    return gap_dict

def Z_score(gaps):
    Outlier_1 = min(gaps, key=gaps.get)
    Mean_ahead = statistics.mean(gaps.values())
    s_dev = statistics.stdev(gaps.values())
    if s_dev == 0:
        Zscore = ('St. Dev == 0, crisp alignment at N-ter')
        return (Zscore)
    Zscore = (Mean_ahead - (gaps[Outlier_1])) / s_dev
    return Zscore

### IMPORTANT NOTE: The AA's are in computational counting ie 0 = 1, so they don't line up to a visual MSA precisely

# Seqs should be a list variable
# RE variables . ^ $ * + ? { } [ ] \ | ( )
# Find the first Methionine
def regx_firM(seqs):
    positions = {}
    sorted_positions = {}
    for line in seqs:
        N_ter = re.match("(([\-,A-L,N-Z]*M))", line)
        if N_ter is None:
            key = len(line)
        else:
            #print(N_ter.group(0))
            key = len(N_ter.group(0))
        if key in positions.keys():
            positions[key] += 1
        else:
            positions[key] = 1
    #print("\nPositions of First Methionines")
    for key in sorted(positions.keys()):
        sorted_positions[key] = positions[key]
        #print(key, positions[key])

    return sorted_positions

# Creating an object that groups nearby M values to account for minor shifts
class triangle():
    def __init__(self, center, frame):
        triangle.self = self
        triangle.center = center
        triangle.frame = frame

# Pretty confusing because positions_score is the same thing and sorted_positions
def build_triangles(positions_score, seq_len):
    one_perc = seq_len//100
    if one_perc % 2 == 0:
        one_perc += 1
    #break out of function as it doesn't make sense to window with these values
    elif one_perc == 1 or one_perc == 0:
        return positions_score
    for key in positions_score.keys():
        tri = triangle(key, [])
        if one_perc == 1:
            tri = triangle(key, [key])
            continue
        else:
            leg = one_perc//2
            for i in range((key-leg), key):
                #print(i)
                if i in positions_score.keys():
                    tri.frame.append(i)
                else:
                    continue
            #print("center",key)

            tri.frame.append(key)
            for i in range((key+1), (key+leg+1)):
                #print(i)
                if i in positions_score.keys():
                    tri.frame.append(i)
                else:
                    continue
        new_score = F2C_distance(positions_score, tri)
        #print(tri.frame)
        #print(tri.center, new_score)
            #print("END LOOP")
        positions_score[tri.center] = new_score

    return (positions_score)
'''
    for key in sorted(positions_score.keys()):
        print(key, positions_score[key])'''

# Frame to Center distance
# This function will take all the frames and then using an equation that takes score*(1/(2^n))
# and adds all the scores in the fram where n is the distance from the center of the triangle
def F2C_distance(positions_score, tri):
    score = 0
    for key in tri.frame:
        #print(tri.center, "-", key)
        F2C = abs(tri.center - key)

        multiplier = 1/((2**F2C))
        #print("multiplier", multiplier)
        score += multiplier*positions_score[key]
        #score += multiplier*1
    return score

def stat_vis(gaps, seqs):
    stat_string = ""
    #print("\n----------------\n---Statistics---\n----------------\n")
    stat_string+=f"\n------------------\n----Statistics----\n------------------\n"
    #print("Number of sequences                       :", len(seqs))
    stat_string+=f"Number of sequences                       : {len(seqs)}\n"
    #print("Length of individual sequence             :", len(seqs[0]), "\n")
    stat_string+=f"Length of each sequence                   : {len(seqs[0])}\n"
    Outlier_1 = min(gaps, key=gaps.get)
    Mean_ahead = statistics.mean(gaps.values())
    s_dev = statistics.stdev(gaps.values())
    if s_dev == 0:
        #Zscore = ('St. Dev == 0, crisp alignment at N-ter')
        return True
    Zscore = (Mean_ahead-(gaps[Outlier_1]))/s_dev
    #print("Average number of gaps", int(Mean_ahead))
    stat_string+=f"Average number of gaps                    : {int(Mean_ahead)}\n"
    #print("Outlier: Sequence", Outlier_1, ", Gaps:", gaps[Outlier_1])
    stat_string+=f"Outlier Sequence                          : {Outlier_1}\n"
    stat_string+=f"Gaps                                      : {gaps[Outlier_1]}"
    #stat_string+=f"{gaps[Outlier_1]}"
    #print("\n------------------\n------St-Dev------\n", s_dev, "\n------------------")
    stat_string+=f"\n------------------\n------St-Dev------\n{s_dev} \n------------------"
    #print("-----","Z-score-----\n",Zscore, "\n------------------\n")
    stat_string+=f"\n-----Z-score-----\n{Zscore}\n------------------\n"
    return stat_string


# run option modes: {0:debug}. {1:Full stats}, {2:Trimmed stats}, {3: .csv format}

def process(filename, outname, run_option, Z_choice):

    output = ""

    Z_choice = int(Z_choice)

    run_option = int(run_option)

# Length checker

    seqs = fasta_parser(filename)

    PTs = ProtIDs(filename)

    seq_len = len(seqs[0])

    sorted_positions = regx_firM(seqs)

    triangle_groups = build_triangles(sorted_positions, seq_len)

    best_M_start = (max(triangle_groups, key=triangle_groups.get)-1)

    #print("\nBEST M", best_M_start,"\n")

    N_ter = gaps(seqs, best_M_start)

    Outlier = min(N_ter, key=N_ter.get)

    Zscore = Z_score(N_ter)

    if type(Zscore) == str:
        return Zscore

    # Recursively build a list of all the bad alignments

    if Zscore > Z_choice:
        pt = PTs[Outlier]
        pt = pt[1:]
        if run_option == 0:
            print("Debugging mode:")
            seq_len = len(seqs[0])
            print("\nSEQUENCES")
            for i in seqs:
                print(i)
            print("\nN_TERM")
            for key in sorted_positions.keys():
                print(key, sorted_positions[key])
            seq_no = 1
            print("\nGAPS")
            for j in N_ter:
                print(seq_no, j)
                seq_no += 1
            print("\nTRIANGLES:")
            tria = build_triangles(sorted_positions, seq_len)
            for k in tria:
                print(k)
            print("\nSTATS")
            print(stat_vis(N_ter, seqs))
            print("\nZ-score\n")
            print(Zscore)
            print("\nOUTLIER")
            print("Outlier Seq Number:", Outlier, "\nOutlier Protein ID:", pt)
            return 'Debug Mode'
        elif run_option == 1:
            output+= f"Alignment {outname}\n------------------"
            output+= f"\nPositions of First Methionines [AA]:[Count]"
            for key in sorted_positions.keys():
                output+= f"\n{key}:{int(sorted_positions[key])}"
            output+=f"{stat_vis(N_ter, seqs)}"
            output+=f"\nOutlier ID: {pt}"
            #print("Outlier Seq Number:", Outlier, "\nOutlier Protein ID:", pt)
            return output
        elif run_option == 2:
            output += f"Alignment {outname}\n"
            output += f"Outlier Seq Number: {Outlier} \nOutlier Protein ID: {pt}"
            return output
        elif run_option == 3:
            #header = "Alignment,Outlier_Seq_No,Z-score,PT_ID"
            # print(header)
            Outlier = str(Outlier)
            Zscore = str(Zscore)
            body = outname + "," + Outlier + "," + pt + "," + Zscore
            print(body)
            return body
    else:
        return True