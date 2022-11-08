import statistics
import re

# Import our Fasta file of MSA
# Add a ">" character so that the for loop will append final sequence

if __name__ == "__main__":
    print("Running in the library")

def txt_to_seq(filename, show):
    addition = open(filename, "a")
    addition.write(" ")
    addition.close()

    file_name = open(filename, 'r')  # ("Nog_seqs.txt", 'r')
    # Stats on this file: 102 individual sequences
    # Proteins of length 1070

    # Clean up the files
    sequences = []
    saved_line = ""
    for line in file_name:
        line = line.strip('\n')
        if line[0] == ">" or line [0] == " ":
            if saved_line == "":
                continue
            else:
                sequences.append(saved_line)
                if show == True:
                    print(saved_line)
                saved_line = ""
                continue
        else:
            saved_line = saved_line+line
    file_name.close()
    return sequences

def ProtIDs(filename):
    Prots = ["index_place_holder"]
    file_name = open(filename, 'r')
    for line in file_name:
        line = line.strip("\n")
        if ">" in line:
            Prots.append(line)
        else:
            continue
    return Prots


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

def stat_vis(gaps, seqs):
    print("\n----------------\n---Statistics---\n----------------\n")
    print("Number of sequences                       :", len(seqs))
    print("Length of individual sequence             :", len(seqs[0]), "\n")
    Outlier_1 = min(gaps, key=gaps.get)
    Mean_ahead = statistics.mean(gaps.values())
    s_dev = statistics.stdev(gaps.values())
    if s_dev == 0:
        Z_score = ('St. Dev == 0, crisp alignment at N-ter')
        return (Z_score)
    Z_score = (Mean_ahead-(gaps[Outlier_1]))/s_dev
    print("Average number of gaps", int(Mean_ahead))
    print("Outlier: Sequence", Outlier_1, ", Gaps:", gaps[Outlier_1])
    print("\n------------------\n------St-Dev------\n", s_dev, "\n------------------")
    print("-----","Z-score-----\n",Z_score, "\n------------------\n")
    return Z_score

def Z_score(gaps):
    Outlier_1 = min(gaps, key=gaps.get)
    Mean_ahead = statistics.mean(gaps.values())
    s_dev = statistics.stdev(gaps.values())
    if s_dev == 0:
        Z_score = ('St. Dev == 0, crisp alignment at N-ter')
        return (Z_score)
    Z_score = (Mean_ahead - (gaps[Outlier_1])) / s_dev
    return Z_score

### IMPORTANT NOTE: The AA's are in computational counting ie 0 = 1, so they don't line up to a visual MSA precisely

# Seqs should be a list variable
# RE variables . ^ $ * + ? { } [ ] \ | ( )
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

# Frame to Center disctance
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

# run option modes: {0:debug}. {1:Full stats}, {2:Trimmed stats}, {3: .csv format}
def summary(seqs, pt, sorted_positions, gaps, Outlier, sys_argv, Z_score, run_option):
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
        for j in gaps:
            print(seq_no, j)
            seq_no += 1
        print("\nTRIANGLES:")
        tria = build_triangles(sorted_positions, seq_len)
        for k in tria:
            print(k)
        print("\nSTATS")
        print(stat_vis(gaps, seqs))
        print("\nZ-score")
        print(Z_score)
        print("\nOUTLIER")
        print("Outlier Seq Number:", Outlier, "\nOutlier Protein ID:", pt[1:])
    elif run_option == 1:
        print("Alignment", sys_argv)
        print("\nPositions of First Methionines")
        for key in sorted_positions.keys():
            print(key, sorted_positions[key])
        stat_vis(gaps, seqs)
        print("Outlier Seq Number:", Outlier, "\nOutlier Protein ID:", pt[1:])
    elif run_option == 2:
        print("Alignment", sys_argv)
        print("Outlier Seq Number:", Outlier, "\nOutlier Protein ID:", pt[1:])
    elif run_option == 3:
        header = "Alignment,Outlier_Seq_No,Z_score,PT_ID"
        print("Alignment", sys_argv)
        print(header)
        Outlier = str(Outlier)
        Z_score = str(Z_score)
        body = sys_argv+","+Outlier+","+pt[1:]+","+Z_score
        print(body)
