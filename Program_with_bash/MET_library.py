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
    Prots = []
    file_name = open(filename, 'r')
    for line in file_name:
        line = line.strip("\n")
        if ">" in line:
            Prots.append(line)
        else:
            continue
    return Prots


# Going over the files to identify how many M's are in each position
### IMPORTANT NOTE: The AA's are in computational counting ie 0 = 1, so they don't line up to a visual MSA precisely
def met_pos_dict(sequences, show):
    AAs_noM = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "F", "P", "S", "T", "W", "Y", "V"]
    sorted_positionals = {}
    positionals = {}
    for seq in sequences:
        index = 1
        upstream_penalty = 0
        M_score = 2*(len(seq)//100)
        for aa in seq:
            if aa == "M":
                if index in positionals:
                    # Comment out second line if you don't want to penalize
                    positionals[index] += 1
                    #positionals[index] += (M_score - upstream_penalty)
                    #upstream_penalty = 0
                    index +=1
                else:
                    positionals[index] = 1
                    index += 1
            elif aa in AAs_noM:
                index += 1
                upstream_penalty += 1
            else:
                index += 1
                continue

    for key in sorted(positionals.keys()):
        sorted_positionals[key] = positionals[key]
        if show == True:
            print(key, positionals[key])

    return sorted_positionals

## Data Visualization Section for the first part of the pipeline
def first_step_vis(sequences):
    print("/\n/\n/\n/\n/\n/\n/\n/\n/\n/\n/\n/\n")

    print("Number of sequences                            :", len(sequences))
    print("Length of individual sequence                  :", len(sequences[0]))
    print("/\n/\n/\n/\n/\n/\n/\n/\n/\n/\n/\n/\n")

# Moving to identify outliers
# which Methionine position is the most conserved and closer to the left hand side

# Concatenating the sequences to the beginning portion upstream the most common M
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
        print(seq_no, i)
        gap_dict[seq_no] = gaps
        seq_no += 1

    return gap_dict

def second_step_vis(gaps, seqs):
    Outlier_1 = min(gaps, key=gaps.get)
    Mean_ahead = statistics.mean(gaps.values())
    s_dev = statistics.stdev(gaps.values())
    Z_score = (Mean_ahead-(gaps[Outlier_1]))/s_dev
    print("/\n/\n/\n/\n/\n/\n/\n/\n/\n/\n/\n/\n")
    print("Average number of gaps", Mean_ahead)
    print("Outlier: Sequence", Outlier_1, ", Gaps:", gaps[Outlier_1])
    print("\n------------------\n------St-Dev------\n", s_dev, "\n------------------")
    print("-----","Z-score-----\n",Z_score, "\n------------------\n")
    return Z_score
'''    if Z_score > 4.0:
        print("Flagged", Outlier_1, "Removed")
        indexable = Outlier_1 - 1
        seqs.pop(indexable)
        return main(seqs)
    else:
        print("Clean")'''

# Seqs should be a list variable
# RE variables . ^ $ * + ? { } [ ] \ | ( )
def regx_firM(seqs):
    seq_no = 0
    positions = {}
    sorted_positions = {}
    for line in seqs:
        seq_no+=1
        #print(line)
        N_ter = re.match("(([\-,A-L,N-Z]*M))", line)
        #print(seq_no,"NTER", N_ter.group(0))
        key = len(N_ter.group())
        if key in positions.keys():
            positions[key] += 1
        else:
            positions[key] = 1
    for key in sorted(positions.keys()):
        sorted_positions[key] = positions[key]
        print(key, positions[key])

    return sorted_positions

# Creating an object that groups nearby M values to account for minor shifts
class triangle():
    def __init__(self, center, frame):
        triangle.self = self
        triangle.center = center
        triangle.frame = frame

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
        #print("FFFF",F2C, "2up", 2**F2C)

        multiplier = 1/((2**F2C))
        #print("multiplier", multiplier)
        score += multiplier*positions_score[key]
        #score += multiplier*1
    return score