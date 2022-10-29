import statistics

# Import our Fasta file of MSA
# Add a ">" character so that the for loop will append final sequence


def txt_to_seq(filename):
    addition = open(filename, "a")
    addition.write(">")
    addition.close()

    file_name = open(filename, 'r')  # ("Nog_seqs.txt", 'r')
    # Stats on this file: 102 individual sequences
    # Proteins of length 1070

    # Clean up the files
    sequences = []
    saved_line = ""
    for line in file_name:
        line = line.strip('\n')
        if line[0] == ">":
            if saved_line == "":
                continue
            else:
                sequences.append(saved_line)
                print(saved_line)
                saved_line = ""
                continue
        else:
            saved_line = saved_line+line

    return sequences

# Going over the files to identify how many M's are in each position
### IMPORTANT NOTE: The AA's are in computational counting ie 0 = 1, so they don't line up to a visual MSA precisely
def met_pos_dict(sequences):
    AAs_noM = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "F", "P", "S", "T", "W", "Y", "V"]
    positionals = {}
    for seq in sequences:
        #print("Cowboy", (len(seq)//100))
        index = 0
        upstream_penalty = 0
        for aa in seq:
            # print(index, aa)
            if aa == "M":
                if index in positionals:
                    positionals[index] += (((len(seq)//100)) - upstream_penalty)
                    upstream_penalty = 0
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
    return positionals

def total_M(position_dict):
    total_M = 0
    for key in sorted(position_dict.keys()):
        #print(key, position_dict[key])
        total_M += position_dict[key]

    return(total_M)


## Data Visualization Section for the first part of the pipeline
def first_step_vis(sequences, total_M):
    print("/\n/\n/\n/\n/\n/\n/\n/\n/\n/\n/\n/\n")

    print("Total # of Methionines across all sequences    :", total_M)
    print("Number of sequences                            :", len(sequences))
    print("Length of individual sequence                  :", len(sequences[0]))
    print("Average # of Methionines in a sequence         :", (total_M/len(sequences)))
    print("/\n/\n/\n/\n/\n/\n/\n/\n/\n/\n/\n/\n")

# Moving to identify outliers
# which Methionine position is the most conserved and closer to the left hand side

def sorts(positionals):
    # Let's make a new dictionary:
    sorted_dict = {}
    # work on balancing this
    for key in sorted(positionals.keys()):
        print(key, positionals[key])
        # The weighing equation
        #weighted_value = (positionals[key])*((seq_len-key)/seq_len)
        #print(key, "Weights", weighted_value)
        #sorted_dict[key] = weighted_value
        sorted_dict[key] = positionals[key]

    return sorted_dict

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
    return
'''    if Z_score > 4.0:
        print("Flagged", Outlier_1, "Removed")
        indexable = Outlier_1 - 1
        seqs.pop(indexable)
        return main(seqs)
    else:
        print("Clean")'''

# Perform a Z-score analysis #

# for square in

### Now we have a cleaner dictionary to window ###

### Current issues, grabbing an M from way too late, gotta weight this. Also what if there are multiple M's only 1 position off???? Combine the scroe maybe



class window_net:
    def __init__(self, left, right, innards, score):
        self.left = left
        self.right = right
        self.innards = innards
        self.score = score

def window_clean(weighted_position_dict, seq_len):
    max_pos = max(weighted_position_dict, key=weighted_position_dict.get)
    print("Max", max_pos, weighted_position_dict[max_pos])
    ### we want the window to only be 2% of the total sequence ###
    ### Arbitrarily set 2% by me ###
    cut_off = (weighted_position_dict[max_pos]) // int((seq_len) * 0.02)
    #print(cut_off)

    for key in sorted(weighted_position_dict.keys()):
        if weighted_position_dict[key] < cut_off:
            del weighted_position_dict[key]
        else:
            #print(key, weighted_position_dict[key])
            continue

    return (weighted_position_dict)

def window_generator(sorted_dict, seq_len):
    # We'll clean this in the future
    max_pos = max(sorted_dict, key=sorted_dict.get)
    cut_off = (sorted_dict[max_pos]) // int((seq_len) * 0.02)

    object_list = []

# Make this step recursive so we can collapse these objects down
    for key in sorted_dict.keys():
        curr_obj = window_net(0, None, [], 0)
        curr_obj.left = key
        curr_obj.score = sorted_dict[key]
        curr_obj.innards.append(key)
        for key in sorted(sorted_dict.keys()):
            if int(key) > int(curr_obj.left) and int(key) < (int(curr_obj.left) + cut_off):
                curr_obj.score += sorted_dict[key]
                curr_obj.innards.append(key)
                curr_obj.right = key
                print(curr_obj.innards)
            else:
                continue
        object_list.append(curr_obj)

    return object_list

def window_scoring(object_list, seq_len):
    #print(object_list)
    curr_max = 0
    max_pos = 0
    for obj in object_list:
        #print("left", obj.left)
        #print('old', obj.score)
        obj.score = (obj.score*(((seq_len-obj.left))**2))
        #print("new", obj.score)
        #We can now weight the scores to the left side of the window
        #print(obj.left, obj.right, obj.innards, obj.score)
        if obj.score > curr_max:
            curr_max = obj.score
            max_pos = obj.left
    print("M start after weighting and windowing:", max_pos, curr_max)

    return max_pos

filename = "FASTA_oma/HOG-B0579281.txt"
seqs = txt_to_seq(filename)

def main(seqs):
    met_dictionary = met_pos_dict(seqs)
    first_step_vis(seqs, total_M(met_dictionary))
    sorted_M_dict = sorts(met_dictionary)
    windows = window_generator(window_clean(sorted_M_dict, len(seqs[0])), len(seqs[0]))
    M_start_V2 = window_scoring(windows, len(seqs[0]))
    start_gaps = gaps(seqs, M_start_V2)
    second_step_vis(start_gaps, seqs)
    #print(int(len(seqs[0])*0.02))
    M_start = max(sorted_M_dict, key=sorted_M_dict.get)
    #print("Generic M start:", M_start)
    #print("Weighed M start:", M_start_V2)

main(seqs)