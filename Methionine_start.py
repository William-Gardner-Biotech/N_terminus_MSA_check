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
    positionals = {}
    for seq in sequences:
        index = 0
        for aa in seq:
            # print(index, aa)
            if aa == "M":
                if index in positionals:
                    positionals[index] += 1
                    index +=1
                else:
                    positionals[index] = 1
                    index += 1
            else:
                index += 1
                continue
    return positionals

def total_M(position_dict):
    total_M = 0
    for key in sorted(position_dict.keys()):
        print(key, position_dict[key])
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

def weights(positionals, seq_len):
    # Let's make a new dictionary:
    weighted_dict = {}
    # work on balancing this
    for key in sorted(positionals.keys()):
        # Length of sequence
        # 1 removed division by 0
        weighted_value = (positionals[key])*((seq_len-key)/seq_len)
        print(key, "Weights", weighted_value)
        weighted_dict[key] = weighted_value

    Consensus_M = max(weighted_dict, key=weighted_dict.get)

    return weighted_dict

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

def second_step_vis(gaps):
    Outlier_1 = min(gaps, key=gaps.get)
    Mean_ahead = statistics.mean(gaps.values())
    print("/\n/\n/\n/\n/\n/\n/\n/\n/\n/\n/\n/\n")
    print("Average number of gaps", Mean_ahead)
    print("Outlier: Sequence", Outlier_1, ", Gaps:", gaps[Outlier_1])
    return

# Perform a Z-score analysis #

# for square in

### Now we have a cleaner dictionary to window ###

### Current issues, grabbing an M from way too late, gotta weight this. Also what if there are multiple M's only 1 position off???? Combine the scroe maybe



class window_net:
    def __init__(self, left, right, score):
        self.left = left
        self.right = right
        self.score = score

def window_clean(weighted_position_dict, seq_len):
    max_pos = max(weighted_position_dict, key=weighted_position_dict.get)
    print("Max", max_pos, weighted_position_dict[max_pos])
    ### we want the window to only be 2% of the total sequence ###
    ### Arbitrarily set 2% by me ###
    cut_off = (weighted_position_dict[max_pos]) // int((seq_len) * 0.02)
    print(cut_off)

    for key in sorted(weighted_position_dict.keys()):
        if weighted_position_dict[key] < cut_off:
            del weighted_position_dict[key]
        else:
            print(key, weighted_position_dict[key])
            continue

    return (weighted_position_dict)

def window_generator(weighted_dict, seq_len):
    # We'll clean this in the future
    max_pos = max(weighted_dict, key=weighted_dict.get)
    cut_off = (weighted_dict[max_pos]) // int((seq_len) * 0.02)

    object_list = []

# Make this step recursive so we can collapse these objects down
    for key in weighted_dict.keys():
        curr_obj = window_net(0, None, 0)
        curr_obj.left = key
        curr_obj.score = weighted_dict[key]
        for key in sorted(weighted_dict.keys()):
            if int(key) > int(curr_obj.left) and int(key) < (int(curr_obj.left) + cut_off):
                curr_obj.score += weighted_dict[key]
                curr_obj.right = key
            else:
                continue
        object_list.append(curr_obj)

    print(object_list)
    for obj in object_list:
        print(obj.left, obj.right, obj.score)

def main():
    filename = "FASTA_oma/HOG-B0579281.txt"
    seqs = txt_to_seq(filename)
    met_dictionary = met_pos_dict(seqs)
    first_step_vis(seqs, total_M(met_dictionary))
    weighed_M = weights(met_dictionary, len(seqs[0]))
    #print(seqs)
    #print("Yo, display methionine positions", met_dictionary)
    M_start = max(weighed_M, key=weighed_M.get)
    print("Flag",M_start)
    start_gaps = gaps(seqs, M_start)
    second_step_vis(start_gaps)
    window_generator(window_clean(weighed_M, len(seqs[0])), len(seqs[0]))

    #print(int(len(seqs[0])*0.02))

main()