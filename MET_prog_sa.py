import MET_lib_sa
import sys
import os
import re
import shutil

### First checkpoint ###

if len(sys.argv) < 4:
    exit("Insufficient Arguments")

if len(sys.argv) > 4:
    exit("Too many arguments given")

run_option = sys.argv[2]

if run_option not in "0123":
    exit("Invalid Run Option provided, please select [1, 2, 3]")

run_dir = os.getcwd()

processed_path = (run_dir+"/processed/")

try:
    os.mkdir(processed_path)
except OSError as error:
    overwrite = str(input("A processed directory already exists, would you like to delete it? [Y/N]"))
    overwrite = overwrite.upper()
    #print(overwrite)
    if overwrite == "Y":
        shutil.rmtree(processed_path)
        os.mkdir(processed_path)
    #maybe make an option here for renaming processed directory
    elif overwrite == "N":
        print("Program terminated")
        exit()
    else:
        print("Unrecognized option, terminating program")
        exit()


if sys.argv[2] == '3':
    results_csv = open(processed_path + "Results_table.csv", "w")
    results_csv.write("Alignment,Outlier_Seq_No,Z-score,PT_ID")
    results_csv.close()

#print(sys.argv)
### sys.argv[0] = MET_program, [1] = infile, [2] = Development choice, [3] = Z-score.
# choosing display option
# Look into splitting PT_ID into two parts as the first part is with egnog id's i believe

for i in os.listdir(sys.argv[1]):
    outname = ""
    # run the program through
    filename = i
    filename = sys.argv[1]+i
    run_option = sys.argv[2]
    outname = re.match("[\w]*", i)
    outname = outname.group(0)
    #### we can check len of sequences ####
    if len(MET_lib_sa.txt_to_seq(filename)) < 10 or \
        MET_lib_sa.process(filename, outname, sys.argv[2], sys.argv[3]) == True: continue

    # CSV format
    if run_option == '3':
        body = MET_lib_sa.process(filename, outname, sys.argv[2], sys.argv[3])
        results_csv = open(processed_path+"Results_table.csv", "a")
        if body != None:
            results_csv.write("\n")
            results_csv.write(body)
            results_csv.close()
        else:
            results_csv.close()
            continue

    # Would make sense to check that they entered a valid run option

    else:
        outfile = (processed_path+outname+"_out.txt")
        contents = MET_lib_sa.process(filename, outname, sys.argv[2], sys.argv[3])
        if contents == True:
            continue
        else:
            output = open(outfile, "w")
            output.write(contents)
            output.close()
        #sys.stdout.close()

    # Allows ability to run in terminal and from IDE

    #show variable lets us see the sequences/you can also just write "T/F" in the functions#
    print(outname)