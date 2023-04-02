import MET_lib_sa
import os
import re
import shutil
import argparse

### Building argparse

parser = argparse.ArgumentParser(description='Program to find improperly annotated or outlier gene N-terminus regions using homologous groups.')

parser.add_argument('-f', '--filename', required=True, type=str, default=None,
	metavar='<str>', help='Location of folder containing FASTA (.faa) file(s) you want to process. [Must be unzipped]')

parser.add_argument('-o', '--outfmt', required=False, type=int, default=None,
	metavar='<int>', help='Output file format. (1 = full detail, 2 = shortened, 3 = .csv [Default = 2])')

parser.add_argument('-z', '--z_choice', required=False, type=int, default=None,
	metavar='<int>', help='Z_score for determining outlier cutoff. NOTE: Lower z-scores will increase run time. [Default = 3, Max value = 9]')

arg = parser.parse_args()

# Code in my default argument values

infile = arg.filename

print("BEFORE", arg.outfmt)
if arg.outfmt or arg.outfmt == 0:
    if str(arg.outfmt) not in '0123' or len(str(arg.outfmt)) > 1:
        exit('Invalid Format Option provided, please select [1, 2, 3]')
    else:
        run_option = arg.outfmt
else:
    print('WRONG WAY')
    run_option = 2

if arg.z_choice:
    if arg.z_choice >= 1 and arg.z_choice <= 9:
        z_choice = arg.z_choice
    else:
        exit('Chosen Z-score is not within 1-9')
else:
    z_choice = 3

## Could use the built in default argument ad make the z_choice a float to increase precision

### First checkpoint ###

print('Arguments', arg.outfmt, run_option, arg.z_choice, z_choice)

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


if run_option == 3:
    results_csv = open(processed_path + "Results_table.csv", "w")
    results_csv.write("Alignment,Outlier_Seq_No,Z-score,PT_ID")
    results_csv.close()

# [0] = MET_program, [1] = infile, [2] = Development choice, [3] = Z-score.
# choosing display option
# Look into splitting PT_ID into two parts as the first part is with egnog id's i believe

for i in os.listdir(infile):
    outname = ""
    # run the program through
    filename = infile+i
    outname = re.match("[\w]*", i)
    outname = outname.group(0)
    #### we can check len of sequences ####
    if len(MET_lib_sa.fasta_parser(filename)) < 10 or \
        MET_lib_sa.process(filename, outname, run_option, z_choice) == True: continue

    # CSV format
    if run_option == 3:
        body = MET_lib_sa.process(filename, outname, run_option, z_choice)
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
        contents = MET_lib_sa.process(filename, outname, run_option, z_choice)
        output = open(outfile, "w")
        output.write(contents)
        output.close()

    # Allows ability to run in terminal and from IDE

    #show variable lets us see the sequences/you can also just write "T/F" in the functions#
    print(outname)