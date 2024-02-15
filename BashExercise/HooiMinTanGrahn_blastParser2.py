"""
Title: Blast Parser Exercise
Created on 2024 - 01 -22
Author: Hooi Min Tan Grahn

Description:
    This script contains a possible solution for the Running Exercise III -
        Part 2.
    This solution is only an example of how to achieve the expected result, the
        exercise may be solved in many different ways.

    The script contains two independent programs:
        Identity score - Convert the identity score measures to similarity
        (genetic distances) between the individuals.
        Alignment score - Convert the alignment score measures to similarity
        (genetic distances) between the individuals.

    All two programs accept str file (in whihc contained identity and alignment
        scores from each parental choromosome) to run, the output files
        are in tab-deliminated files, if not provided the program generates
        an output file with a default name.

Procedure:
    1. Parse all arguments from the command line using argparse.
    2. Check if the input files exist and are in the corresponded formats.
    3. Run the block of code according to the program chosen by the user.
        a. Identity score - Checks if the user is using Identity score program.
        The column of identity score will be computed, then output as a
            2-dimensional matrix.
        b. Alignment score - Checks if the user is using Alignment score program.
        The column of alignment score will be computed, then output as a
            2-dimensional matrix.

List of functions: Operation actions of add, sub and div

Usage:
    python HooiMinTanGrahn_blastParser.py input_file HooiMinTanGrahn_output.txt

"""
#%%
import argparse, sys, os, re

usage = '''Tool box containing one program require a str input file to run):\n'''
program = """\
\tMost_similar\tConvert the (genetic distances) to ranked list.
"""

# Create the argparse functionality, the description will be displayed in the
# help option. formatter_class is used to read the string in raw format,
# allowing the use of \n and \t.
parser = argparse.ArgumentParser(description=usage + program,
                                 formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('-i', dest='infile', type=str, required=True,
                    help="Input file in str format")

# Add optional argument with flag to output file.
parser.add_argument('-o', dest='outfile', type=str, required=False,
                    default='output.txt',
                    help="Output file name")


# Parse the arguments and save it in a variable, this creates a class object.
args = parser.parse_args()

# Account for error
try:  # try to read the infile
    if not os.path.exists(args.infile):
        raise FileNotFoundError
except FileNotFoundError:
    print("The file {} was not found!".format(args.infile))
    sys.exit(1)


def parse_blast_output(file):
    list_t = [["#query", "target", "e-value", "score", "identity(%)"]]
    temp_list = []
    with open(file, 'r') as f:
        query = None
        target = ""
        evalue = ""
        score = ""
        identity = ""
        last_match = ""
        for line in f:
            if line.startswith("Query= "):
                query = line.split()[1]
                temp_list = [query]  # start a new temp_list with the new query
                last_match = "query"
            elif line.startswith(">"):
                if last_match == "identity":
                    temp_list = [query]
                target = line[1:].strip()
                temp_list.append(target)
                last_match = "target"
            elif "Score =" in line:
                if last_match == "identity":
                    temp_list = [query]
                    temp_list.append(target)
                pattern = r"Score = (.+) bits.+Expect = (.+),"
                match = re.search(pattern, line)
                score = match.group(1)
                evalue = match.group(2)
                # Append evalue first, then score
                temp_list.append(evalue)
                temp_list.append(score)
                last_match = "score"
            elif "Identities =" in line:
                identity = line.split("Identities =")[1].strip().split("(")[1].split(")")[0].replace('%', '')
                temp_list.append(identity)
                list_t.append(temp_list)
                last_match = "identity"
            elif "***** No hits found *****" in line:
                temp_list = [query, " ", " ", " ", " "]  # start a new temp_list with the new query
                list_t.append(temp_list)
        
        # add the last temp_list to list_t after the loop ends
        if temp_list:
            list_t.append(temp_list)
        
    list_out = []    
    for i in list_t:
        fix_identity = i[4]
        fix_score = i[3]
        i[3] = fix_identity
        i[4] = fix_score
        list_out.append(i)
            
    return list_out  # return the list_t
result = parse_blast_output(args.infile)



output = open(args.outfile, "w") # Specify the path user want to open
try: 
    output.write('\n'.join('\t'.join(str(item) for item in sublist) for sublist in result))  
except:
    print("Something went wrong when writing to the output file") # Error message
finally:
    output.close()
                
