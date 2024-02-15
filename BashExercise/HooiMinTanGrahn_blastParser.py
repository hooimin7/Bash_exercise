"""
Title: Blast Parser Exercise
Created on 2024 - 01 -22
Author: Hooi Min Tan Grahn

Description:
    This script contains a possible solution for the Blast Parser Exercise.
    This solution is only an example of how to achieve the expected result, the
        exercise may be solved in many different ways.

    The script contains one program:
        Blast Parser - Convert the output (option -outfmt) to tab-delimited text file.
    

    This programs accept str file (an optional BLAST output file as input in 
        which contained query, targets, e-value, identity and scores from each 
        blast output) to run, the output files are in tab-deliminated files, if
        not provided the program generates an output file with a default name.

Procedure:
    1. Parse all arguments from the command line using argparse.
    2. Check if the input files exist and are in the corresponded formats.
    3. Run the block of code according to the program.

List of functions: parse the blast output

Usage:
    python HooiMinTanGrahn_blastParser.py -i input_file -o HooiMinTanGrahn_output.txt

"""
#%%
import argparse, sys, os, re # Import the modules needed

usage = '''Tool box containing one program require a str input file to run):\n'''
program = """\
\tBlastParser\tConvert the (str input) to text-delimited text file.
"""

# Create the argparse functionality, the description will be displayed in the
# help option. formatter_class is used to read the string in raw format,
# allowing the use of \n and \t.
parser = argparse.ArgumentParser(description=usage + program, 
                                 formatter_class=argparse.RawTextHelpFormatter) 
# Add required argument with flag to input file.
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

def parse_blast_output(file): # define a function to parse the blast output
    # create a list to store the output
    list_t = [["#query", "target", "e-value", "score", "identity(%)"]] 
    temp_list = [] # create a temp_list to store the information of each query
    with open(file, 'r') as f: # open the file
        query = None # initialize the query
        target = "" # initialize the target
        evalue = "" # initialize the evalue
        score = "" # initialize the score
        identity = "" # initialize the identity
        last_match = "" # initialize the last_match
        for line in f: # loop through each line
            if line.startswith("Query= "): 
                query = line.split()[1] # get the query
                temp_list = [query]  # start a new temp_list with the new query
                last_match = "query" # update the last_match
            elif line.startswith(">"): # get the target
                if last_match == "identity": # if the last_match is identity, start a new temp_list
                    temp_list = [query] # start a new temp_list with the new query
                target = line[1:].strip() # get the target
                temp_list.append(target) # append the target to the temp_list
                last_match = "target" # update the last_match
            elif "Score =" in line: # get the score and evalue
                if last_match == "identity": # if the last_match is identity, start a new temp_list
                    temp_list = [query] # start a new temp_list with the new query
                    temp_list.append(target) # append the target to the temp_list
                pattern = r"Score = (.+) bits.+Expect = (.+)," # use regular expression to get the score and evalue
                match = re.search(pattern, line) # search the pattern in the line
                score = match.group(1) # get the score
                evalue = match.group(2) # get the evalue
                # Append evalue first, then score
                temp_list.append(evalue)
                temp_list.append(score)
                last_match = "score" # update the last_match
            elif "Identities =" in line: # get the identity
                identity = line.split("Identities =")[1].strip().split("(")[1].split(")")[0].replace('%', '') 
                temp_list.append(identity) # append the identity to the temp_list
                list_t.append(temp_list) # append the temp_list to the list_t
                last_match = "identity" # update the last_match
            elif "***** No hits found *****" in line: # if no hits found, start a new temp_list
                temp_list = [query, " ", " ", " ", " "]  # start a new temp_list with the new query
                list_t.append(temp_list) # append the temp_list to the list_t
        
        
    list_out = []   # create a list to store the output
    for i in list_t: # loop through each list in list_t
        fix_identity = i[4] # get the identity
        fix_score = i[3] # get the score
        i[3] = fix_identity # fix the identity
        i[4] = fix_score # fix the score
        list_out.append(i) # append the list to the list_out
            
    return list_out  # return the list_t
result = parse_blast_output(args.infile) # call the function and save the result in a variable
                
try:
    with open(args.outfile, "w") as output:
        output.write('\n'.join('\t'.join(str(item) for item in sublist) for sublist in result))  
        
except Exception as e:
    print("Something went wrong when writing to the output file")
    print(e) 