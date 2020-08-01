#tbox_pipeline_master.py
#By Merrick Pierson Smela
#Reads INFERNAL output data and calculates T-box features
#Also calculates thermodynamic parameters (code by Thomas Jordan)

import sys
import re
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import subprocess

#Function to read an INFERNAL output file and extract sequence names, metadata, structure, and sequence
#Metadata is as described in the INFERNAL manual
def read_INFERNAL(file):
    metadata = []
    Tbox_start = -1
    Tbox_end = -1
    sp = ""
    sp_start = ""

    metadataLine = -1
    structLine = -1
    seqLine = -1

    Tboxes = {'Name':[], 'Rank':[], 'E_value':[], 'Score':[], 'Bias':[], 'Tbox_start':[],'Tbox_end':[],
               'Strand':[], 'CM_accuracy':[], 'GC':[], 'Sequence':[], 'Structure':[]}

    with open(file) as f:
        for lineCount, line in enumerate(f):
            if re.match(">>", line) is not None: #We found a match!
                Tboxes['Name'].append(line.split(" ")[1])

                metadataLine = lineCount + 3
                structLine = lineCount + 6
                seqLine = lineCount + 9

            if lineCount == metadataLine:
                metadata = list(filter(None, line.split(' '))) #Splits by spaces, and strips empty strings

                Tboxes['Rank'].append(int(metadata[0][1:len(metadata[0])-1])) #All except first and last characters
                Tboxes['E_value'].append(float(metadata[2])) #third entry
                Tboxes['Score'].append(float(metadata[3])) #fourth
                Tboxes['Bias'].append(float(metadata[4])) #fifth
                Tbox_start = metadata[9]
                Tboxes['Tbox_start'].append(int(Tbox_start))
                Tbox_end = metadata[10]
                Tboxes['Tbox_end'].append(int(Tbox_end))
                Tboxes['Strand'].append(metadata[11])
                Tboxes['CM_accuracy'].append(float(metadata[13]))
                Tboxes['GC'].append(float(metadata[15][0:4])) #ignore the \n at the end

            if lineCount == structLine:
                sp = list(filter(None,line.split(' ')))
                structure = sp[0]
                Tboxes['Structure'].append(structure) #Take the second-to-last one

            if lineCount == seqLine:
                seq_line = list(filter(None, line.strip().split(' ')))
                sequence = ' '.join(seq_line[2:len(seq_line) - 1])

                #Fallback method (adapted from Thomas Jordan's code):
                if sequence not in line: #This can happen if there are two consecutive spaces in the middle of the sequence.
                    rsp = line.rsplit(Tbox_end, maxsplit = 1)
                    lsp = rsp[0].split(Tbox_start + ' ')
                    sequence = lsp[len(lsp) - 1].strip()
                    print("Fallback on line %d" % lineCount)

                #Do a sanity check
                if len(sequence) != len(structure): #This is an error!
                    print("\nParsing ERROR occured on line %d" % lineCount)
                    print(line) #For debug purposes
                    print(sequence)
                    print(structure)
                    print(seq_line)
                    #sequence = "ERROR"
                    
                Tboxes['Sequence'].append(sequence)
                
    return pd.DataFrame(Tboxes)

#Function to find features of a T-box given the secondary structure
#Parameters:
#seq is the sequence containing the T-box
#struct is the secondary structure
#start is the position of the T-box within the structure (default = 1)

#Output:
#Stem 1 start, stem 1 specifier loop, codon, stem 1 end, antiterminator start, discriminator, antiterminator end

def tbox_features(seq, struct, offset = 0):
    warnings = ""
    codon_region = ""
    
    if len(seq)!=len(struct):
        raise RuntimeError("Sequence length (%d) is not equal to structure length (%d)" % (len(seq), len(struct)))
        
    if not (struct.startswith(':') or struct.startswith('<')):
        warnings += "BAD_SEC_STRUCT;"
    
    #Find the Stem 1 start
    m = re.search('<',struct)
    if m is not None:
        s1_start = m.start() #Gets the first '<'
    else:
        s1_start = -1
        warnings += "NO_STEM1_START;"
        
    #Find the Stem 1 end
    m = re.search('>\.*,',struct)
    if m is not None:
        s1_end = m.start() #Gets the first '>,' or possibly '>.,'
    else:
        s1_end = -1
        warnings += "NO_STEM1_END;"
    
    #Find the Stem 1 specifier loop start
    matches = [m.start() for m in re.finditer('>\.*-', struct)]
    if len(matches) > 1:
        s1_loop_start = matches[1] #Gets the second occurrence of '>-' or possibly '>.-'
        #Find the Stem 1 specifier end
        matches = [m.end() for m in re.finditer('-\.*>', struct)]
        if len(matches) > 1:
            s1_loop_end = matches[1] #Gets the second occurrence of '->' or possibly '-.>'
        else:
            s1_loop_end = -1
            warnings += "NO_SPEC_END;"
    else:
        s1_loop_start = -1
        warnings += "NO_SPEC_START;"
        #Use fallback method to find Stem 1 specifier end
        s1_loop_end = -1
        for m in re.finditer('-\.*>', struct):
            end = m.end()
            if end > s1_loop_start and end < s1_end - 1: #The last loop STRICTLY before the stem 1 end
                s1_loop_end = end
        if s1_loop_end == -1: warnings += "NO_SPEC_END;"

    #Check to see if the Stem 1 has truncations
    if '~' in struct[s1_start:s1_end+1]: #there is a truncation
        warnings += "TRUNCATED_STEM_1;"
        #Recalculate Stem 1 loop end
        matches = [m.end() for m in re.finditer('-\.*>', struct[:s1_end+1])]
        if len(matches) > 1: #there should be at least 2
            s1_loop_end = matches[-2] #get the second to last
        if s1_loop_end == -1:
            warnings += "NO_SPEC_END;"
        else: #Recalculate Stem 1 loop start
            matches = [m.start() for m in re.finditer('>\.*-', struct[:s1_loop_end+1])]
            if len(matches) >= 1:
                s1_loop_start = matches[-1] #get the last one before the Stem 1 loop end

    if s1_loop_end > s1_loop_start: #Sanity check
        #Read the codon
        codon = seq[s1_loop_end - 5: s1_loop_end - 2]
        #Check the codon
        if re.search('[^AaCcGgUu]', codon) is not None:
            warnings += "BAD_CODON;"
        else: #Assign the codon region
            minus_one_pos = s1_loop_end - 6
            while minus_one_pos > 0 and re.match('[^AaCcGgUu]', seq[minus_one_pos]) is not None:
                minus_one_pos -= 1 #Look for the first ACGU character before the codon
            plus_one_pos = s1_loop_end - 2
            while plus_one_pos < len(seq) - 1 and re.match('[^AaCcGgUu]', seq[plus_one_pos]) is not None:
                plus_one_pos += 1 #Look for the first ACGU character after the codon
            codon_region = seq[minus_one_pos] + codon + seq[plus_one_pos] #Get the surrounding region too, for +1/-1 predictions
    else:
        codon = ""
        warnings += "NO_CODON;"
        
    #Find the antiterminator start
    antiterm_list = [m.start() for m in re.finditer(',\.*<', struct)] #Makes a list of all occurences of ',<'
    if len(antiterm_list) > 0:
        antiterm_start = antiterm_list[len(antiterm_list)-1] #Gets the last one
        discrim_start = struct.find('---', antiterm_start + 3) #Find the loop containing the discriminator
        discrim_end = discrim_start + 4
        discrim = seq[discrim_start:discrim_end]
    else:
        antiterm_start = -1
        discrim_start = -1
        discrim_end = -1
        discrim = ""
        warnings += "NO_ANTITERM_START;"
    
    #Check the discriminator
    if not discrim.startswith('UGG') or re.search('[^AaCcGgUu]', discrim) is not None:
        warnings += "BAD_DISCRIM;"
    
    #Find the antiterminator
    match = re.search('>\.*:',struct)
    if match is not None: #Sometimes the antiterminator end is missing from the sequence
        antiterm_end = match.start()
    else: #Simply get the last '>'
        end_list = [m.start() for m in re.finditer('>', struct)]
        if len(end_list) > 0:
            antiterm_end = end_list[len(end_list)-1]
        else:
            antiterm_end = -1
            warnings += "NO_ANTITERM_END;"
    
    #Adjust values based on offset
    s1_start += offset
    s1_loop_start += offset + 1
    s1_loop_end += offset - 1
    s1_end += offset
    antiterm_start += offset + 1
    antiterm_end += offset
    discrim_start += offset
    discrim_end += offset - 1
    
    #Return a tuple with the features identified
    return (s1_start, s1_loop_start, s1_loop_end, codon, s1_end, antiterm_start, discrim, antiterm_end, codon_region, warnings, discrim_start, discrim_end)

#Convert between position in INFERNAL output and fasta sequence
#Count the gaps and adjust for them
#Returns a mapping between INFERNAL position and fasta position
def map_fasta(seq, fasta, offset = 0, allowed = 'AaGgCcUu'):
    #Initialize counters
    count_fasta = offset
    parse_buffer = ""
    parsing = False
    mapping = []

    for c in seq:
        mapping.append(count_fasta)
        
        if parsing:
            if c == ']': #The end of the numerical gap
                parsing = False
                count_fasta += int(parse_buffer.strip()) #Parse the value
                #print(parse_buffer) #debug
                parse_buffer = "" #reset buffer for the next gap
            else:
                parse_buffer += c #Add the character to the parse buffer
            
        elif c in allowed:
            count_fasta += 1
            
        elif c == '[': # The start of a numerical gap
            parsing = True
        
    return mapping

#Function to find the end of the terminator (last occurence of TTTTT in the fasta sequence)
def term_end(sequence, start, pattern = 'TTTTT'):
    match = sequence.rfind(pattern, start) #get the last occurence of the pattern, after the start
    if match > 0:
        return match + len(pattern) #- 1
    return len(sequence) #fallback: return the end

#Function to compute derived T-box features from the prediction
#Note: tboxes must contain fasta sequences!
def tbox_derive(tboxes):
    #Derive more features for visualization
    #ALSO: Handle negative-strand T-boxes

    for i in range(len(tboxes['FASTA_sequence'])):
        fasta = tboxes['FASTA_sequence'][i]
        print('Mapping ' + tboxes['Name'][i]) #debug
        seq = tboxes['Sequence'][i]
        if isinstance(seq, str): #sanity check. Skip NaN
            #Check if the T-box is on the NEGATIVE strand
            if tboxes['Tbox_start'][i] > tboxes['Tbox_end'][i]:
                print("Converting – strand to +: " + tboxes['Name'][i])
                #Convert name
                split_name = tboxes['Name'][i].split(':')
                seq_start = split_name[1].split('-')[0]
                seq_end = split_name[1].split('-')[1]
                tboxes.at[tboxes.index[i], 'Name'] = split_name[0] + ':' + seq_end + '-' + seq_start
                #Convert FASTA sequence
                sequence = Seq(fasta, generic_dna)
                tboxes.at[tboxes.index[i], 'FASTA_sequence'] = str(sequence.reverse_complement())
                #Convert T-box start and end (since these are FASTA-relative)
                #Other features, which are INFERNAL-relative, should not be converted yet
                tboxes.at[tboxes.index[i], 'Tbox_start'] = len(fasta) - tboxes['Tbox_start'][i] + 1
                tboxes.at[tboxes.index[i], 'Tbox_end'] = len(fasta) - tboxes['Tbox_end'][i] + 1
                print("Conversion complete. New name is: " + tboxes['Name'][i])
                
            #Create mapping between INFERNAL sequence and FASTA sequence
            mapping = map_fasta(seq, tboxes['FASTA_sequence'][i], offset = tboxes['Tbox_start'][i] - 1)
            #Update the positions of existing features.
            s1_start = int(tboxes['s1_start'][i])
            if s1_start > 0:
                tboxes.at[tboxes.index[i], 's1_start'] = mapping[s1_start]

            s1_loop_start = int(tboxes['s1_loop_start'][i])
            if s1_loop_start > 0:
                tboxes.at[tboxes.index[i], 's1_loop_start'] = mapping[s1_loop_start]

            s1_loop_end = int(tboxes['s1_loop_end'][i])
            if s1_loop_end > 0:
                if s1_loop_end < len(mapping):
                    tboxes.at[tboxes.index[i], 's1_loop_end'] = mapping[s1_loop_end]
                    #Calculate codon range
                    tboxes.at[tboxes.index[i], 'codon_start'] = mapping[s1_loop_end - 4]
                    tboxes.at[tboxes.index[i], 'codon_end'] = mapping[s1_loop_end - 2]
                else:
                    print("Warning: mapping error for s1_loop_end:")
                    print(s1_loop_end)
                    print(len(mapping))
                    print(mapping)

            s1_end = int(tboxes['s1_end'][i])
            if s1_end > 0:
                tboxes.at[tboxes.index[i], 's1_end'] = mapping[s1_end]

            aterm_start = int(tboxes['antiterm_start'][i])
            if aterm_start > 0:
                tboxes.at[tboxes.index[i], 'antiterm_start'] = mapping[aterm_start]
            
            #Calculate discriminator range
            discrim_start = int(tboxes['discrim_start'][i])
            if discrim_start > 0:
                tboxes.at[tboxes.index[i], 'discrim_start'] = mapping[discrim_start]
                tboxes.at[tboxes.index[i], 'discrim_end'] = mapping[discrim_start + 3] #+3 because inclusive

            aterm_end = int(tboxes['antiterm_end'][i])
            if aterm_end > 0:
                aterm_end = min(aterm_end, len(mapping)-1)
                tboxes.at[tboxes.index[i], 'antiterm_end'] = mapping[aterm_end]
                #Calculate terminator end
                tboxes.at[tboxes.index[i], 'term_end'] = term_end(tboxes['FASTA_sequence'][i], int(tboxes['antiterm_end'][i]))
    return tboxes

def term_end_regex(sequence, start, pattern = '[T]{3,}[ACGT]{,1}[T]{1,}[ACGT]{,1}[T]{1,}'):
    print(start)
    if pd.isna(start):
        return len(sequence)
    matches = [m.end() for m in re.finditer(pattern, sequence[start:])]
    if len(matches)>0:
        return start + matches[0]
    return len(sequence) #fallback: return the end

#RNAfold on target sequence
def get_fold(sequence):
    #Initialize outputs
    structure = ""
    energy = ""
    errors = ""

    vienna_args = ['RNAfold', '-T', '37', '--noPS'] # arguments used to call RNAfold at 37 degrees
    vienna_input = str(sequence) # the input format
    vienna_call = subprocess.run(vienna_args, stdout = subprocess.PIPE, stderr = subprocess.PIPE, input = vienna_input, encoding = 'ascii')
    
    output = vienna_call.stdout.split('\n')
    if len(output) > 1: # if there is a result
        output = output[-2]
        output = output.split()
        structure = output[0] # gets last line's structure (always will be first element when sliced)
        energy = output[-1].replace(')', '').replace('(', '') # get energy (always will be last element in slice)
    errors = vienna_call.stderr.replace('\n',' ')
    return structure, energy, errors

#RNAfold on target sequence, with constraints
def get_fold_constraints(sequence, structure):
    #Initialize outputs
    energy = ""
    errors = ""
    structure_out = ""

    vienna_args = ['RNAfold', '-T', '37', '-C', '--noPS'] # arguments used to call RNAfold at 37 degrees with constraints
    vienna_input = str(sequence) + '\n' + str(structure) # the input format
    vienna_call = subprocess.run(vienna_args, stdout = subprocess.PIPE, stderr = subprocess.PIPE, input = vienna_input, encoding = 'ascii')
    
    output = vienna_call.stdout.split('\n')
    if len(output) > 1: # if there is a result
        output = output[-2]
        output = output.split()
        structure_out = output[0] # gets last line's structure (always will be first element when sliced)
        energy = output[-1].replace(')', '').replace('(', '') # get energy (always will be last element in slice)
    errors = vienna_call.stderr.replace('\n','') #changed from ' '
    return structure_out, energy, errors

#Make antiterminator constraints for folding
def make_antiterm_constraints(sequence, structure):
    if pd.isna(sequence) or pd.isna(structure):
        return None, None, None
    
    # search for antiterm end
    antiterm_end = term_end(structure, 0, pattern = ')') + 1 #search for last ')', which is the antiterm end

    # search for bulge of 4 or more bases
    match = re.search('[(][\.]{4,}[(]', structure)
    if(match):
        #Find TGGN
        discriminator = re.search('[T][G][G][ATGC]', sequence[(match.start() + 1):match.end()])
        # make hard constraint for unpaired 3' from antiterm end
        constraints = structure[:antiterm_end] + 'x' * len(structure[antiterm_end:])
        # make hard constraint for discriminator
        if(discriminator):
            constraints = constraints[:(discriminator.start() + match.start() + 1)] + 'xxxx' + constraints[(discriminator.end() + match.start() + 1):]
            structure_out, energy, errors = get_fold_constraints(sequence, constraints)
            pre_UGG = discriminator.start() + match.start()
            if (pre_UGG < 1) or (structure_out[pre_UGG] != '('): #check if base immediately before UGGN is paired
                if structure_out[0:4] != '((((': 
                    errors += "BAD_ANTITERM_STEM"
            return structure_out, energy, errors
    return None, None, None

#RNAeval to get energy
def get_energy(sequence, structure):
    if pd.isna(sequence) or pd.isna(structure):
        return None, ""

    vienna_args = ['RNAeval', '-T', '37'] # arguments that are used to call vienna RNAeval at T 37 degrees
    vienna_input = str(sequence + "\n" + structure) # the input format
    vienna_call = subprocess.run(vienna_args, stdout = subprocess.PIPE, stderr = subprocess.PIPE, input = vienna_input, encoding = 'ascii')
    # calls the subprocess with the vienna_input as input for the program
    energy = vienna_call.stdout.split()
    if energy: # if there actually is a result
        energy = energy[-1].replace('(', '').replace(')', '')
    errors = vienna_call.stderr.replace('\n',' ')
    return energy, errors

def get_sequence(fasta, start, end):
    if pd.isna(fasta) or pd.isna(start) or pd.isna(end):
        return None
    start, end = int(start), int(end)
    sequence = fasta[start:end]
    return sequence

def get_whole_structs(antiterm_start, whole_antiterm_structure, terminator_structure):
    if pd.isna(antiterm_start) or pd.isna(whole_antiterm_structure) or pd.isna(terminator_structure):
        return None, None
    start = int(antiterm_start) - 1 #NOTE: this is 1-indexed! if it's 0 indexed it will be off by one
    if start > 0:
        whole_term = whole_antiterm_structure[:start] + terminator_structure
        term_start = whole_term.find('(',start) + 1 #convert to 1-indexed
        return whole_term, term_start
    return None, None

def transcribe(sequence):
    # transcribes DNA to RNA
    tab = str.maketrans("ACGTacgt", "ACGUacgu")
    return sequence.translate(tab)

def replace_characters(structure):
    # replaces the weird infernal characters with normal dot brackets
    chars = ":,.~_-" # all characters representing unpaired structs, INCLUDING truncations (assumption)
    for c in chars:
        if c in structure:
            structure = structure.replace(c, ".") # creates a normal dot instead
    chars = "<{[(" # all characters representing paired bases
    for c in chars:
        if c in structure:
            structure = structure.replace(c, "(") # creates a normal open bracket
    chars = ">}])" # all characters representing paired bases
    for c in chars:
        if c in structure:
            structure = structure.replace(c, ")") # creates a normal closed bracket
    return structure

# removes truncations and adds the missing nucleotides from the FASTA
def resolve_truncations(fasta_sequence, processed_sequence, processed_structure):

    processed_length = len(processed_sequence)
    
    for truncations in re.finditer('\*\[.+?\][\*>]', processed_sequence):
        # finds the characters representing the truncation in the sequence
        # truncations have the form *[ N]* where N is a number of nucleotides
        # regexp : finds *[ followed by optional space, then number, then ]* or ]> if end of sequence
        #print(truncations.group(),truncations.start())
        
        difference = len(processed_sequence) - processed_length
        # as we are modifying the length of the processed_sequence inside the for loop, we need to adjust
        # our indices obtained from the regexp by how much the sequence length has changed
        
        missing = int(truncations.group().replace('*', '').replace('[', '').replace(']','').replace(' ', ''))
        # obtain the number of missing nucleotides by stripping away the useless characters
        
        match_start, match_end = truncations.start() + difference, truncations.end() + difference
        
        sequence_gaps = processed_sequence[:match_start].count('-')
        # number of gaps from start of the sequence to the truncation itself, important for alignment
        
        # obtain the positions of the truncations (adjusted by the difference in length
        # adds the missing n nucleotides to the sequence by getting them from the fasta sequence
        
        if(missing > 0):
            processed_sequence = processed_sequence[:match_start] + transcribe(fasta_sequence[(match_start - sequence_gaps):(match_start - sequence_gaps + missing)]) + processed_sequence[match_end:]
            processed_structure = processed_structure[:match_start] + "~" * (missing) + processed_structure[match_end:]
        else:
            processed_sequence = processed_sequence[:match_start] + processed_sequence[match_end:]
            processed_structure = processed_structure[:match_start] + processed_structure[match_end:]
    # appends the start of the T-box fasta if it is not in the infernal sequence
    return processed_sequence, processed_structure


#Makes a list of (lists of length 4) containing pair information
# in pairs, the first element of each entry is the character type and the second its index
# the third is the index of the bracket it is matched with and the fourth the matching bracket
def make_pairs(processed_structure, dot_structure):
    struct_list = list(dot_structure)
    length_list = len(struct_list) # length of the list struct_list
    pairs = [None] * length_list # creates a list of length dot structure filled with None
    
    for index, element in enumerate(struct_list):
        if(element == '('):
            pairs[index] = ['(', index, None, ''] # creates entries for the left brackets in pairs
        elif(element == ')'):
            pairs[index] = [')', index, None, ''] # creates entries for the right brackets in pairs

    for index, element in enumerate(pairs):
        l_count = 0
        r_count = 0
        i = index + 1
        if(element != None):
            if(element[0] == '('):
                l_count += 1
                while(i < length_list):
                    if(pairs[i] != None):
                        if('(' == pairs[i][0]):
                            l_count += 1
                        elif(')' == pairs[i][0]):
                            r_count += 1
                            l_count -= 1
                        if((r_count > 0) and (l_count == 0)):
                            element[2] = pairs[i][1]
                            element[3] = processed_structure[pairs[i][1]]
                            pairs[i][2] = element[1]
                            pairs[i][3] = processed_structure[element[1]]
                            break
                    i += 1
    return pairs

def replace_gaps(processed_sequence, processed_structure, dot_structure, pairs):
    structure_changed = False # check whether structure was modified by program
    corrected_brackets = False # check whether brackets were adjusted
    
    my_dot_structure = dot_structure # creates a copy of dot_structure that does not get changed
    my_processed_sequence = processed_sequence # creates a copy of processed_sequence that does not get changed
    
    for match in re.finditer('-', my_processed_sequence):
        index = match.start() # gets the index of each gap '-' in the sequence
        
        if(my_dot_structure[index] == '.'):
            processed_structure = processed_structure[:index] + 'x' + processed_structure[(index + 1):]
            dot_structure = dot_structure[:index] + 'x' + dot_structure[(index + 1):]
            # replaces the secondary structure with an x to tag it for deletion
            processed_sequence = processed_sequence[:index] + 'x' + processed_sequence[(index + 1):]
            # replaces the nucleotide in the sequence with an x, tagging it for deletion
        
        elif((my_dot_structure[index] == '(') or (my_dot_structure[index] == ')')):
            processed_structure = processed_structure[:index] + 'x' + processed_structure[(index + 1):]
            dot_structure = dot_structure[:index] + 'x' + dot_structure[(index + 1):]
            # replaces the open/closed bracket in the secondary structures by x to tag it for deletion
            processed_sequence = processed_sequence[:index] + 'x' + processed_sequence[(index + 1):]
            # replaces the corresponding nucleotide in the sequence by x to tag it for deletion
            
            structure_changed = True # we have changed the sequence
            
            if(isinstance(pairs[index][2], int) == True): # if there is a matching bracket
                paired = pairs[index][2] # the position of the paired closed/open bracket
                if(my_processed_sequence[paired] != '-'):
                # if the matching closed/open bracket is a real nucleotide in the sequence and not '-'
                    processed_structure = processed_structure[:paired] + '.' + processed_structure[(paired + 1):]
                    dot_structure = dot_structure[:paired] + '.' + dot_structure[(paired + 1):]
                    # replace the matching closed/open bracket by a '.'
                    corrected_brackets = True
    return processed_sequence, processed_structure, dot_structure, structure_changed, corrected_brackets

def parse_structures(fasta_sequence, tbox_start, tbox_end, sequence, structure):
    # Fills truncations, strips gaps, and aligns secondary structure
    # tbox_start and tbox_end are the positions of the tbox within the FASTA
    
    #Catch Nones and NaNs
    if any(pd.isna(arg) for arg in [fasta_sequence, tbox_start, tbox_end, sequence, structure]):
        return None, None, None
    
    tbox_start = int(tbox_start) # for some reason needed cuz its a float otherwise
    tbox_end = int(tbox_end)
    messages = '' # a string containing all messages that are not errors
    warnings = '' # a string containing all of the error messages
    processed_sequence = sequence # initializes the sequence for processing
    processed_structure = structure # initializes the structure for processing
    
    # count whether there is an unequal number of open and closed brackets
    count_open = structure.count('<')
    count_closed = structure.count('>')
    if(count_open != count_closed):
        warnings += ' unequal_number_of_brackets_in_input'
    
    #appends the start of the T-box fasta if it is not in the infernal sequence
    if(tbox_start > 1):
        processed_sequence = transcribe(fasta_sequence[:(tbox_start -1)]) + processed_sequence
        processed_structure = "~" * (tbox_start - 1) + processed_structure
        
    # appends the end of the T-box fasta if it is not in the infernal sequence
    if tbox_end < len(fasta_sequence):
        processed_sequence = processed_sequence + transcribe(fasta_sequence[tbox_end:])
        processed_structure = processed_structure + "~" * (len(fasta_sequence) - tbox_end)
    

    # Resolve internal truncations by filling with corresponding nucleotides from FASTA
    processed_sequence, processed_structure = resolve_truncations(fasta_sequence, processed_sequence, processed_structure)
    if(len(processed_sequence) != len(processed_structure)):
        warnings += ' unequal_lengths_structure_and_sequence'
        
    processed_sequence = processed_sequence.upper()
    dot_structure = replace_characters(processed_structure) # make a normal dot bracket structure
    
    pairs = make_pairs(processed_structure, dot_structure) # creates a list with all paired
    
    processed_sequence, processed_structure, dot_structure, structure_changed, corrected_brackets = replace_gaps(processed_sequence, processed_structure, dot_structure, pairs)
    if(structure_changed == True):
        messages += ' bracket_structure_was_changed'
    if(corrected_brackets == True):
        messages += ' changed_some_brackets_to_dots'
    
    #Delete gaps
    processed_sequence = processed_sequence.replace('x', '')
    processed_structure = processed_structure.replace('x', '')
    dot_structure = dot_structure.replace('x', '')
    
    if(transcribe(fasta_sequence) != processed_sequence):
        warnings += ' fasta_and_sequence_do_not_match'
    
    #Fix loops
    dot_structure = dot_structure.replace('(..)', '....')
    dot_structure = dot_structure.replace('(.)', '...')
    dot_structure = dot_structure.replace('(())', '....')
    
    # count whether there is an unequal number of open and closed brackets
    count_open = dot_structure.count('(')
    count_closed = dot_structure.count(')')
    if(count_open != count_closed):
        warnings += ' unequal_number_of_brackets'

    # make pairs again to process the secondary structure features of the Tbox
    # Find stems
    my_pairs = make_pairs(processed_structure, dot_structure)
    index = 0
    structures = []
    while(index < len(dot_structure) - 1):
        if(my_pairs[index] != None) and not None in my_pairs[index]:
            if(my_pairs[index][0] == '('):
                # this is a secondary structure element (a stem)
                #use 1-index
                structures.append([ my_pairs[index][1] + 1, my_pairs[index][2] + 1])
                index = my_pairs[index][2] # go to the paired bracket
        index += 1
    
    #Stem 3 is the second-to-last. Only calculate if there also a stem 1 (total stems >2)
    stem_three = None
    if(len(structures) > 2):
        stem_three = structures[-1]
    
    #Return: dot structure, list of secondary structural elements, messages and warnings
    return dot_structure, structures, (messages + warnings)

def local_fold(sequence, length):
    vienna_args = ['RNALfold', '-T', '37', '-L', str(length), '−−noClosingGU'] # https://academic.oup.com/nar/article/40/12/5215/2414626
    # 100bp length will be good for our purpose
    vienna_input = str(sequence)
    vienna_call = subprocess.run(vienna_args, stdout = subprocess.PIPE, stderr = subprocess.PIPE, input = vienna_input, encoding = 'ascii')
    
    output = vienna_call.stdout.strip().split('\n')
    errors = vienna_call.stderr

    return output, errors

def term_local_fold(sequence, term_struct, term_energy, antiterm_struct):
    polyU_regex = '[T]{3,}[ACGT]{,1}[T]{1,}[ACGT]{,1}[T]{1,}'
    #print(sequence) #debug
    
    #search for last ')'
    if pd.isna(antiterm_struct):
        return term_struct, term_energy, 'no_good_structure_found;no_aterm'
    antiterm_end = antiterm_struct.rfind(')')
    if antiterm_end < 0: #if ')' not found
        antiterm_end = 0
        #Possibly return here
    
    found_poly_u = " NO_POLY_U"
    offset = 10 # how many chars after antiterm end is start of poly U search
    poly_u = re.search(polyU_regex, sequence[(antiterm_end + offset):])
    # search for a poly U starting at the 10th nucleotide after the antiterminators end
    
    if(poly_u):
        search_start = poly_u.start() + offset + antiterm_end
        search_sequence = sequence[:search_start] # truncate sequence to poly U
        #print(sequence[:search_start], 'poly U search')
        found_poly_u = "" # flag that a poly U sequence was found
    else:
        search_sequence = sequence[:-5] # use entire sequence except last 5 bp
    
    #print(search_sequence)
    
    lfold_window = [100, 50, 40, 30, 20] # window sizes for the lfold algorithm
    for window in lfold_window:
        output, errors = local_fold(search_sequence, window)
        correct = False
        loops = [] # all structures from RNALfold get dumped here
    
        if(output):
            output_lines = output[::-1] # reverse order
            index = 2 # skip first line (last line in original output)
            maximum = len(output_lines)
            while(index < maximum):
                # parse all of the lines containing structures
                split_line = output_lines[index].split()
                lfold_structure, lfold_energy, lfold_start = split_line[0], split_line[-2], int(split_line[-1])
                
                if((lfold_start < (antiterm_end)) ):
                # We take any structures that start before the end of the antiterminator
                    brackets = lfold_structure.count('(') + lfold_structure.count(')')
                    loops.append([lfold_structure, lfold_start, lfold_energy, brackets])
                index += 1
        # reject any structures with multiple loops or that start before 5nt
        good_loops = []
        for element in loops:
            #print(element)
            if((re.search('[\)][\.]+[\(]', element[0]) == None) and (re.search('[\)][\(]', element[0]) == None)):
                if(element[1] >= 5): # must start at least after the fifth nucleotide
                    if(element[3] > 11): # at least 12 brackets in the structure
                        #print(sequence[(len(element[0]) + element[1] - 1):(len(element[0]) + element[1]) + 7])
                        if(re.search(polyU_regex, sequence[(len(element[0]) + element[1] - 1):(len(element[0]) + element[1]) + 7])):
                        # search for a poly T near the end of the loop
                            good_loops.append(element)
                            correct = True
                            #print('Good!')
        if(correct):
            break
    
    if(correct):
        best_loop = max(good_loops, key=lambda x: x[3]) # get loop with largest number of brackets
        lfold_start, lfold_structure = best_loop[1], best_loop[0]
        term_foldseq = sequence[lfold_start - 1:(lfold_start + len(lfold_structure) - 1)]
        new_term_structure = '.' * (lfold_start - 1) + lfold_structure + '.' * (len(sequence) - len(lfold_structure) - lfold_start + 1)
        
        energy, errors = get_energy(sequence, new_term_structure)
        return new_term_structure, energy, errors+found_poly_u

    # fallback option takes original term structure etc in case nothing works
    return term_struct, term_energy, 'no_good_structure_found'+found_poly_u

def run_thermo(tboxes):
    #Don't need: fasta_antiterm_start, fasta_antiterm_end
    #From RUN:
    tboxes[['whole_antiterm_structure', 'other_stems', 'whole_antiterm_warnings']] = tboxes.apply(lambda x: parse_structures(x['FASTA_sequence'], x['Tbox_start'], x['Tbox_end'], x['Sequence'], x['Structure']), axis = 'columns', result_type = 'expand')
    print('Structure correction complete.')
    #Get the terminator sequence as a string
    tboxes['term_sequence'] = tboxes.apply(lambda x: get_sequence(x['FASTA_sequence'], x['discrim_end'], x['term_end']), axis = 'columns', result_type = 'expand')
    print('Terminators found.')
    # fold the terminator sequence obtained above to get the secondary structure
    tboxes[['term_structure', 'terminator_energy', 'term_errors']] = tboxes.apply(lambda x: get_fold(x['term_sequence']), axis = 'columns', result_type = 'expand')
    print('Terminators folded.')
    # get the sequence from antiterm start to term end, will be used for comparing conformations with equal length
    tboxes['antiterm_term_sequence'] = tboxes.apply(lambda x: get_sequence(x['FASTA_sequence'],x['antiterm_start'] - 1, x['term_end']), axis = 'columns', result_type = 'expand')
    print('Antiterminator sequence found.')
    # get the antiterminator structure
    tboxes['infernal_antiterminator_structure'] = tboxes.apply(lambda x: get_sequence(x['whole_antiterm_structure'],x['antiterm_start'] - 1, x['term_end']), axis = 'columns', result_type = 'expand')
    print('Antiterminator INFERNAL structure found.')
    # get energy using RNAeval
    # this is DEPRECATED because we'll do constrained folding refinement
    #tboxes[['infernal_antiterminator_energy', 'infernal_antiterminator_errors']] = tboxes.apply(lambda x: get_energy(x['antiterm_term_sequence'], x['infernal_antiterminator_structure']), axis = 'columns', result_type = 'expand')
    print('Antiterminator energy done.')
    # refold the antiterminator using RNAfold with hard constraints
    tboxes[['vienna_antiterminator_structure', 'vienna_antiterminator_energy', 'vienna_antiterminator_errors']] = tboxes.apply(lambda x: make_antiterm_constraints(x['antiterm_term_sequence'], x['infernal_antiterminator_structure']), axis = 'columns', result_type = 'expand')
    print('Antiterminator re-folding done.')
    # get the terminator structure and energy, structure is term structure with dots in front so that structure
    # spans from antiterm_start to term_end
    tboxes['terminator_structure'] = tboxes.apply(lambda x: ('.' * 8 + x['term_structure']), axis = 'columns', result_type = 'expand')
    print('Terminator structure generated.')
    tboxes[['terminator_energy', 'terminator_errors']] = tboxes.apply(lambda x: get_energy(x['antiterm_term_sequence'], x['terminator_structure']), axis = 'columns', result_type = 'expand')
    print('Terminator energy calculated.')
    
    #Re-calculate the terminator structure
    tboxes[['new_term_structure', 'new_term_energy', 'new_term_errors']] = tboxes.apply(lambda x: term_local_fold(x.antiterm_term_sequence, x.terminator_structure, x.terminator_energy, x.vienna_antiterminator_structure), axis = 'columns', result_type = 'expand')
    print('Terminator refined')
    
    # make the tbox terminator structure. Use the refined structure.
    tboxes[['whole_term_structure', 'term_start']] = tboxes.apply(lambda x: get_whole_structs(x['antiterm_start'], x['whole_antiterm_structure'], x['new_term_structure']), axis = 'columns', result_type = 'expand')
    #make the Vienna antiterm structure for the whole T-box
    tboxes['folded_antiterm_structure'] = tboxes.apply(lambda x: get_whole_structs(x['antiterm_start'], x['whole_antiterm_structure'], x['vienna_antiterminator_structure'])[0], axis = 'columns', result_type = 'expand')
    print('Thermodynamic calculations complete.')
    return tboxes

#Generate trimmed structures and sequences
def trim(seq_df):
    seq_df['Trimmed_sequence']=""
    seq_df['Trimmed_antiterm_struct']=""
    seq_df['Trimmed_term_struct']=""

    counter = 0

    for i in range(0,len(seq_df)):
        print('Trimming ' + seq_df['Name'][i])
        seq = seq_df['FASTA_sequence'][i]
        t_struct = seq_df['whole_term_structure'][i]
        a_struct = seq_df['folded_antiterm_structure'][i]
    
        if pd.isna(seq):
            continue
        if pd.isna(t_struct) or pd.isna(a_struct) or pd.isna(seq_df['antiterm_end'][i]):
            #Fallback: only trim sequence
            try:
                tbox_start = int(seq_df['Tbox_start'][i])
                #Slice fasta sequence to get the tbox sequence
                seq = seq[tbox_start-1:]
                seq_df.at[seq_df.index[i], 'Trimmed_sequence'] = seq
                counter += 1
            except:
                pass
            continue
        
        tbox_start = int(seq_df['Tbox_start'][i])
        term_end = min(term_end_regex(seq,int(seq_df['antiterm_end'][i]) + 10),int(seq_df['term_end'][i]))
        seq_df.at[seq_df.index[i], 'term_end'] = term_end #update term_end
    
        if (tbox_start < -1) or term_end > len(seq): #note these are 1-indexed
            print('Error')
            continue
    
        #Slice fasta sequence to get the tbox sequence
        seq = seq[tbox_start-1:term_end]
    
        if len(seq) != len(a_struct):
            a_struct += '.'*(len(seq)-len(a_struct))
        
        if len(seq) != len(t_struct):
            t_struct += '.'*(len(seq)-len(t_struct))
    
        #Also slice the structures
        t_struct = t_struct[tbox_start-1:term_end]
        a_struct = a_struct[tbox_start-1:term_end]
    
        seq_df.at[seq_df.index[i], 'Trimmed_sequence'] = seq
        seq_df.at[seq_df.index[i], 'Trimmed_antiterm_struct'] = a_struct
        seq_df.at[seq_df.index[i], 'Trimmed_term_struct'] = t_struct
        counter += 1

    print("Trimmed sequences: " + str(counter))
    return seq_df

def wobblepair(base1, base2):
    base1 = base1.upper().replace('T','U')
    base2 = base2.upper().replace('T','U')
    if base1 == 'A':
        return base2 == 'U'
    if base1 == 'C':
        return base2 == 'G'
    if base1 == 'G':
        return base2 == 'C' or base2 == 'U'
    if base1 == 'U':
        return base2 == 'A' or base2 == 'G'
    if base1 == 'N' or base2 == 'N':
        return True
    return False

def clean(dot_structure, seq):
    if pd.isna(dot_structure) or pd.isna(seq):
        return None
    str_clean = ['.']*len(dot_structure)
    
    str_len = len(dot_structure)
    
    for i in range(str_len - 1):
        if dot_structure[i] == '(':
            l_count = 1 #counting the original '(' that we're looking for
            for j in range(i+1, str_len):
                if dot_structure[j] == '(':
                    l_count += 1
                if dot_structure[j] == ')':
                    l_count -= 1
                if l_count == 0: #We found the pair!
                    if wobblepair(seq[i],seq[j]): #But is it good?
                        str_clean[i] = '('
                        str_clean[j] = ')'
                    break
    return "".join(str_clean) #Convert to string from list


def clean_sequences(predseq):
    predseq["codon_region"]=predseq["codon_region"].str.upper()
    predseq["codon"]=predseq["codon"].str.upper()
    predseq["FASTA_sequence"]=predseq["FASTA_sequence"].str.upper()
    
    predseq[["Trimmed_antiterm_struct"]] = predseq.apply(lambda x: clean(x['Trimmed_antiterm_struct'], x['Trimmed_sequence']), axis = 'columns', result_type = 'expand')
    predseq[["Trimmed_term_struct"]] = predseq.apply(lambda x: clean(x['Trimmed_term_struct'], x['Trimmed_sequence']), axis = 'columns', result_type = 'expand')
    
    return predseq

def add_thermocalc(predseq):
    deltadeltaG = [None] * len(predseq)
    
    for i in range(0, len(predseq)):
        print('Thermocalculation for index - '+str(i))
        try:
            term_e=float(predseq['new_term_energy'].iloc[i])
            anti_e=float(predseq['vienna_antiterminator_energy'].iloc[i])
            ddG=float(term_e)-float(anti_e)
            deltadeltaG[i]=str(ddG)

        except (ValueError, TypeError):
            pass
    
    predseq['deltadelta_g']=deltadeltaG
    
    return predseq

#The main function to predict T-boxes
#For TRANSCRIPTIONAL T-boxes only (RF00230)
def tbox_predict(INFERNAL_file, predictions_file, fasta_file = None, score_cutoff = 15):
    score_cutoff = int(score_cutoff) #makes it an int, if it was passed as a string
    downstream_bases = 100 #for the fasta processing, number of downstream bases to include after antiterminator end
    
    #Read the input file into a dataframe
    tbox_all_DF = read_INFERNAL(INFERNAL_file)
    #Initialize the dataframe columns for prediction output
    tbox_all_DF['s1_start'] = -1
    tbox_all_DF['s1_loop_start'] = -1
    tbox_all_DF['s1_loop_end'] = -1
    tbox_all_DF['s1_end'] = -1
    tbox_all_DF['antiterm_start'] = -1
    tbox_all_DF['antiterm_end'] = -1
    tbox_all_DF['term_start'] = -1
    tbox_all_DF['term_end'] = -1
    tbox_all_DF['codon_start'] = -1
    tbox_all_DF['codon_end'] = -1
    tbox_all_DF['codon'] = ""
    tbox_all_DF['codon_region'] = ""
    tbox_all_DF['discrim_start'] = -1
    tbox_all_DF['discrim_end'] = -1
    tbox_all_DF['discriminator'] = ""
    tbox_all_DF['warnings'] = ""
    tbox_all_DF['type'] = "Transcriptional"
    tbox_all_DF['source'] = fasta_file.split('/')[-1]
        
    #Predict the t-boxes
    for i, name in enumerate(tbox_all_DF['Name']):
        
        #Predict the features. Use offset of 1 to convert 0-indexed Python format to standard sequence format
        tbox = tbox_features(tbox_all_DF['Sequence'][i], tbox_all_DF['Structure'][i], offset = 1)
        
        #Assign output to dataframe
        tbox_all_DF.at[tbox_all_DF.index[i], 's1_start'] = tbox[0]
        tbox_all_DF.at[tbox_all_DF.index[i], 's1_loop_start'] = tbox[1]
        tbox_all_DF.at[tbox_all_DF.index[i], 's1_loop_end'] = tbox[2]
        tbox_all_DF.at[tbox_all_DF.index[i], 'codon'] = tbox[3]
        tbox_all_DF.at[tbox_all_DF.index[i], 's1_end'] = tbox[4]
        tbox_all_DF.at[tbox_all_DF.index[i], 'antiterm_start'] = tbox[5]
        tbox_all_DF.at[tbox_all_DF.index[i], 'discriminator'] = tbox[6]
        tbox_all_DF.at[tbox_all_DF.index[i], 'antiterm_end'] = tbox[7]
        tbox_all_DF.at[tbox_all_DF.index[i], 'codon_region'] = tbox[8]
        tbox_all_DF.at[tbox_all_DF.index[i], 'warnings'] += tbox[9]
        tbox_all_DF.at[tbox_all_DF.index[i], 'discrim_start'] = tbox[10]
        tbox_all_DF.at[tbox_all_DF.index[i], 'discrim_end'] += tbox[11]
        
        #Check the score
        if tbox_all_DF.at[tbox_all_DF.index[i], 'Score'] < score_cutoff:
            tbox_all_DF.at[tbox_all_DF.index[i], 'warnings'] += "LOW_SCORE;"
    
    #Add the fasta sequences
    for i in range(len(tbox_all_DF['Name'])):
    
        with open(fasta_file) as f:

            for fasta in SeqIO.parse(f,'fasta'):
                if fasta.id == tbox_all_DF['Name'][i]:
                    if tbox_all_DF['Strand'][i] == '+':
                        tbox_all_DF.at[tbox_all_DF.index[i], 'FASTA_sequence'] = str(fasta.seq[tbox_all_DF['Tbox_start'][i] - 1:tbox_all_DF['Tbox_end'][i] + downstream_bases])
                    elif tbox_all_DF['Strand'][i] == '-':
                        tbox_all_DF.at[tbox_all_DF.index[i], 'FASTA_sequence'] = str(fasta.seq[tbox_all_DF['Tbox_end'][i] - 1:tbox_all_DF['Tbox_start'][i] + downstream_bases].reverse_complement())

    #Reset the T-box start and end
    tbox_all_DF['Name'] = tbox_all_DF['Name']+':'+tbox_all_DF['Tbox_start'].astype(str)+'-'+tbox_all_DF['Tbox_end'].astype(str)
    
    tbox_all_DF['Tbox_start'] = 1
    tbox_all_DF['Tbox_end'] = tbox_all_DF['FASTA_sequence'].str.len()
    
    tbox_all_DF['Tbox_end'] = tbox_all_DF['Tbox_end'] - downstream_bases #adjust for downstream bases

    #Convert positions from INFERNAL-relative to FASTA-relative
    merged = tbox_derive(tbox_all_DF)
    
    print('Feature derivation complete. Running thermodynamics.')
    thermo = run_thermo(merged)
    thermo = add_thermocalc(thermo)
    
    print('Removing duplicate T-boxes')
    thermo.drop_duplicates(subset = 'Sequence', keep = 'first', inplace = True) #Drop duplicate T-boxes
    
    print('Trimming structures and sequences')
    thermo.to_csv(predictions_file, index = False, header = True)
    thermo = trim(thermo)
    thermo = clean_sequences(thermo)
    
    #Write output
    thermo.to_csv(predictions_file, index = False, header = True)
    return 0

#Get arguments from command line and run the prediction
if len(sys.argv) > 5 or len(sys.argv) < 3:
    print("Error: incorrect number of arguments: %d" % len(sys.argv))
    print(sys.argv)
else:
    tbox_predict(*sys.argv[1:len(sys.argv)])
