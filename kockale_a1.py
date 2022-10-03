#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 21:22:27 2021

@author: logankocka
"""

#get the arguments from the command line call, save in list to reference later
import sys
    
# ------------ MODULES -------------------
def readFASTA(fileName, seqChoice):
    f = open(fileName, 'r')
    # allLines = f.readlines()
    # name is just accession number, description is everything else
    
    sequence = ""
    seq_name_list = []
    sequence_list =[]
    description_list = []
    
    if "fasta" in fileName:
        for line in f:
            if line[0] == ">":
                
                seq_name = line.split(" ",1)[0]
                seq_name = seq_name.strip('>')
                seq_name_list.append(seq_name)
                description = line.split(" ",1)[1]
                description = description.replace('\n','')
                description_list.append(description)
                if sequence:
                    sequence = sequence.replace("\n","")
                    sequence_list.append(sequence)
                    sequence = ''
            else:
                sequence += line
        
        sequence_list.append(sequence.replace("\n",""))
        sequence = sequence_list[seqChoice-1]
        seq_name = seq_name_list[seqChoice-1]
        description = description_list[seqChoice-1]
        
    # if it is a fastq file run this 
    if "fastq" in fileName:
        allLines = open(fileName, 'r').readlines()
        for i in range(0,len(allLines),4):
            seq = allLines[i:i+4]
            l = seq[:2]
            seq_name = l[0].split(" ")[0]
            description = l[0].split(" ")[1]
            sequence = l[1].replace("\n","")
   
    
    return seq_name,description,sequence

# screen 2
def printlnFASTA(seqName, seqDescription, sequence, n):
    # print the name and description
    print(">" + seqName + " " + seqDescription)
    # FASTA_LINE_NUM = n
    # for i in range(0,len(sequence),FASTA_LINE_NUM):
    #     print(sequence[i:i+FASTA_LINE_NUM])
    for i in range(0, len(sequence),n):
        print(sequence[i:i+n])
    return

# screens 3 and 4
def printWithRuler(seqName, sequence, seqDescription, spacerChoice, n):
    NON_FASTA_LINE_NUM = n
    if spacerChoice == "Y":
        spacer = " "
    if spacerChoice == "N":
        spacer = ""
    
    # print first line and numbers with spaces
    print('>' + seqName + ' ' + seqDescription + '\n')
    print(15*' ', end="")
    nums = list(range(1,11))
    for i in nums:
        print(nums[i-1], end="")
        if i < 9:
            print(9*' ' + spacer, end="")
        if i == 9:
            print(8*' ' + spacer, end="")
    print('\n' + " Line", end=" ")
    for i in nums:
        print("1234567890" + spacer, end="")
        if i == 10: print()
    
    # print each line with spacer every ten characters
    if spacerChoice == "N":
        # format the list
        NON_FASTA_LINE_NUM = 100 # split every 100 characters
        l_N = [sequence[i:i+NON_FASTA_LINE_NUM] for i in range(1, len(sequence), NON_FASTA_LINE_NUM)]
        #print them
        for count,val in enumerate(l_N):
            print('{:>5} {:<100}'.format(count+1,val))
    
    if spacerChoice == "Y":
        # format the list
        l_Y = [sequence[i:i+10] for i in range(1, len(sequence), 10)]
        l_Y2 = []; start = 0
        for i in l_Y:
            str = " ".join(l_Y[start:start+10])
            l_Y2.append(str)
            start += 10
        l_Y3 = list(filter(None, l_Y2))
        # print them
        for count,val in enumerate(l_Y3):
            print('{:>5} {:<110}'.format(count+1,val))

    return


def nucleotideCounter(sequence):
    # initialize
    counts = {'A': 0, 'G': 0, 'C': 0, 'T': 0}
  
    for i in sequence: # count values
        if i in counts:
            counts[i] += 1
        else:
            counts[i] = 1
    
    print("Nucleotide Counts: ") # print the key value pairs
    if 'N' in counts: print("A=[{A}] T=[{T}] G=[{G}] C=[{C}] N=[{N}]".format(**counts))
    else: print("A=[{A}] T=[{T}] G=[{G}] C=[{C}]".format(**counts))
    # some don't have N, some do...
    return counts


def gcContent(sequence):
    # find percentage of all nucleotides are G or C
    counts = {'A': 0, 'G': 0, 'C': 0, 'T': 0}
    for i in sequence: # count values
        if i in counts:
            counts[i] += 1
        else:
            counts[i] = 1
    total = sum(counts.values())
    
    GC_Counts = int(counts.get('G')) + int(counts.get('C'))
    
    return print("GC content: ", round(GC_Counts/total*100,2), "%")


def diCounter(sequence):
    # count freq of all nucleotide pairings
    counts = {'AA': 0, 'AT': 0, 'AG': 0, 'AC': 0, 'TA': 0, 'TT': 0, 'TG': 0, 'TC': 0,
              'GA': 0, 'GT': 0, 'GG': 0, 'GC': 0, 'CA': 0, 'CT': 0, 'CG': 0, 'CC': 0}
    for i in range(len(sequence)):
        if i != (len(sequence)-1):
            slice = sequence[i:i+2:1]
            if slice in counts:
                counts[slice] += 1
            else:
                counts[slice] = 1
   
    print("Di-nucleotide Counts: ")
    print("AA={AA} AT={AT} AG={AG} AC={AC}".format(**counts))
    print("TA={TA} TT={TT} TG={TG} TC={TC}".format(**counts))
    print("GA={GA} GT={GT} GG={GG} GC={GC}".format(**counts))
    print("CA={CA} CT={CT} CG={CG} CC={CC}".format(**counts))
    print()
    return counts


def inquiry(sequence, frag):
    start = int(frag.partition('::')[0])
    end = int(frag.partition('::')[2])
    
    #print fragment length
    print("The fragment you selected has a length of ", end-start, " nucleotides:")
    numDashes = (end-start) - 1 - (int(len(str(start)))) - (int(len(str(end))))
    print("<",start,"-"*numDashes,end,">", sep="")
    seqFrag = sequence[start:end+1]
    print("|"*len(seqFrag))
    print(seqFrag) # keep this
    nucleotideCounter(seqFrag)
    if 'G' in seqFrag and 'C' in seqFrag:
        gcContent(seqFrag)
    elif 'G' in seqFrag or 'C' in seqFrag:
        print("GC content: ", round((1/len(seqFrag)),2), "%")
    else:
        print("GC content: 0%")
    
    diCounter(seqFrag)
    
    return 



# ------- MAIN --------------------------------
def main():
    args = sys.argv
    fileName = args[1]
    user = args[0].split("_")[0]
    
    #print the welcome message
    print("\nWelcome Sequence Viewer! Programmer: ", user)
    #print sequences detected
    num = len([1 for line in open(fileName) if line.startswith(">")])
    if num == 1:
        print("There is 1 sequence detected in the file: ", fileName)
    else: 
        print("There are ", num, " sequences detected in the file: ", fileName)

    print("Which sequence do you want to examine ", [*range(1, num+1)], "?")
    seqChoice = int(input())

    print("\n[Part I]: Display Mode \n")
    
    # get viewing format choice
    while True:
        try:
            formatChoice = input("Do you want to view the sequence in FASTA format? (Y/N) ")
            if formatChoice == "Y" or formatChoice == "N":
                break;
            else: 
                print("Invalid. Enter Y or N.")
        except:
            continue

    # read the file
    new_lst = readFASTA(fileName, seqChoice)

    # if yes, print Screen 2
    if formatChoice == "Y":
        printlnFASTA(new_lst[0], new_lst[1], new_lst[2], 60)

    #if N, ask about spacer option
    if formatChoice == "N":
        while True:
            try:
                spacerChoice = input("Do you need a spacer for viewing nucleotide positions? (Y/N) ")
                if spacerChoice == "Y" or spacerChoice == "N":
                    break;
                else: 
                    print("Invalid. Enter Y or N.")
            except:
                continue
    
    if formatChoice == "N":
        # print with or without spacer
        printWithRuler(new_lst[0], new_lst[2], new_lst[1], spacerChoice, 100)
        
    print("\n[Part II]: Analysis Mode \n")
    nucleotideCounter(new_lst[2])
    gcContent(new_lst[2])
    diCounter(new_lst[2])
    
    print("[Part III]: Inquiry Mode \n")
    print("You can extract a DNA fragment from the given sequence.")
    
    #loop for getting the right kind of input for viewing a fragment
    while True:
        try:
            frag = input("Please enter the start and end positions (e.g., 19::48): ")
            # if frag == "N":
            #     sys.exit()
            if "::" not in frag:
                print("Separate start and end value with ""::""")
                continue
            start = int(frag.partition('::')[0])
            end = int(frag.partition('::')[2])
            # if a fragment is entered, call the function 
        except ValueError:
            print("Start and end values must be integers.")
            continue
        if start < 0:
            print("Start value must be nonnegative.")
            continue
        elif end > len(new_lst[2]):
            print("End value can't be larger than sequence length.")
            continue
        elif start-end > 0:
            print("End value must be higher than start value.")
            continue
        else:
            inquiry(new_lst[2], frag)
            another = input("Do you want to extract another fragment? ")
            if another == "N":
                sys.exit()
            continue
    
    return


if __name__=='__main__':
    main()
    
    
    