#!/usr/bin/env python

# Import required modules
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import regex
from itertools import product
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from lru import LRU

## Store primers as string
#cyst_ITSF = "CTGCTGCTGGATCATTACCC"
#cyst_ITSR = "CGCCGGCTGTCTCTCCAA"

# Convert reverse primer sequence string to seq object,
# then convert to reverse complement and store this as string
cyst_ITSR_seq = Seq(cyst_ITSR, generic_dna)
cyst_ITSR_rc = str(cyst_ITSR_seq.reverse_complement())

# Compile primers as regular expressions allowing
# a positional mismatch (substitution) of one or less
cyst_ITSF = regex.compile("(CTGCTGCTGGATCATTACCC){s<=1}")
cyst_ITSR = regex.compile("(CGCCAGCACAGCCGTTAG){s<=1}")

# Store list of barcode sequences as variable

barcodes_AF = ["CTACCT", "GGATCT", "GTCGTA", "GAGTCA", "CCTTCT", "GACTTC",
               "GGAGAA", "TCGACT"]

barcodes_BF = ["GTAACC", "GATGAC", "CCATAC", "CTCACA", "GAACGT", "CCGTTA",
               "ACCTCA", "CAACTC"]

barcodes_AR = ["CTACCT", "GGATCT", "GTCGTA", "GAGTCA", "CCTTCT", "GACTTC",
               "GGAGAA", "TCGACT", "GTAACC", "GATGAC", "CCATAC", "CTCACA"]

barcodes_BR = ["GAACGT", "CCGTTA", "ACCTCA", "CAACTC", "CTTGGT", "GTTCCA",
               "CTGAAC", "AGGACA", "CGTCAT", "TTCCGT", "ACCGAT", "TGTGTC"]

# Create an empty list to store reverse complement of barcodes
Abar_rev = []
Bbar_rev = []

# Using Biopython generate a list of reverse complement barcodes
for bar in barcodes_AR:
    # Convert string to seq object
    dna_bar = Seq(bar, generic_dna)
    # Use .reverse_complement function to genereate
    # reverse complement of seq object
    rcbarseq = dna_bar.reverse_complement()
    # Convert reverse complement seq objects back to string
    rcbarstring = str(rcbarseq)
    # Append rctags as string to list
    Abar_rev.append(rcbarstring)

for bar in barcodes_BR:
    # Convert string to seq object
    dna_bar = Seq(bar, generic_dna)
    # Use .reverse_complement function to genereate
    # reverse complement of seq object
    rcbarseq = dna_bar.reverse_complement()
    # Convert reverse complement seq objects back to string
    rcbarstring = str(rcbarseq)
    # Append rctags as string to list
    Bbar_rev.append(rcbarstring)

# Generate a list of all forward and reverse barcode combinations
A_barcode_combos = list(product(barcodes_AF, Abar_rev))
B_barcode_combos = list(product(barcodes_BF, Bbar_rev))

# Start with A barcodes -------------------------------------------------------

# Define sample file handle
input_handle = open("all-cyst-amp-maxeefiltered.fastq")
# Open all the output handles at the start, and keep them in a dictionary
cache = LRU(999)

filenames = {}
for x, (fbar, rbar) in enumerate(A_barcode_combos):
    f = "cyst-%05iA-%s-%s.fastq" % (x, fbar, rbar)
    filenames[(fbar, rbar)] = f
    outfile = open(f, "w")
    outfile.close()

oddities = open('cyst-odditiesA.txt', 'w')

fragmentary = 0
fr_tallies = dict()

# Using SimpleFastaParser search each sequence for the 6nt seq before and
# after primers
for title, seq, qual in FastqGeneralIterator(input_handle):
    seqlen = len(seq)
    fprimersear = cyst_ITSF.search(seq)
    rprimersear = cyst_ITSR.search(seq)
    fstart = fprimersear.start()
    rstart = rprimersear.start()
    fend = fprimersear.end()
    rend = rprimersear.end()
    frame = seq[fstart:rend]
    fbar = seq[fstart-6:fstart]
    rbar = seq[rend:rend + 6]
    # Store the amplified region between primers
    noprimframe = seq[fend:rstart]
    noprimqual = qual[fend:rstart]
    # Hoping to find pair of 6bp known barcodes.
    #
    # Checking against the expected pairs via the dictionary ensures
    # will only write this out once, without needing a for loop.
    #
    # Might have only partial sequences, e.g. 'CTGA' and 'GG'
    # Might have pcr or sequencing erors returning barcodes not on our list
    if len(fbar) != 6 or len(rbar) != 6:
        print("Ignoring %s %s" % (fbar, rbar))
        fragmentary += 1
    else:
        # Right length, first count the barcode pair
        if (fbar, rbar) in fr_tallies:
            fr_tallies[(fbar, rbar)] += 1
        else:
            # First time to see it, count it
            fr_tallies[(fbar, rbar)] = 1
        if (fbar, rbar) in filenames and (fbar, rbar) in cache:
            print("Wanted  %s %s" % (fbar, rbar))
            cache[(fbar, rbar)].write("@%s\n%s\n+\n%s\n" %
                                      (title, noprimframe, noprimqual))
        elif (fbar, rbar) in filenames:
            name = filenames[(fbar, rbar)]
            cache[(fbar, rbar)] = open(name, "a")
            cache[(fbar, rbar)].write("@%s\n%s\n+\n%s\n" %
                                      (title, noprimframe, noprimqual))
        else:
            print("Unexpected %s %s" % (fbar, rbar))
            oddities.write("%s\t%s\t%s\t%s\n" % (fbar, rbar, title, seq))
        # Do we already have an output file ready for this pair?
        # There is a limit to how many files we can open at once...
        # if (fbar, rbar) not in out_handles:
        #    out_handles[(fbar, rbar)] =
        #    open("Pas-" + fbar + "-" + rbar + ".fasta", "w")
        # out_handles[(fbar, rbar)].write(">%s\n%s\n" % (title, seq))

oddities.close()
input_handle.close()

# Print the tallies for each expected barcode pair
print("Observed barcode pairs, tally count, wanted or not?")
for (fbar, rbar) in fr_tallies:
    print("%s %s count %i %r" %
          (fbar, rbar, fr_tallies[(fbar, rbar)],
           (fbar, rbar) in fr_combo_final))

print("In total %i full length barcode pairs, and %i fragments" %
      (sum(fr_tallies.values()), fragmentary))


# Start again for B barcodes --------------------------------------------------
input_handle = open("all-cyst-amp-maxeefiltered.fastq")
oddities = open('cyst-odditiesA.txt', 'w')
cache = LRU(999)

filenames = {}
for x, (fbar, rbar) in enumerate(B_barcode_combos):
    f = "cyst-%05iB-%s-%s.fastq" % (x, fbar, rbar)
    filenames[(fbar, rbar)] = f
    outfile = open(f, "w")
    outfile.close()

oddities = open('cyst-odditiesB.txt', 'w')

fragmentary = 0
fr_tallies = dict()

# Using SimpleFastaParser search each sequence
# for the 6nt seq before and after primers
for title, seq, qual in FastqGeneralIterator(input_handle):
    seqlen = len(seq)
    fprimersear = cyst_ITSF.search(seq)
    rprimersear = cyst_ITSR.search(seq)
    fstart = fprimersear.start()
    rstart = rprimersear.start()
    fend = fprimersear.end()
    rend = rprimersear.end()
    frame = seq[fstart:rend]
    fbar = seq[fstart-6:fstart]
    rbar = seq[rend:rend + 6]
    # Store the amplified region between primers
    noprimframe = seq[fend:rstart]
    noprimqual = qual[fend:rstart]
    # Hoping to find pair of 6bp known barcodes.
    #
    # Checking against the expected pairs via the dictionary ensures
    # will only write this out once, without needing a for loop.
    #
    # Might have only partial sequences, e.g. 'CTGA' and 'GG'
    # Might have pcr or sequencing erors returning barcodes not on our list
    if len(fbar) != 6 or len(rbar) != 6:
        print("Ignoring %s %s" % (fbar, rbar))
        fragmentary += 1
    else:
        # Right length, first count the barcode pair
        if (fbar, rbar) in fr_tallies:
            fr_tallies[(fbar, rbar)] += 1
        else:
            # First time to see it, count it
            fr_tallies[(fbar, rbar)] = 1
        if (fbar, rbar) in filenames and (fbar, rbar) in cache:
            print("Wanted  %s %s" % (fbar, rbar))
            cache[(fbar, rbar)].write("@%s\n%s\n+\n%s\n" %
                                      (title, noprimframe, noprimqual))
        elif (fbar, rbar) in filenames:
            name = filenames[(fbar, rbar)]
            cache[(fbar, rbar)] = open(name, "a")
            cache[(fbar, rbar)].write("@%s\n%s\n+\n%s\n" %
                                      (title, noprimframe, noprimqual))
        else:
            print("Unexpected %s %s" % (fbar, rbar))
            oddities.write("%s\t%s\t%s\t%s\n" % (fbar, rbar, title, seq))
        # Do we already have an output file ready for this pair?
        # There is a limit to how many files we can open at once...
        # if (fbar, rbar) not in out_handles:
        #    out_handles[(fbar, rbar)] =
        #    open("Pas-" + fbar + "-" + rbar + ".fasta", "w")
        # out_handles[(fbar, rbar)].write(">%s\n%s\n" % (title, seq))

oddities.close()
input_handle.close()

# Print the tallies for each expected barcode pair
print("Observed barcode pairs, tally count, wanted or not?")
for (fbar, rbar) in fr_tallies:
    print("%s %s count %i %r" %
          (fbar, rbar, fr_tallies[(fbar, rbar)],
           (fbar, rbar) in fr_combo_final))

print("In total %i full length barcode pairs, and %i fragments" %
      (sum(fr_tallies.values()), fragmentary))
