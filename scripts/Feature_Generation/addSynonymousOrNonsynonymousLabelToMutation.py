# script to figure out if a MAPT mutation is synonymous or non-synonymous
# mutations need to have coordinated with respect to the exon-intron junction
# exonic mutations will have negative coordinates while intronic mutations will have 
# positive coordinates. 
# Script takes in two inputs: mutation coordinate file and a write file nameß

from Bio.Data import CodonTable
import math
import sys

mut_file = sys.argv[1]
write_file = sys.argv[2]

with open("../data/MAPT_exon10intron10withDMSdata.fa") as f:
    mapt_seq = [line.strip() for line in f][1]
    
mapt_seq_exon10 = mapt_seq[:93]


standard_table = CodonTable.unambiguous_rna_by_id[1]

with open(mut_file) as f:
    muts = [line.strip().split("\t") for line in f]


labels_for_mut = []
for mut in muts:
    if mut[0]=="WT":
        labelToAdd = "WT"
    else:
        mut_pos = int(mut[1])
        mut_val = mut[3]
        if mut_pos > 0:
            labelToAdd = "Intronic"
        else:
            # Index position of mutation w.r.t exon seq
            mut_pos_e = (93+mut_pos)+1
        
            # Figure out codon start position
            if (mut_pos_e%3) == 0:
                start_codon_pos = mut_pos_e - 2
                pos_on_codon = 2
            elif (mut_pos_e+1)%3 == 0:
                start_codon_pos = mut_pos_e - 1
                pos_on_codon = 1
            else:
                start_codon_pos = mut_pos_e
                pos_on_codon = 0
        
            current_codon = mapt_seq_exon10[start_codon_pos-1:start_codon_pos+2]
            changed_codon = current_codon[0:pos_on_codon] + mut_val + current_codon[pos_on_codon+1:]
        
            if changed_codon in standard_table.stop_codons:
                labelToAdd = "NonSynonymous"
            elif standard_table.forward_table[current_codon]==standard_table.forward_table[changed_codon]:
                labelToAdd = "Synonymous"
            else:
                labelToAdd = "NonSynonymous"
        
    mut_toAdd = mut[0:]
    mut_toAdd.append(labelToAdd)
    labels_for_mut.append(mut_toAdd)

with open(write_file,"w") as fw:
    for mutToWrite in labels_for_mut:
        fw.write("\t".join(mutToWrite))
        fw.write("\n")