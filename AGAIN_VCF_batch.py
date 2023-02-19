# python3.8
__author__ = "Peng Zhang"
__copyright__ = "Copyright 2023, " \
                "Laboratory of Human Genetics of Infectious Diseases, " \
                "The Rockefeller University"
__license__ = "CC BY-NC-ND 4.0"
__version__ = "02-01-2023"

import os
import re
import time
import argparse

print('********************************************')
print('       #    ###        #   #####   #    #   ')
print('      ##   #          ##     #     ##   #   ')
print('     # #   #  ##     # #     #     #  # #   ')
print('    ####   #   #    ####     #     #   ##   ')
print('   #   #    ###    #   #   #####   #    #   ')
print(' Detecting human intronic AG-gain variants  ')
print('      in the region between BP and ACC      ')
print('*********************************************\n')

###
# input parameters
###

parser = argparse.ArgumentParser(description="AGAIN")
parser.add_argument("-d", "--dir", help="directory of vcf files")
parser.add_argument("-s", "--sample", help="sample list")
parser.add_argument("-o", "--output", help="output filename, tab-delimited")
parser.add_argument("-g", "--genome_ref", type=str, default='GRCh37', choices=['GRCh37', 'GRCh38'],
                    help="human reference genome assembly")
parser.add_argument("-t", "--transcript_ref", type=str, default='all', choices=['all', 'canonical'],
                    help="all/canonical transcripts")

args = parser.parse_args()
directory = args.dir
filename_sample = args.sample
filename_output = args.output
genome_ref = args.genome_ref
transcript_ref = args.transcript_ref

filename_coding_map = 'Data_AGAIN_'+genome_ref+'_coding_map.txt'
if transcript_ref == 'all':
    filename_AGAIN_detection = 'Data_AGAIN_'+genome_ref+'_detection_all.bed'
else:
    filename_AGAIN_detection = 'Data_AGAIN_'+genome_ref+'_detection_canonical.bed'

file_out = open(filename_output, 'w')
file_out.write('SAMPLE\tCHR\tPOS\tID\tREF\tALT\tSTR\tTYPE\tGENE\tTRANSCRIPT_IVS\tCANONICAL\t'
               'AGAIN_ZONE\tAGAIN_YAG\tAGAIN_BP_DIST\tAGAIN_ACC_DIST\tAGAIN_HIGHRISK\tAGAIN_SCORE\t'
               'PROT_SEQ_WT\tPROT_SEQ_NEW_ACC\tHGVS_NEW_ACC\tPROT_SEQ_EXON_SKIP\tHGVS_EXON_SKIP\n')

###
# global variables, and load coding map
###

BASE_PAIRING = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}

GENTIC_CODE = {
    'AAA':'K', 'AAC':'N', 'AAG':'K', 'AAT':'N', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AGA':'R', 'AGC':'S', 'AGG':'R', 'AGT':'S', 'ATA':'I', 'ATC':'I', 'ATG':'M', 'ATT':'I',
    'CAA':'Q', 'CAC':'H', 'CAG':'Q', 'CAT':'H', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'GAA':'E', 'GAC':'D', 'GAG':'E', 'GAT':'D', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'TAA':'*', 'TAC':'Y', 'TAG':'*', 'TAT':'Y', 'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TGA':'*', 'TGC':'C', 'TGG':'W', 'TGT':'C', 'TTA':'L', 'TTC':'F', 'TTG':'L', 'TTT':'F',}

AA_CODE = {
    'A':'Ala', 'C':'Cys', 'D':'Asp', 'E':'Glu', 'F':'Phe', 'G':'Gly', 'H':'His', 'I':'Ile',
    'K':'Lys', 'L':'Leu', 'M':'Met', 'N':'Asn', 'P':'Pro', 'Q':'Gln', 'R':'Arg', 'S':'Ser',
    'T':'Thr', 'V':'Val', 'W':'Trp', 'Y':'Tyr', '*':'*'}

file_coding_map = open(filename_coding_map, 'r')
file_coding_map.readline()
coding_map_transcript_start_codon = dict()
coding_map_transcript_start_exon = dict()
coding_map_transcript_exon_count = dict()
coding_map_transcript_prot_seq = dict()
coding_map_transcript_exon_region_seq = dict()

for eachline in file_coding_map:
    item = eachline.strip().split('\t')
    transcript = item[1]
    start_codon = int(item[4])
    start_exon = int(item[5])
    exon_count = int(item[6])
    exon_region_seq = item[7]
    exon_region_seq_list = exon_region_seq.split(',')
    prot_seq = item[8]

    coding_map_transcript_start_codon[transcript] = start_codon
    coding_map_transcript_start_exon[transcript] = start_exon
    coding_map_transcript_exon_count[transcript] = exon_count
    coding_map_transcript_prot_seq[transcript] = prot_seq
    coding_map_transcript_exon_region_seq[transcript] = dict()
    for exon in range(1, exon_count+1):
        coding_map_transcript_exon_region_seq[transcript][exon] = exon_region_seq_list[exon-1]


###
# AGAIN functions
###

# function to reverse sequence
def rev_seq(fwd_seq):
    output_seq = ''
    for fwd_pos in range(0, len(fwd_seq)):
        output_seq += BASE_PAIRING[fwd_seq[len(fwd_seq) - fwd_pos - 1]]
    return output_seq


# translation nucleotide sequence to protein sequence
def translation(nucl_seq):
    prot_seq = ''
    codon_list = [nucl_seq[i:i+3] for i in range(0, len(nucl_seq), 3)]
    for codon in codon_list:
        if len(codon) != 3:
            prot_seq += '?'
            break
        if 'N' in codon:
            prot_seq = 'N'
            break
        if GENTIC_CODE[codon] == '*':
            prot_seq += '*'
            break
        prot_seq += GENTIC_CODE[codon]
    return prot_seq


# new acceptor site in intron: protein sequence & HGVS
def AGAIN_prot_new_acc(ivs, strand, start_codon, start_exon, exon_count,
                       exon_region_seq_dict, prot_seq_wt, ins_seq):
    next_exon = ivs + 1

    # AGAIN variant before 5'UTR exon (no change to protein sequence)
    if next_exon < start_exon:
        prot_seq_mt = prot_seq_wt
        hgvs = '5UTR_exon_extended:no_change'

    # AGAIN variant before START exon
    elif next_exon == start_exon:
        start_exon_item = exon_region_seq_dict[start_exon].split('-')
        start_exon_start = int(start_exon_item[0])
        start_exon_end = int(start_exon_item[1])
        start_exon_seq = start_exon_item[2]

        if strand == '+':
            before_start_seq = ins_seq + start_exon_seq[0: start_codon - start_exon_start]
        else:
            before_start_seq = ins_seq + start_exon_seq[0: start_exon_end - start_codon]

        # no START within before_start_seq (no change to protein sequence)
        if 'ATG' not in before_start_seq:
            prot_seq_mt = prot_seq_wt
            hgvs = 'START_exon_extended:no_change'

        # new START within before_start_seq
        else:
            new_start_seq = before_start_seq[before_start_seq.find('ATG'):]
            new_start_prot_seq = translation(new_start_seq)

            # new STOP within new_start_seq (like uORF added before START)
            if '*' in new_start_prot_seq:
                prot_seq_mt = '(' + new_start_prot_seq + ') ' + prot_seq_wt
                hgvs = 'uORF_before_START:no_change'

            # no STOP within new_start_seq
            else:
                # new_start_seq is 3x (in-frame insertion, extension)
                if len(new_start_seq) % 3 == 0:
                    prot_seq_mt = new_start_prot_seq + prot_seq_wt
                    hgvs = 'p.Met1ext-' + str(len(new_start_prot_seq))

                # new_start_seq is NOT 3x (frameshift before START)
                else:
                    prot_seq_mt = '.'
                    hgvs = 'fs_before_START:total_change'

    # AGAIN variant after START exon
    else:
        # join all coding exons, from START codon to LAST exon, with the inserted seq
        coding_seq = ''
        for i in range(1, exon_count+1):
            exon_region_seq_item = exon_region_seq_dict[i].split('-')
            exon_start = int(exon_region_seq_item[0])
            exon_end = int(exon_region_seq_item[1])
            exon_seq = exon_region_seq_item[2]

            if i < start_exon:
                pass
            elif i == start_exon:
                if strand == '+':
                    coding_seq += exon_seq[start_codon - exon_start: ]
                else:
                    coding_seq += exon_seq[exon_end - start_codon: ]
            elif i == next_exon:
                coding_seq += ins_seq + exon_seq
            else:
                coding_seq += exon_seq

        # translate the coding sequence after adding ins_seq, and ins_seq alone
        prot_seq_mt = translation(coding_seq)
        ins_prot_seq = translation(ins_seq)
        hgvs = 'p.'

        # position of the first mismatch in protein sequence, wt vs mt
        j = 0
        while (len(prot_seq_wt[j:]) > 0) and (prot_seq_wt[j] == prot_seq_mt[j]):
            j += 1

        # the first mismatch is a STOP (stop-gain)
        if prot_seq_mt[j] == '*':
            hgvs += AA_CODE[prot_seq_wt[j]] + str(j+1) + '*'

        # inframe-insertion, with STOP in ins_seq
        elif (len(ins_seq) % 3 == 0) and ('*' in ins_prot_seq):
            ins_prot_seq_THREE = ''
            for each_AA in ins_prot_seq:
                ins_prot_seq_THREE += AA_CODE[each_AA]
            hgvs += AA_CODE[prot_seq_wt[j-1]] + str(j) + '_' + \
                    AA_CODE[prot_seq_wt[j]] + str(j+1) + 'ins' + ins_prot_seq_THREE

        # inframe-insertion, without STOP
        elif len(ins_seq) % 3 == 0:
            # remove STOP (*) from prot_seq_mt, for scanning the sequence inserted in mt
            prot_seq_mt = prot_seq_mt[0:-1]
            k = j
            while (len(prot_seq_mt[k:]) > 0) and (prot_seq_wt[j:] != prot_seq_mt[k:]):
                k += 1
            hgvs += AA_CODE[prot_seq_wt[j-1]] + str(j) + '_' + \
                    AA_CODE[prot_seq_wt[j]] + str(j+1) + 'ins' + str(k-j)
            # add STOP (*) back to prot_seq_mt
            prot_seq_mt += '*'

        # frameshift-insertion, with premature STOP, the fs_length includes STOP
        elif '*' in prot_seq_mt:
            fs_length = str(len(prot_seq_mt[j:]))
            hgvs += AA_CODE[prot_seq_wt[j]] + str(j+1) + \
                    AA_CODE[prot_seq_mt[j]] + 'fs*' + fs_length

        # frameshift-insertion, without premature STOP
        else:
            hgvs += AA_CODE[prot_seq_wt[j]] + str(j+1) + \
                    AA_CODE[prot_seq_mt[j]] + 'fs*?'

    return prot_seq_mt, hgvs


# complete skipping of the next exon: protein sequence & HGVS
def AGAIN_prot_exon_skip(ivs, strand, start_codon, start_exon, exon_count,
                         exon_region_seq_dict, prot_seq_wt):
    skipped_exon = ivs + 1

    # AGAIN variant before 5'UTR exon (no change to protein sequence)
    if skipped_exon < start_exon:
        prot_seq_mt = prot_seq_wt
        hgvs = '5UTR_exon_skipped:no_change'

    # AGAIN variant before START exon (totally different protein sequence)
    elif skipped_exon == start_exon:
        prot_seq_mt = '.'
        hgvs = 'START_exon_skipped:total_change'

    # AGAIN variant before LAST exon (exon not skipped, no change to protein sequence)
    elif skipped_exon == exon_count:
        prot_seq_mt = prot_seq_wt
        hgvs = 'LAST_exon_not_skipped:no_change'

    # AGAIN variant between START and LAST exon
    else:
        # join all coding exons, from START codon to LAST exon, skip variant's next exon
        coding_seq = ''
        for i in range(1, exon_count+1):
            exon_region_seq_item = exon_region_seq_dict[i].split('-')
            exon_start = int(exon_region_seq_item[0])
            exon_end = int(exon_region_seq_item[1])
            exon_seq = exon_region_seq_item[2]
            if i < start_exon:
                pass
            elif i == start_exon:
                if strand == '+':
                    coding_seq += exon_seq[start_codon - exon_start: ]
                else:
                    coding_seq += exon_seq[exon_end - start_codon: ]
            elif (i > start_exon) and (i != skipped_exon):
                coding_seq += exon_seq
            elif (i > start_exon) and (i == skipped_exon):
                pass

        # translate the coding sequence after exon skipping
        prot_seq_mt = translation(coding_seq)
        hgvs = 'p.'

        # position of the first mismatch in protein sequence, wt vs mt
        j = 0
        while (len(prot_seq_wt[j:]) > 0) and (prot_seq_wt[j] == prot_seq_mt[j]):
            j += 1

        # check the length of skipped exon, and determine its consequences and hgvs
        skipped_exon_seq = exon_region_seq_dict[skipped_exon].split('-')[2]

        # skipped exon is 3x (in-frame deletion)
        if len(skipped_exon_seq) % 3 == 0:
            # remove STOP (*) from prot_seq_mt, for scanning the sequence inserted in mt
            prot_seq_mt = prot_seq_mt[0:-1]
            k = j
            while (len(prot_seq_wt[k:]) > 0) and (prot_seq_wt[k:] != prot_seq_mt[j:]):
                k += 1
            hgvs += AA_CODE[prot_seq_wt[j]] + str(j+1) + '_' + \
                    AA_CODE[prot_seq_wt[k-1]] + str(k) + 'del'
            # add STOP (*) back to prot_seq_mt
            prot_seq_mt += '*'

        # skipped exon is NOT 3x (frameshift)
        else:
            # the first mismatch is a STOP (stop-gain)
            if prot_seq_mt[j] == '*':
                hgvs += AA_CODE[prot_seq_wt[j]] + str(j+1) + '*'

            # found premature STOP in prot_seq_mt, the fs_length includes STOP
            elif '*' in prot_seq_mt[j:]:
                fs_length = str(len(prot_seq_mt[j:]))
                hgvs += AA_CODE[prot_seq_wt[j]] + str(j+1) + \
                        AA_CODE[prot_seq_mt[j]] + 'fs*' + fs_length

            # not-found premature STOP in prot_seq_mt
            else:
                hgvs += AA_CODE[prot_seq_wt[j]] + str(j+1) + \
                        AA_CODE[prot_seq_mt[j]] + 'fs*?'

    return prot_seq_mt, hgvs


###
# AGAIN MAIN
###

file_sample = open(filename_sample, 'r')
for eachsample in file_sample:
    sample = eachsample.strip()
    filename_vcf = sample + '.vcf'
    filename_bed = sample + '.bed'
    filename_mapping = sample + '.mapping'
    start_time = time.time()
    print(sample)

    try:
        ###
        # read variants, convert variants from VCF to BED
        ###
        file_var_vcf = open(directory + filename_vcf, 'r')
        file_var_bed = open(filename_bed, 'w')
        input_var_count = 0

        for eachline in file_var_vcf:
            if not eachline.startswith('#'):
                input_var_count += 1
                item = eachline.strip().split('\t')
                chrom = item[0] if 'chr' in item[0] else 'chr' + item[0]
                pos = item[1]
                var_id = item[2]
                ref = item[3]
                alt = item[4]
                var_info = chrom+'*'+pos+'*'+var_id+'*'+ref+'*'+alt
                start = end = var_type = '.'

                if len(ref) == 1 and len(alt) == 1:
                    var_type = 'snv'
                    start = str(int(pos) - 1)
                    end = pos

                elif len(ref) == 1 and len(alt) > 1:
                    var_type = 'ins-' + str(len(alt) - 1) + 'nt'
                    start = pos
                    end = str(int(pos) + len(alt) - 1)

                elif len(ref) > 1 and len(alt) == 1:
                    var_type = 'del-' + str(len(ref) - 1) + 'nt'
                    start = pos
                    end = str(int(pos) + len(ref) - 1)

                if (start != '.') and (end != '.') and (var_type != '.'):
                    file_var_bed.write(chrom+'\t'+start+'\t'+end+'\t'+var_info+'\t.\t+\t'+var_type+'\n')
                    file_var_bed.write(chrom+'\t'+start+'\t'+end+'\t'+var_info+'\t.\t-\t'+var_type+'\n')

        file_var_vcf.close()
        file_var_bed.close()


        ###
        # mapping variants to BP-ACC region
        # identify AG-gain variant with AGAIN annotations
        # generate consequent protein sequence and hgvs
        ###
        os.system("bedtools intersect -wo -s"
                  " -a " + filename_bed +
                  " -b " + filename_AGAIN_detection +
                  " > " + filename_mapping)

        file_mapping = open(filename_mapping, 'r')
        AGAIN_var_set = set()
        AGAIN_var_highrisk_count = 0

        for eachline in file_mapping:
            try:
                item = eachline.strip().split('\t')
                chrom = item[0]
                var_start = int(item[1])
                var_end = int(item[2])
                var_info = item[3]
                strand = item[5]
                var_type = item[6]
                bp_acc_start = int(item[8])  # (+): BP-1, (-): ACC-1
                bp_acc_end = int(item[9])    # (+): ACC,  (-): BP
                gene = item[13]
                bp_name = item[14]
                bp_acc_dist = item[15]
                bp_order = item[16]
                bp_total = item[17]
                first_bp_pos = int(item[18])
                transcript_ivs = item[19]
                canonical = item[20]
                bp_acc_seq_wt = item[21]
                exist_AG = item[22]
                overlap = item[23]

                var_info_item = var_info.split('*')
                var_pos = int(var_info_item[1])
                var_id = var_info_item[2]
                ref = var_info_item[3]
                alt = var_info_item[4]

                nucleotide_list = list(bp_acc_seq_wt)
                variant_pass1 = False
                variant_pass2 = False

                AGAIN_zone = '.'
                # variants in first_BP-ACC
                if bp_order == 'FIRST':
                    AGAIN_zone = 'ZONE1'
                    variant_pass1 = True
                # variants in second_BP-ACC, only keep those within second_BP-first_BP
                if bp_order == 'SECOND':
                    AGAIN_zone = 'ZONE2'
                    if strand == '+' and (var_pos < first_bp_pos):
                        variant_pass1 = True
                    if strand == '-' and (var_pos > first_bp_pos):
                        variant_pass1 = True

                if variant_pass1:
                    # snv
                    if var_type == 'snv':
                        variant_pass2 = True
                        if strand == '+':
                            var_index = var_pos - bp_acc_start - 1
                            nucleotide_list[var_index] = alt
                        else:
                            var_index = bp_acc_end - var_pos
                            nucleotide_list[var_index] = BASE_PAIRING[alt]

                    # insertion (length < 20)
                    elif ('ins' in var_type) and (len(alt) < 20):
                        variant_pass2 = True
                        if strand == '+':
                            var_index = var_pos - bp_acc_start - 1
                            nucleotide_list[var_index] = alt
                        else:
                            var_index = bp_acc_end - var_pos
                            nucleotide_list[var_index] = rev_seq(alt)

                    # deletion (within BP-ACC region)
                    elif ('del' in var_type) and (bp_acc_end - len(ref) > var_pos > bp_acc_start):
                        variant_pass2 = True
                        if strand == '+':
                            var_index = var_pos - bp_acc_start - 1
                            for temp_index in range(var_index + 1, var_index + len(ref)):
                                nucleotide_list[temp_index] = ''
                        else:
                            var_index = bp_acc_end - var_pos
                            for temp_index in range(var_index - len(ref) + 1, var_index):
                                nucleotide_list[temp_index] = ''

                bp_acc_seq_mt = '.'
                AGAIN_YAG = 'NO'
                AGAIN_bp_dist = '.'
                AGAIN_acc_dist = '.'
                AGAIN_highrisk = 'NO'
                AGAIN_score = 1

                if variant_pass2:
                    bp_acc_seq_mt = ''.join(nucleotide_list)

                    # screen [BP+1, ACC-4] for AG
                    # output end+1 position of AG, e.g., (XXXAG):5
                    bp_acc_seq_wt_check = bp_acc_seq_wt[1:-3]
                    AG_hit_wt = set(hit.end() for hit in re.finditer('AG', bp_acc_seq_wt_check))

                    bp_acc_seq_mt_check = bp_acc_seq_mt[1:-3]
                    AG_hit_mt = set(hit.end() for hit in re.finditer('AG', bp_acc_seq_mt_check))

                    AGAIN_hit_set = AG_hit_mt - AG_hit_wt

                    # for insertion, shift AGAIN_hit by ins-size to check with pre-exisited AG
                    if 'ins' in var_type:
                        AGAIN_hit_set_new = set()
                        for each_pos in AGAIN_hit_set:
                            if (each_pos-len(alt)+1) not in AG_hit_wt:
                                AGAIN_hit_set_new.add(each_pos)
                        AGAIN_hit_set = AGAIN_hit_set_new

                    # AGAIN output
                    if AGAIN_hit_set and (var_info not in AGAIN_var_set):
                        AGAIN_hit_list = list(AGAIN_hit_set)
                        AGAIN_hit_list.sort()
                        AGAIN_hit = AGAIN_hit_list[-1]
                        AGAIN_bp_dist = AGAIN_hit
                        AGAIN_acc_dist = AGAIN_bp_dist - len(bp_acc_seq_mt) + 1
                        if bp_acc_seq_mt[AGAIN_hit - 2] in ['C', 'T']:
                            AGAIN_YAG = 'YES'

                        if AGAIN_zone == 'ZONE1':
                            AGAIN_score += 1
                        if AGAIN_YAG == 'YES':
                            AGAIN_score += 1
                        if AGAIN_bp_dist >= 8:
                            AGAIN_score += 1
                        if exist_AG == 'NO':
                            AGAIN_score += 1
                        if AGAIN_zone == 'ZONE1' and AGAIN_bp_dist >= 8:
                            AGAIN_highrisk = 'YES'
                            AGAIN_var_highrisk_count += 1

                        # AGAIN_prot
                        AGAIN_acc_seq = bp_acc_seq_mt[AGAIN_acc_dist:]
                        transcript_ivs_list = transcript_ivs.split(';')
                        transcript_ivs_item = transcript_ivs_list[0].split('_IVS')  # the first transcript
                        transcript = transcript_ivs_item[0]
                        ivs = int(transcript_ivs_item[1])

                        start_codon = int(coding_map_transcript_start_codon[transcript])
                        start_exon = int(coding_map_transcript_start_exon[transcript])
                        exon_count = int(coding_map_transcript_exon_count[transcript])
                        prot_seq_wt = coding_map_transcript_prot_seq[transcript]
                        exon_region_seq_dict = coding_map_transcript_exon_region_seq[transcript]

                        AGAIN_prot_new_acc_result = AGAIN_prot_new_acc(ivs, strand, start_codon,
                                                    start_exon, exon_count, exon_region_seq_dict,
                                                    prot_seq_wt, AGAIN_acc_seq)
                        AGAIN_prot_new_acc_seq = AGAIN_prot_new_acc_result[0]
                        AGAIN_prot_new_acc_hgvs = AGAIN_prot_new_acc_result[1]

                        AGAIN_prot_exon_skip_result = AGAIN_prot_exon_skip(ivs, strand, start_codon,
                                                      start_exon, exon_count, exon_region_seq_dict, prot_seq_wt)
                        AGAIN_prot_exon_skip_seq = AGAIN_prot_exon_skip_result[0]
                        AGAIN_prot_exon_skip_hgvs = AGAIN_prot_exon_skip_result[1]

                        # AGAIN output
                        AGAIN_var_set.add(var_info)
                        output = sample+'\t'+chrom+'\t'+str(var_pos) +'\t'+var_id+'\t'+ref+'\t'+alt+'\t'+strand+'\t'+\
                                 var_type+'\t'+gene+'\t'+transcript_ivs+'\t'+canonical+'\t'+\
                                 AGAIN_zone+'\t'+AGAIN_YAG+'\t'+str(AGAIN_bp_dist)+'\t'+str(AGAIN_acc_dist)+'\t'+\
                                 AGAIN_highrisk+'\t'+str(AGAIN_score) +'\t'+prot_seq_wt+'\t'+\
                                 AGAIN_prot_new_acc_seq+'\t'+AGAIN_prot_new_acc_hgvs+'\t'+\
                                 AGAIN_prot_exon_skip_seq+'\t'+AGAIN_prot_exon_skip_hgvs+'\n'
                        file_out.write(output)
            except:
                pass

        end_time = time.time()
        time_cost = int(end_time - start_time)
        file_mapping.close()
        file_out.close()
        os.remove(filename_bed)
        os.remove(filename_mapping)

        # log
        print('# Input variants:', str(input_var_count))
        print('# AGAIN variants:', str(len(AGAIN_var_set)))
        print('# AGAIN variants (high-risk):', str(AGAIN_var_highrisk_count))
        print('Time cost:', str(time_cost), 'seconds\n')

    except:
        print('Error occured.')
