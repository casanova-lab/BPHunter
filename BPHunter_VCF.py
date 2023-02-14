# python3.8
__author__ = "Peng Zhang"
__copyright__ = "Copyright 2023, " \
                "Laboratory of Human Genetics of Infectious Diseases, " \
                "The Rockefeller University"
__license__ = "CC BY-NC-ND 4.0"
__version__ = "ver-2, 02-10-2023"

import os
import time
import argparse


print('**********************************************************')
print('  ####   #####  #   #  #   #  #   #  #####  #####  ####   ')
print('  #   #  #   #  #   #  #   #  ##  #    #    #      #   #  ')
print('  ####   #####  #####  #   #  # # #    #    ####   ###    ')
print('  #   #  #      #   #  #   #  #  ##    #    #      #  #   ')
print('  ####   #      #   #   ###   #   #    #    #####  #   #  ')
print(' Detecting human intronic variants disrupting branchpoint ')
print('**********************************************************\n')


###
# input parameters
###
parser = argparse.ArgumentParser(description="BPHunter")
parser.add_argument("-i", "--input",
                    help="input vcf file of variants with 5 mandatory fields: CHROM, POS, ID, REF, ALT")
parser.add_argument("-g", "--genome_ref", type=str, default='GRCh37', choices=['GRCh37', 'GRCh38'],
                    help="Human reference genome assembly")
parser.add_argument("-t", "--transcript_ref", type=str, default='all', choices=['all', 'canonical'],
                    help="all/canonical transcripts?")

args = parser.parse_args()
filename_vcf = args.input
genome_ref = args.genome_ref
transcript_ref = args.transcript_ref

if transcript_ref == 'all':
    filename_BPHunter_detection = 'Data_BPHunter_'+genome_ref+'_detection_all.bed'
else:
    filename_BPHunter_detection = 'Data_BPHunter_'+genome_ref+'_detection_canonical.bed'

timestamp = str(int(time.time()))
filename = filename_vcf[0:filename_vcf.rfind('.')]
filename_bed = filename+'.bed'
filename_mapping = filename+'.mapping'
filename_output = filename+'_BPHunter_output_'+timestamp+'.txt'


###
# BPHunter MAIN
###

try:
    ###
    # read variants, convert variants from VCF to BED
    ###
    file_var_vcf = open(filename_vcf, 'r')
    file_var_bed = open(filename_bed, 'w')
    input_var_count = 0

    for eachline in file_var_vcf:
        if not eachline.startswith('#'):
            input_var_count += 1
            item = eachline.strip().split('\t')
            chrom = item[0] if 'chr' in item[0] else 'chr' + item[0]
            pos = item[1].replace(',', '')
            var_id = item[2]
            ref = item[3]
            alt = item[4]
            var_name = chrom+'*'+pos+'*'+var_id+'*'+ref+'*'+alt
            var_start = var_end = var_type = '.'

            if len(ref) == 1 and len(alt) == 1:
                var_type = 'snv'
                var_start = str(int(pos) - 1)
                var_end = pos

            elif len(ref) == 1 and len(alt) > 1:
                var_type = 'ins-' + str(len(alt) - 1) + 'nt'
                var_start = pos
                var_end = str(int(pos) + len(alt) - 1)

            elif len(ref) > 1 and len(alt) == 1:
                var_type = 'del-' + str(len(ref) - 1) + 'nt'
                var_start = pos
                var_end = str(int(pos) + len(ref) - 1)

            if (var_start != '.') and (var_end != '.') and (var_type != '.'):
                file_var_bed.write(chrom+'\t'+var_start+'\t'+var_end+'\t'+var_name+'\t.\t+\t'+var_type+'\n')
                file_var_bed.write(chrom+'\t'+var_start+'\t'+var_end+'\t'+var_name+'\t.\t-\t'+var_type+'\n')

    file_var_vcf.close()
    file_var_bed.close()

    ###
    # mapping variants to BP position
    # identify BP variant with BPHunter annotations
    ###
    os.system("bedtools intersect -wo -s"
              " -a " + filename_bed +
              " -b " + filename_BPHunter_detection +
              " > " + filename_mapping)

    file_mapping = open(filename_mapping, 'r')
    file_out = open(filename_output, 'w')
    file_out.write('CHROM\tPOS\tID\tREF\tALT\tSTRAND\tVAR_TYPE\tGENE\tTRANSCRIPT_IVS\tCANONICAL\t'
                   'BP_NAME\tBP_ACC\tBP_RANK\tBP_TOTAL\tBP_HIT\tBP_SOURCE\tCONSENSUS\t'
                   'BP_GERP\tBP_PHYL\t\tBP2_GERP\tBP2_PHYL\tBPHUNTER_HIGHRISK\tBPHUNTER_SCORE\n')

    BPHunter_var_set = set()
    BPHunter_var_highscore_count = 0
    BPHunter_var_highrisk_count = 0
    for eachline in file_mapping:
        try:
            item = eachline.strip().split('\t')
            chrom = item[0]
            var_start = int(item[1])
            var_end = int(item[2])
            var_info = item[3]
            strand = item[5]
            var_type = item[6]
            bp_motif_start = int(item[8])
            bp_motif_end = int(item[9])
            bp_name = item[10]
            gene = item[13]
            transcript_ivs = item[14]
            canonical = item[15]
            bp_acc = item[16]
            bp_rank = item[17]
            bp_total = item[18]
            bp_source = item[19]
            consensus = item[20]
            bp_maf = item[21]
            bp_gerp = item[22]
            bp_phyl = item[23]
            bp2_maf = item[24]
            bp2_gerp = item[25]
            bp2_phyl = item[26]

            var_info_list = var_info.split('*')
            var_pos = var_info_list[1]
            var_id = var_info_list[2]
            var_ref = var_info_list[3]
            var_alt = var_info_list[4]

            hit_pos_list = list()
            if strand == '+':
                bp_hit_checking = -2
                for i in range(bp_motif_start+1, bp_motif_end+1):
                    if var_start+1 <= i <= var_end:
                        hit_pos_list.append(bp_hit_checking)
                    bp_hit_checking += 1
            else:
                bp_hit_checking = 0
                for i in range(bp_motif_start+1, bp_motif_end+1):
                    if var_start+1 <= i <= var_end:
                        hit_pos_list.append(bp_hit_checking)
                    bp_hit_checking -= 1
            hit_pos_list.sort()
            bp_hit = '|'.join([str(i) for i in hit_pos_list])

            BPHunter_score = 0
            if bp_rank == '1':
                BPHunter_score += 1
            if bp_total == '1':
                BPHunter_score += 1
            if '1:' in consensus:
                BPHunter_score += 1
            if int(bp_source) > 1:
                BPHunter_score += 1

            if 0 in hit_pos_list:
                if bp_maf == '.':
                    BPHunter_score += 1
                else:
                    if float(bp_maf) < -2:
                        BPHunter_score += 1
                if bp_gerp != '.':
                    if float(bp_gerp) > 0:
                        BPHunter_score += 1
                if bp_phyl != '.':
                    if float(bp_phyl) > 0:
                        BPHunter_score += 1

            if -2 in hit_pos_list:
                if bp2_maf == '.':
                    BPHunter_score += 1
                else:
                    if float(bp2_maf) < -2:
                        BPHunter_score += 1
                if bp2_gerp != '.':
                    if float(bp2_gerp) > 0:
                        BPHunter_score += 1
                if bp2_phyl != '.':
                    if float(bp2_phyl) > 0:
                        BPHunter_score += 1

            if var_info not in BPHunter_var_set:
                BPHunter_var_set.add(var_info)

                BPHunter_highrisk = 'NO'
                if (bp_rank in ['1', '2']) and (bp_total in ['1', '2']) and ('1:' in consensus) and int(bp_source) > 1:
                    BPHunter_highrisk = 'YES'
                    BPHunter_var_highrisk_count += 1

                if BPHunter_score >= 3:
                    BPHunter_var_highscore_count += 1

                output = chrom+'\t'+var_pos+'\t'+var_id+'\t'+var_ref+'\t'+var_alt+'\t'+strand+'\t'+var_type+'\t'+\
                         gene+'\t'+transcript_ivs+'\t'+canonical+'\t'+bp_name+'\t'+bp_acc+'\t'+bp_rank+'\t'+bp_total+'\t'+\
                         bp_hit+'\t'+bp_source+'\t'+consensus+'\t'+bp_gerp+'\t'+bp_phyl+'\t'+bp2_gerp+'\t'+bp2_phyl+'\t'+\
                         BPHunter_highrisk+'\t'+str(BPHunter_score)
                file_out.write(output + '\n')

        except:
            pass

    file_mapping.close()
    file_out.close()
    os.remove(filename_bed)
    os.remove(filename_mapping)

    ###
    # log
    ###
    print('Input file:', filename_vcf)
    print('Output file:', filename_output)
    print('Genome:', genome_ref)
    print('Transcripts:', transcript_ref, '\n')
    print('# Input variants:', str(input_var_count))
    print('# BP variants:', str(len(BPHunter_var_set)))
    print('# BP variants (score >= 3):', str(BPHunter_var_highscore_count))
    print('# BP variants (high-risk):', str(BPHunter_var_highrisk_count), '\n')

except:
    print('Error occured. Please check your input file, or contact the developer: pzhang@rockefeller.edu.\n')
