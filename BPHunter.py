# python3.8
import os
import argparse

###
# input parameters
###
print('**********************************************************')
print('  ####   #####  #   #  #   #  #   #  #####  #####  ####   ')
print('  #   #  #   #  #   #  #   #  ##  #    #    #      #   #  ')
print('  ####   #####  #####  #   #  # # #    #    ####   ###    ')
print('  #   #  #      #   #  #   #  #  ##    #    #      #  #   ')
print('  ####   #      #   #   ###   #   #    #    #####  #   #  ')
print(' Detecting human variants disrupting intronic branchpoint ')
print('**********************************************************\n')

parser = argparse.ArgumentParser(description="BPHunter")
parser.add_argument("-i", "--input", help="input vcf file of variants with 5 mandatory fields: CHROM, POS, ID, REF, ALT")
parser.add_argument("-g", "--genome", type=str, default='hg19', choices=['hg19', 'hg38'], help="Human reference genome assembly")
parser.add_argument("-c", "--canon", type=str, default='no', choices=['no', 'yes'], help="canonical transcripts?")

args = parser.parse_args()
filename_vcf = args.input
genome = args.genome
canonical = args.canon

filename = filename_vcf[0:filename_vcf.rfind('.')]
filename_bed = filename+'.bed'
filename_mapping = filename+'.bphunter.mapping'
filename_out = filename+'.bphunter.txt'

if canonical == 'no':
    filename_bphunter_ref = 'Data_BPHunter_'+genome+'_Detection_All.bed'
elif canonical == 'yes':
    filename_bphunter_ref = 'Data_BPHunter_'+genome+'_Detection_Canonical.bed'

try:
    ###
    # convert input variants from VCF to BED with addtional var_type
    ###
    file_var_vcf = open(filename_vcf, 'r')
    file_var_bed = open(filename_bed, 'w')
    count_var_input = 0
    for eachline in file_var_vcf:
        if not eachline.startswith('#'):
            try:
                item = eachline.strip().split('\t')
                chrom = item[0]
                pos = item[1]
                ref = item[3]
                alt = item[4]
                if 'chr' not in chrom:
                    chrom = 'chr' + chrom
                var_name = chrom + '*' + pos + '*' + ref + '*' + alt
                count_var_input += 1

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

            except:
                pass

    file_var_vcf.close()
    file_var_bed.close()

    ###
    # overlapping the input variants with BPHunter reference data
    ###
    os.system("bedtools intersect -wo -s -a " + filename_bed + " -b " + filename_bphunter_ref + " > " + filename_mapping)

    ###
    # processing the variants overlapped with BPHunter
    ###
    file_mapping = open(filename_mapping, 'r')
    file_BPHunter_out = open(filename_out, 'w')
    file_BPHunter_out.write('CHROM\tPOS\tREF\tALT\tSTRAND\tVAR_TYPE\tGENE\tBP_NAME\tBP_TYPE\tRANK\tHIT_POS\tDIST_3SS\t')
    file_BPHunter_out.write('MAF\tGERP\tPHYLOP\tENERGY\tCONSENSUS\tEV_SCORE\tSOURCE\tINTRON_TYPE\tTRANSCRIPT_INTRON#\n')
    BPHunter_var_set = set()
    BPHunter_output_list = list()
    for eachline in file_mapping:
        item = eachline.strip().split('\t')
        chrom = item[0]
        var_start = int(item[1])
        var_end = int(item[2])
        var_name = item[3]
        strand = item[5]
        var_type = item[6]

        bp_motif_start = int(item[8])
        bp_motif_end = int(item[9])
        bp_name = item[10]
        score = item[11]
        source = item[13]
        gene = item[14]
        transcript_intron = item[15]
        intron_type = item[16]
        dist_3ss = item[17]
        rank_total = item[18]
        energy = item[19]
        maf = item[20]
        gerp = item[21]
        phylop = item[22]
        consensus = item[23]
        overlap = item[24]

        bp_type = bp_name.split('_')[0]
        var_info_list = var_name.split('*')
        var_pos = var_info_list[1]
        var_ref = var_info_list[2]
        var_alt = var_info_list[3]

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
        hit_pos = '|'.join([str(i) for i in hit_pos_list])

        BPHunter_output = chrom+'\t'+var_pos+'\t'+var_ref+'\t'+var_alt+'\t'+strand+'\t'+var_type+'\t'+gene+'\t'
        BPHunter_output += bp_name+'\t'+bp_type+'\t'+rank_total+'\t'+hit_pos+'\t'+dist_3ss+'\t'
        BPHunter_output += maf+'\t'+gerp+'\t'+phylop+'\t'+energy+'\t'+consensus+'\t'+score+'\t'+source+'\t'
        BPHunter_output += intron_type+'\t'+transcript_intron
        BPHunter_output_list.append(BPHunter_output)
        BPHunter_var_set.add(var_name)

    os.remove(filename_bed)
    os.remove(filename_mapping)

    ###
    # BPHunter output
    ###
    count_var_output = len(BPHunter_var_set)
    count_var_bp = len(BPHunter_output_list)
    BPHunter_output_list.sort()
    file_BPHunter_out.write('\n'.join(BPHunter_output_list) + '\n')

    print('Input file:', filename_vcf)
    print('Output file:', filename_out)

    if genome == 'hg19':
        print('Genome assembly: hg19/GRCh37')
    elif genome == 'hg38':
        print('Genome assembly: hg38/GRCh38')

    if canonical == 'no':
        print('Transcripts: all')
    elif canonical == 'yes':
        print('Transcripts: canonical')

    print('')
    print('# Input variants:', str(count_var_input))
    print('# BPHunter-detected variants:', str(count_var_output))
    print('# Variant-BP associations:', str(count_var_bp))
    print('')

except:
    print('Error occured. Please check your input file, or contact the developer.')
