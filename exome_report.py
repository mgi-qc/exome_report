#script for generating report and SSheet for exome
#1.0 == 0.4 W/O TEMP LOCAL READ-INS

__author__ = 'Thomas Antonacci'


import csv
import os
import sys
import glob
import datetime
import argparse
import subprocess
from string import Template


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

mm_dd_yy = datetime.datetime.now().strftime("%m%d%y")

#Check for metrics file;
if not glob.glob('*.cwl.metrics.*.tsv'):
    sys.exit('cwl.metrics file not found')
else:
    metrics_files = glob.glob('*.cwl.metrics.*.tsv')

#Check, open, and create template file using Template;
if not os.path.isfile('/gscmnt/gc2783/qc/GMSworkorders/reports/exome_report_template.txt'):
    sys.exit('Template file not found.')

with open('/gscmnt/gc2783/qc/GMSworkorders/reports/exome_report_template.txt', 'r', encoding='utf-8') as fh:
    template = fh.read()
    template_file = Template(template)


filename_list = []

for file in metrics_files:

    PCT_20X_check = False

    #Ini. outgoing files;
    file_name = file.split('.')[0]
    SSheet_outfile = '{}.cwl.results.{}.tsv'.format(file_name, mm_dd_yy)
    report_outfile = '{}.cwl.report.{}.txt'.format(file_name, mm_dd_yy)

    # Ini. dicts
    results = {}
    template_file_dict = {}

    # Interested metrics;
    metrics_tracked = ['PASS_SAMPLES', 'ALN_FAIL']

    #Additional Metrics prompts and input
    ad_met_check = False
    mt_check = False
    #Haploid prompt

    while not ad_met_check:
        ad_met_in = input("Would you like to require additional metrics for {}? y/n: ".format(file_name))

        if ad_met_in is 'y':
            ad_met_check = True

            #Mean Target prompt
            while not mt_check:
                mt_in = input('\nPlease enter a MEAN_TARGET_COVERAGE value for {}, or enter 0 to skip metric: '
                              .format(file_name))

                if not is_number(mt_in) or float(mt_in) < 0:  # not in constratints
                    print('\nPlease enter a value >= 0 ')
                elif float(mt_in) > 0:
                    mt_value = float(mt_in)
                    mt_check = True
                    print('MEAN_TARGET_COVERAGE minimum set to {}'.format(mt_value))
                    add_met = 'MEAN_TARGET_COVERAGE (minimum requirement) : {}'.format(mt_value)
                    metrics_tracked.extend(['MEAN_TAR_COV_PASS', 'MEAN_TAR_COV_FAIL'])
                else:
                    print('Skipping MEAN_TARGET_COVERAGE')
                    mt_check = True
        elif ad_met_in is 'n':
            print('Skipping addtional metrics')
            add_met = 'No other metric required/reviewed for assignment of QC pass/fail judgement'
            ad_met_check = True
        else:
            print('Please enter y or n')


    for metric in metrics_tracked:
        template_file_dict[metric] = 0

    #Metrics File Open, Check Metrics, Generate 'results', Get Totals;
    with open(file, 'r') as fh,open (SSheet_outfile, 'w') as of:
        metrics_dict = csv.DictReader(fh, delimiter='\t')
        header = metrics_dict.fieldnames


        ofd = csv.DictWriter(of, fieldnames=header, delimiter='\t')
        header.extend(['QC_Status','QC_failed_metrics'])
        ofd.writeheader()

        #Ini. build id and totals for averages in template
        last_succeeded_build_id = []
        tot_pct_tar_bases = tot_pf_aln_bases = tot_pct_usbl_tar = tot_pct_usbl_bait = 0
        tot_mean_tar_cov = tot_pct_exc_off = tot_pct_exc_dup =tot_per_dup = 0

        count = 0
        per_dup_count = 0

        for line in metrics_dict:

            template_file_dict['WOID'] = line['WorkOrder']


            #Check metrics based on Metrics File;
            if 'PCT_TARGET_BASES_20X' in line:
                if float(line['PCT_TARGET_BASES_20X']) < 0.70:
                    line['QC_Status'] = 'Fail'
                    line['QC_failed_metrics'] = 'PCT_TARGET_BASES_20X'
                    template_file_dict['ALN_FAIL'] += 1
                else:
                    line['QC_Status'] = 'Pass'
                    line['QC_failed_metrics'] = 'NA'
                    template_file_dict['PASS_SAMPLES'] += 1

                if 'MEAN_TAR_COV_PASS' in metrics_tracked:
                    if 'MEAN_TARGET_COVERAGE' in line:
                        if float(line['MEAN_TARGET_COVERAGE']) < mt_value:
                            line['QC_failed_metrics'] += 'MEAN_TARGET_COVERAGE'
                            template_file_dict['MEAN_TAR_COV_FAIL'] += 1
                        else:
                            template_file_dict['MEAN_TAR_COV_PASS'] += 1


                tot_pct_tar_bases += float(line['PCT_TARGET_BASES_20X'])
                tot_pct_usbl_tar += float(line['PCT_USABLE_BASES_ON_TARGET'])
                tot_pct_usbl_bait += float(line['PCT_USABLE_BASES_ON_BAIT'])
                tot_mean_tar_cov += float(line['MEAN_TARGET_COVERAGE'])
                tot_pf_aln_bases += float(line['PF_BASES_ALIGNED'])
                tot_pct_exc_off += float(line['PCT_EXC_OFF_TARGET'])
                tot_pct_exc_dup += float(line['PCT_EXC_DUPE'])

                if line['PERCENT_DUPLICATION'] == 'FNF':
                    avg_per_dup = 'NA'
                else:
                    tot_per_dup += float(line['PERCENT_DUPLICATION'])
                    per_dup_count += 1



                last_succeeded_build_id.append(line['last_succeeded_build'])
                #fill line in results;
                ofd.writerow(line)
                count += 1



            else:
                PCT_20X_check = True

        if tot_per_dup > 0:
            avg_per_dup = tot_per_dup / per_dup_count

        if not PCT_20X_check:

            #Prompt for seq notes:
            seq_check = False
            while not seq_check:
                seq_in = input('Would you like to add a SEQUENCING_NOTE? y/n: ')

                if seq_in is 'y':
                    seq_check = True
                    print('See https://confluence.ris.wustl.edu/pages/viewpage.action?spaceKey=AD&title=WorkOrder+{} for notes;'
                          ' enter into console, and enter double return when finished: '.format(file_name))
                    seq_notes = []
                    while True:
                        note_line = input()
                        if note_line:
                            seq_notes.append(note_line)
                        else:
                            break

                elif seq_in is 'n':
                    seq_check = True
                    seq_notes = ['']
                    print('Skipping SEQUENCING_NOTE')
                else:
                    print('Please enter y or n')
            test = '\n'.join(seq_notes)
            #write report
            with open(report_outfile, 'w', encoding='utf-8') as fhr:
                fhr.write(template_file.substitute(WOID=template_file_dict['WOID'],
                                                   ADDITIONAL_METRICS = add_met,
                                                   SAMPLE_NUMBER=count,
                                                   PASS_SAMPLES=template_file_dict['PASS_SAMPLES'],
                                                   ALN_FAIL=template_file_dict['ALN_FAIL'],
                                                   PCT_TARGET_BASES_20X=tot_pct_tar_bases / count,
                                                   PCT_USABLE_BASES_ON_TARGET=tot_pct_usbl_tar / count,
                                                   PCT_USABLE_BASES_ON_BAIT=tot_pct_usbl_bait / count,
                                                   MEAN_TARGET_COVERAGE=tot_mean_tar_cov / count,
                                                   PF_BASES_ALIGNED=tot_pf_aln_bases / count,
                                                   PCT_EXC_OFF_TARGET=tot_pct_exc_off / count,
                                                   PCT_EXC_DUPE=tot_pct_exc_dup / count,
                                                   PERCENT_DUPLICATION=avg_per_dup,
                                                   SEQUENCING_NOTE = test,
                                                   MEAN_TAR_COV_PASS = template_file_dict['MEAN_TAR_COV_PASS'],
                                                   MEAN_TAR_COV_FAIL = template_file_dict['MEAN_TAR_COV_FAIL'],
                                                   RESULTS_SPREADSHEET = SSheet_outfile))
            filename_list.append(file_name)

            builds = ','.join(last_succeeded_build_id)

            with open('{}.Data_transfer_help.txt'.format(template_file_dict['WOID']), 'w') as df:
                df.write('Data Transfer Directory =\ncd to parent data dir\ncd to model_data'
                      '\nmkdir data_transfer/{}\ngenome model cwl-pipeline prep-for-transfer --md5sum'
                      ' --directory=full_path../data_transfer/{}  --builds {}'
                      ' or model_groups.project.id {}'.format(template_file_dict['WOID'], template_file_dict['WOID'],
                                                              builds, template_file_dict['WOID']))
        else:
            print('\nNo report generated for {}; PCT_TARGET_BASES_20X not found.'.format(file))

for file in filename_list:
    print('Report for {} generated.'.format(file))
