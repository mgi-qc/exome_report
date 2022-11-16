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

parser = argparse.ArgumentParser()
parser.add_argument('-nod', help='Turn off directory creation', action='store_true')
args = parser.parse_args()


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

DT_TRANSFER_DIR = '/storage1/fs1/gtac-mgi/Active/Bioinformatics_analysis/mds_contract_gms_compute1/mds_dt'
dt_dir_status = 'NA'

#original function definition
# def data_dir_check(dir_list, woid, date):
def data_dir_check(woid, date):
    """Create and return transfer directory if 'model' found in dir path."""

    # return NA if no transfer dir found
    # transfer_dir = 'NA'

    #variable for new DT directory
    transfer_dir = os.path.join(DT_TRANSFER_DIR, '{}_{}/'.format(woid, date))

    #check if it exists
    #if exists, add note to helper file
    if os.path.isdir(transfer_dir):
        dt_dir_status = 'DT directory exists'
        return transfer_dir
    
    #if doesn't exist, mkdir
    os.mkdir(transfer_dir) 
         
    #orginal function removed for MDS
    # iterate over data dirs
    # for directory in dir_list:
    #     # if model found, create transfer dir and return path
    #     if os.path.isdir(directory) and 'model' in directory:

    #         dir_path_items = directory.split('/')

    #         for no, d in enumerate(dir_path_items):

    #             if 'model' in d:

    #                 model_directory = '/'.join(dir_path_items[:no + 1]) + '/'
    #                 transfer_dir = os.path.join(model_directory, 'data_transfer/{}_{}/'.format(woid, date))

    #                 if os.path.isdir(transfer_dir):
    #                     print('Transfer Directory already exists: {}'.format(transfer_dir))
    #                     return 'NA'

    #                 if os.path.isdir(model_directory) and not os.path.isdir(transfer_dir):
    #                     try:
    #                         os.mkdir(transfer_dir)
    #                     except OSError:
    #                         # raise OSError("Can't create destination directory {}!".format(transfer_dir))
    #                         return 'NA'
    #                     print('Data transfer directory created:\n{}'.format(transfer_dir))
    #                     return transfer_dir

    return transfer_dir


mm_dd_yy = datetime.datetime.now().strftime("%m%d%y")

# Check for metrics file;
if not glob.glob('*.cwl.metrics.*.tsv'):
    sys.exit('cwl.metrics file not found')
else:
    metrics_files = glob.glob('*.cwl.metrics.*.tsv')

template_text_full_path = '/storage1/fs1/gtac-mgi/Active/Bioinformatics_analysis/mds_contract_gms_compute1/' \
                          'mds_analysis/bin/qc/cwl_exome_report/exome_report_template.txt'

# Check, open, and create template file using Template;
if not os.path.isfile(template_text_full_path):
    sys.exit('Template file not found.')

with open(template_text_full_path, 'r', encoding='utf-8') as fh:
    template = fh.read()
    template_file = Template(template)

filename_list = []

for file in metrics_files:

    PCT_20X_check = False
    data_directories = []

    # Ini. outgoing files;
    file_name = file.split('.')[0]
    file_date = file.split('.')[-2]
    SSheet_outfile = '{}.cwl.results.{}.tsv'.format(file_name, file_date)
    report_outfile = '{}.cwl.report.{}.txt'.format(file_name, file_date)

    # Ini. dicts
    results = {}
    template_file_dict = {}

    # Interested metrics;
    metrics_tracked = ['PASS_SAMPLES', 'ALN_FAIL', 'MEAN_TAR_COV_PASS', 'MEAN_TAR_COV_FAIL']

    mt_value = 600

    add_met = 'MEAN_TARGET_COVERAGE (minimum requirement) : {}'.format(mt_value)

    # Additional Metrics prompts and input
    # ad_met_check = False
    # mt_check = False
    # Haploid prompt

    # while not ad_met_check:
    #     ad_met_in = input("\nConfluence link:\nhttps://confluence.ris.wustl.edu/pages/viewpage.action?spaceKey=AD&title=WorkOrder+{}"
    #                       " \nWould you like to require additional metrics for {}? y/n: ".format(file_name, file_name))

    #     if ad_met_in is 'y':
    #         ad_met_check = True

    #         # Mean Target prompt
    #         while not mt_check:
    #             mt_in = 600
    #             #input('\nPlease enter a MEAN_TARGET_COVERAGE value for {}, or enter 0 to skip metric: '.format(file_name))

    #             if not is_number(mt_in) or float(mt_in) < 0:  # not in constratints
    #                 print('\nPlease enter a value >= 0 ')
    #             elif float(mt_in) > 0:
    #                 mt_value = float(mt_in)
    #                 mt_check = True
    #                 metrics_tracked.extend(['MEAN_TAR_COV_PASS', 'MEAN_TAR_COV_FAIL'])
    #                 print('MEAN_TARGET_COVERAGE minimum set to {}'.format(mt_value))
    #                 add_met = 'MEAN_TARGET_COVERAGE (minimum requirement) : {}'.format(mt_value)
    #             else:
    #                 print('Skipping MEAN_TARGET_COVERAGE')
    #                 add_met = 'No other metric required/reviewed for assignment of QC pass/fail judgement'
    #                 mt_check = True
    #     elif ad_met_in is 'n':
    #         print('Skipping addtional metrics')
    #         add_met = 'No other metric required/reviewed for assignment of QC pass/fail judgement'
    #         ad_met_check = True
    #     else:
    #         print('Please enter y or n')

    for metric in metrics_tracked:
        template_file_dict[metric] = 0

    # Metrics File Open, Check Metrics, Generate 'results', Get Totals;
    with open(file, 'r') as fh,open (SSheet_outfile, 'w') as of:
        metrics_dict = csv.DictReader(fh, delimiter='\t')
        header = metrics_dict.fieldnames

        ofd = csv.DictWriter(of, fieldnames=header, delimiter='\t')
        header.extend(['QC_Status','QC_failed_metrics'])
        ofd.writeheader()

        # Ini. build id and totals for averages in template
        last_succeeded_build_id = []
        tot_pct_tar_bases = tot_pf_aln_bases = tot_pct_usbl_tar = tot_pct_usbl_bait = 0
        tot_mean_tar_cov = tot_pct_exc_off = tot_pct_exc_dup = tot_per_dup = 0

        count = 0
        per_dup_count = 0

        for line in metrics_dict:
            failed_metrics = []

            template_file_dict['WOID'] = line['WorkOrder']
            data_directories.append(line['data_directory'])



            # Check metrics based on Metrics File;
            if 'PCT_TARGET_BASES_20X' in line and is_number(line['PCT_TARGET_BASES_20X']):
                if float(line['PCT_TARGET_BASES_20X']) < 0.70:
                    line['QC_Status'] = 'Fail'
                    failed_metrics.append('PCT_TARGET_BASES_20X')
                    template_file_dict['ALN_FAIL'] += 1
                else:
                    line['QC_Status'] = 'Pass'
                    line['QC_failed_metrics'] = 'NA'
                    template_file_dict['PASS_SAMPLES'] += 1

                if 'MEAN_TAR_COV_PASS' in metrics_tracked:
                    if 'MEAN_TARGET_COVERAGE' in line:
                        if float(line['MEAN_TARGET_COVERAGE']) < mt_value:
                            failed_metrics.append('MEAN_TARGET_COVERAGE')
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

                # if there are failed metrics, overwrite NA with failed metrics.
                if len(failed_metrics) > 0:
                    line['QC_failed_metrics'] = ','.join(failed_metrics)
                else:
                    line['QC_failed_metrics'] = 'NA'

                # fill line in results;
                ofd.writerow(line)
                count += 1

            else:
                PCT_20X_check = True

        if tot_per_dup > 0:
            avg_per_dup = tot_per_dup / per_dup_count

        if not PCT_20X_check:

            #sequencing notes no longer needed for MDS
            #remove sequencing notes from template

            # Prompt for seq notes:
            # seq_check = False
            # while not seq_check:
            #     seq_in = input('Would you like to add a SEQUENCING_NOTE? y/n: ')

            #     if seq_in is 'y':
            #         print('Enter "return q return" to exit SEQUENCING_NOTE prompt.')
            #         seq_check = True
            #         seq_notes = []
            #         while True:
            #             note_line = input()
            #             if note_line != 'q':
            #                 seq_notes.append(note_line)
            #             else:
            #                 break

            #     elif seq_in is 'n':
            #         seq_check = True
            #         seq_notes = ['']
            #         print('Skipping SEQUENCING_NOTE')
            #     else:
            #         print('Please enter y or n')

            if 'MEAN_TAR_COV_PASS' in template_file_dict:
                MEAN_TAR_PASS = template_file_dict['MEAN_TAR_COV_PASS']
                MEAN_TAR_FAIL = template_file_dict['MEAN_TAR_COV_FAIL']
            else:
                MEAN_TAR_PASS = 'NA'
                MEAN_TAR_FAIL = 'NA'

            transfer_data_directory = 'NA'
            if not args.nod:
                transfer_data_directory = data_dir_check(template_file_dict['WOID'], mm_dd_yy)

            # write report
            with open(report_outfile, 'w', encoding='utf-8') as fhr:
                fhr.write(template_file.substitute(WOID=template_file_dict['WOID'],
                                                   ADDITIONAL_METRICS=add_met,
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
                                                   MEAN_TAR_COV_PASS=MEAN_TAR_PASS,
                                                   MEAN_TAR_COV_FAIL=MEAN_TAR_FAIL,
                                                #    SEQUENCING_NOTE='\n'.join(seq_notes),
                                                   TRANSFER_DIR=transfer_data_directory,
                                                   RESULTS_SPREADSHEET=SSheet_outfile,
                                                   REPORT_FILE=report_outfile))

            filename_list.append(file_name)

            builds = ','.join(last_succeeded_build_id)

            with open('{}.Data_transfer_help.{}.txt'.format(template_file_dict['WOID'], file_date), 'w') as df:
                df.write('Data Transfer Directory Status = {s}\n'
                      'Data Transfer Directory ={td}\ncd to parent data dir\ncd to model_data'
                      '\nmkdir data_transfer/{w}\nTransfer Commands:\n\ngenome model cwl-pipeline prep-for-transfer --md5sum'
                      ' --directory={td}  --builds {b}\n\n'
                      'genome model cwl-pipeline prep-for-transfer --md5sum'
                      ' --directory={td} model_groups.project.id={w}\n'.format(td=transfer_data_directory, w=template_file_dict['WOID'], s = dt_dir_status, b=builds,))

            if dt_dir_status == 'NA':
                print('Running: genome model cwl-pipeline prep-for-transfer --md5sum')
                subprocess.run(['genome','model', 'cwl-pipeline', 'prep-for-transfer', '--md5sum', f'--directory={transfer_data_directory}', '--builds', builds])
            print('------------------------')
        
        else:
            print('\nNo report generated for {}; PCT_TARGET_BASES_20X not found.'.format(file))
            print('------------------------')
for file in filename_list:
    print('Report for {} generated.'.format(file))
