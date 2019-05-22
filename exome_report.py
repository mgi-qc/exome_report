#script for generating report and SSheet for exome
#1.0 == 0.4 W/O TEMP LOCAL READ-INS



# TODO PRINT TEXT FILE FOR DATA TRANSFER DIRECTORY
__author__ = 'Thomas Antonacci'


import csv
import os
import sys
import glob
import datetime
import argparse
import subprocess
from string import Template

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


    #Ini. dicts
    results = {}
    template_file_dict = {}

    #Interested metrics;
    metrics_tracked = ['PASS_SAMPLES','ALN_FAIL']
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
            with open(report_outfile, 'w', encoding='utf-8') as fh:
                fh.write(template_file.substitute(WOID=template_file_dict['WOID'],
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
                                                  BUILDS=','.join(last_succeeded_build_id),
                                                  RESULTS_SPREADSHEET=SSheet_outfile))
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
