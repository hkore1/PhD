#!/usr/bin/env python

"""
Filter Chimeras with Genome Reference (Captus Only)

This script processes Captus assembly data to identify and filter out chimeric sequences
by mapping them to a reference genome. It analyzes paralogous gene sequences, detects
chimeras based on mapping patterns, and generates reports and visualizations of the results.

The script performs the following main tasks:
1. Creates a mapping reference from a high-quality genome
2. Parses target gene files to extract gene IDs
3. Maps Captus gene sequences to the reference genome
4. Identifies chimeric sequences based on mapping patterns
5. Filters out chimeric sequences based on configurable thresholds
6. Generates reports and visualizations of the results

########################################################################################################################
Additional information:

This script is designed to work with Captus assembly data and requires:
- A high-quality reference genome in FASTA format
- A target file in FASTA format as used for the Captus runs
- A folder containing Captus sample folders with assembled data

The script uses BBmap.sh for genome mapping and generates various reports and visualizations
to help analye the results.
########################################################################################################################
"""

import logging
import sys
import argparse
import os
import socket
import glob
from concurrent.futures.process import ProcessPoolExecutor
from concurrent.futures import wait, as_completed
from concurrent.futures import wait
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from collections import defaultdict
import datetime
import itertools
import subprocess
import traceback
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
import pandas as pd
from collections import Counter
import numpy as np

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

# f-strings will produce a 'SyntaxError: invalid syntax' error if not supported by Python version:
f'Must be using Python 3.6 or higher.'

if sys.version_info[0:2] < (3, 6):
    sys.exit(f'Must be using Python 3.6 or higher. You are using version {sys.version_info[0]}.{sys.version_info[1]}.')

# Set version:
__version__ = '0.0.2'

########################################################################################################################
########################################################################################################################
# Get current working directory and host name:

cwd = os.getcwd()
host = socket.gethostname()


# Configure logger:
def setup_logger(name,
                 log_file,
                 log_directory=None,
                 console_level=logging.INFO,
                 file_level=logging.DEBUG,
                 logger_object_level=logging.DEBUG):
    """
    Function to create a logger instance.

    By default, logs level DEBUG and above to file.
    By default, logs level INFO and above to stdout and file.

    Args:
        name (str): Name for the logger instance
        log_file (str): Filename for log file
        log_directory (str, optional): Name for the log directory to create. Defaults to None.
        console_level (int, optional): Level for logging to console. Defaults to logging.INFO.
        file_level (int, optional): Level for logging to file. Defaults to logging.DEBUG.
        logger_object_level (int, optional): Level for logger object. Defaults to logging.DEBUG.

    Returns:
        logging.Logger: A logger object
    """

    # Get date and time string for log filename:
    date_and_time = datetime.datetime.now().strftime("%Y-%m-%d-%H_%M_%S")

    # Create log directory if supplied:
    if log_directory:
        if not os.path.exists(log_directory):
            os.makedirs(f'{log_directory}')
        log_file_name = f'{log_directory}/{log_file}_{date_and_time}.log'
    else:
        log_file_name = f'{log_file}_{date_and_time}.log'

    # Log to file:
    file_handler = logging.FileHandler(f'{log_file_name}', mode='w')
    file_handler.setLevel(file_level)
    file_format = logging.Formatter('%(asctime)s - %(filename)s - %(name)s - %(funcName)s - %(levelname)s - %('
                                    'message)s')
    file_handler.setFormatter(file_format)

    # Log to Terminal (stdout):
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(console_level)
    console_format = logging.Formatter('%(message)s')
    console_handler.setFormatter(console_format)

    # Setup logger:
    logger_object = logging.getLogger(name)
    logger_object.setLevel(logger_object_level)  # Default level is 'WARNING'

    # Add handlers to the logger
    logger_object.addHandler(console_handler)
    logger_object.addHandler(file_handler)

    return logger_object


# Create logger(s):
logger = setup_logger(__name__,
                      'filter_chimeras_with_genome_ref',
                      log_directory='logs')


########################################################################################################################
########################################################################################################################
# Define general functions:


def createfolder(directory):
    """
    Attempts to create a directory named after the name provided, and provides an error message on failure.

    Args:
        directory (str): Name of directory to create

    Returns:
        str: The created directory path
    """

    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
        return directory
    except OSError:
        logger.info(f'Error: Creating directory: {directory}')


def file_exists_and_not_empty(file_name):
    """
    Check if file exists and is not empty by confirming that its size is not 0 bytes.

    Args:
        file_name (str): Path to the file to check

    Returns:
        bool: True if the file exists and is not empty, False otherwise
    """
    # Check if file exist and is not empty
    return os.path.isfile(file_name) and not os.path.getsize(file_name) == 0


def parse_target_file(target_file_fasta):
    """
    Parse a target FASTA file to extract gene IDs.

    This function reads a FASTA file and extracts gene IDs from sequence names.
    The gene ID is assumed to be the last part of the sequence name after the last hyphen.

    Args:
        target_file_fasta (str): Path to the target FASTA file

    Returns:
        list: A sorted list of unique gene IDs extracted from the FASTA file
    """

    target_file_basename = os.path.basename(target_file_fasta)

    logger.info(f'{"[INFO]:":10} Parsing target file "{target_file_basename}" to recover gene IDs...')

    gene_ids = set()

    for seq in SeqIO.parse(target_file_fasta, 'fasta'):
        gene_id = seq.name.split('-')[-1]
        gene_ids.add(gene_id)

    logger.info(f'{"[INFO]:":10} Number of unique gene IDs: {len(gene_ids)}')

    return sorted(list(gene_ids))


def create_genome_mapping_reference(reference_genome_fasta):
    """
    Create a mapping reference for a genome using BBmap.sh.

    This function creates a mapping reference for a genome using BBmap.sh.
    If the reference already exists, it skips the creation step.

    Args:
        reference_genome_fasta (str): Path to the reference genome FASTA file

    Raises:
        ValueError: If there is an issue running BBmap.sh
    """

    genome_basename = os.path.basename(reference_genome_fasta)

    logger.info(f'{"[INFO]:":10} Creating mapping reference for genome "{genome_basename}" using BBmap.sh...')

    expected_file = f'ref/genome/1/summary.txt'

    try:
        assert file_exists_and_not_empty(expected_file)
        logger.info(f'{"[INFO]:":10} BBmap.sh reference exists, skipping!')

    except AssertionError:
        try:
            command = f'bbmap.sh -Xmx30g k=13 ref={reference_genome_fasta}'

            result = subprocess.run(command,
                                    universal_newlines=True,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    check=True,
                                    shell=True)

            logger.debug(f'BBmap.sh check_returncode() is: {result.check_returncode()}')
            logger.debug(f'BBmap.sh stdout is: {result.stdout}')
            logger.debug(f'BBmap.sh stderr is: {result.stderr}')

        except subprocess.CalledProcessError as exc:
            logger.error(f'BBmap.sh FAILED. Output is: {exc}')
            logger.error(f'BBmap.sh stdout is: {exc.stdout}')
            logger.error(f'BBmap.sh stderr is: {exc.stderr}')

            raise ValueError('There was an issue running BBmap.sh. Check input files!')


def map_captus_stitch_contigs_to_genome_mp(gene_list,
                                           reference_genome_fasta,
                                           captus_sample_parent_folder,
                                           min_seq_length=75,
                                           min_length_percentage=0.80,
                                           min_contig_number_percentage=0.80,
                                           min_samples_threshold=0.75,
                                           min_num_paralogs_per_sample=0,
                                           outdir_reports=None,
                                           outdir_captus_seqs_mapped=None,
                                           outdir_non_chimeras=None,
                                           outdir_all_seqs=None,
                                           pool=1,
                                           threads=1):
    """
    Map Captus stitched contigs to a reference genome using multiprocessing.

    This function maps Captus stitched contigs to a reference genome using BBmap.sh.
    It processes multiple samples in parallel and identifies chimeric sequences.
    The function generates various reports and visualizations of the results.

    Args:
        gene_list (list): List of gene IDs to process
        reference_genome_fasta (str): Path to the reference genome FASTA file
        captus_sample_parent_folder (str): Path to the parent folder containing Captus sample folders
        min_seq_length (int, optional): Minimum sequence length to map. Defaults to 75.
        min_length_percentage (float): Minimum percentage of the paralog sequence length remaining after
                                      filtering contig hits via <min_seq_length> for chimera detection to be performed.
                                      Otherwise, no sequence is retained for this sample/gene. Defaults to 0.80.
        min_contig_number_percentage (float): Minimum percentage of the total number of contig hits remaining for a
                                             given paralog after filtering contig hits via <min_seq_length> for
                                             chimera detection to be performed. Otherwise, no sequence is retained
                                             for this sample/gene. Defaults to 0.80.
        min_samples_threshold (float, optional): For a given gene, the minimum percentage of total samples to have '
                                               '>= <min_num_paralogs_per_sample> non-chimeric sequences for sequences
                                               to be written to file. Defaults to 0.75.
        min_num_paralogs_per_sample (int, optional): For a given gene for a given sample, the minimum number of
                                                   non-chimeric paralog sequences recovered for the sequences to be
                                                   written to file. Defaults to 0.
        outdir_reports (str, optional): Output directory for report files. Defaults to None.
        outdir_captus_seqs_mapped (str, optional): Output directory for mapped Captus sequences. Defaults to None.
        outdir_non_chimeras (str, optional): Output directory for non-chimeric sequences. Defaults to None.
        outdir_all_seqs (str, optional): Output directory for all sequences. Defaults to None.
        pool (int, optional): Number of processes to run concurrently. Defaults to 1.
        threads (int, optional): Number of threads to use for each concurrent process. Defaults to 1.
    """

    # Create output directories:
    outdir_reports = createfolder(outdir_reports)
    outdir_captus_seqs_mapped = createfolder(outdir_captus_seqs_mapped)
    outdir_non_chimeras = createfolder(outdir_non_chimeras)
    outdir_all_seqs = createfolder(outdir_all_seqs)

    genome_basename = os.path.basename(reference_genome_fasta)
    logger.info(f'{"[INFO]:":10} Mapping Captus seqs to genome {genome_basename} using BBmap.sh...')

    captus_samples = sorted(list(glob.glob(f'{captus_sample_parent_folder}/*')))
    captus_samples_number = len(captus_samples)
    logger.info(f'{"[INFO]:":10} Number of Captus samples found in folder "{captus_sample_parent_folder}": '
                f'{captus_samples_number}')

    combined_samples_dict = dict()

    with (ProcessPoolExecutor(max_workers=pool) as pool):
        future_results = [pool.submit(map_captus_stitched_contigs,
                                      sample,
                                      gene_list,
                                      min_seq_length,
                                      min_length_percentage,
                                      min_contig_number_percentage,
                                      outdir_captus_seqs_mapped,
                                      threads=threads)

                          for sample in captus_samples]

        for future in as_completed(future_results):

            try:
                (sample_name,
                 sample_dict) = future.result()

                combined_samples_dict[sample_name] = sample_dict

            except Exception as error:
                print(f'Error raised: {error}')
                tb = traceback.format_exc()
                print(f'traceback is:\n{tb}')
                sys.exit()

        wait(future_results, return_when="ALL_COMPLETED")

    ####################################################################################################################
    # Write filtered fasta files:
    ####################################################################################################################
    total_number_of_samples = len(combined_samples_dict)
    min_threshold = min_samples_threshold * total_number_of_samples
    logger.info(f'{"[INFO]:":10} Total number of samples: {total_number_of_samples}')
    logger.info(f'{"[INFO]:":10} Minimum number of samples threshold: {min_threshold}')

    ####################################################################################################################
    # Get a dictionary of the scipio/captus NUC_genes.fna seqrecords for each sample, for each gene/paralog seq:
    ####################################################################################################################
    sample_to_gene_seqrecord_dict = defaultdict(lambda: defaultdict(dict))

    for sample_path in sorted(captus_samples):
        sample_basename = os.path.basename(sample_path)
        sample_name = sample_basename.split('__')[0]

        combined_gene_fasta_file = f'{sample_path}/01_coding_NUC/NUC_coding_NT.fna'  # contains paralogs too
        fna_sequences = list(SeqIO.parse(combined_gene_fasta_file, 'fasta'))

        for seq in fna_sequences:
            name_split = seq.name.split('__')
            gene = name_split[1]
            if len(name_split) == 3:  # i.e. there are paralogs for this gene, so it has a numbered suffix e.g. 03
                paralog_number = int(name_split[2])
                sample_to_gene_seqrecord_dict[sample_name][gene][f'paralog_{int(paralog_number) + 1}'] = seq
            else:
                sample_to_gene_seqrecord_dict[sample_name][gene][f'paralog_1'] = seq  # only one seq recovered

    ####################################################################################################################
    # For each non-chimeric captus paralog sequence, recover the seqrecord:
    ####################################################################################################################
    non_chimeric_seqs_dict = defaultdict(lambda: defaultdict(list))
    unfiltered_seqs_dict = defaultdict(lambda: defaultdict(list))
    chimera_count_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    all_seqs_count = 0
    chimeric_seqs_count = 0
    non_chimeric_seqs_count = 0
    seqs_filtered_out_count = 0
    seqs_not_tested_count = 0

    for sample, genes_dict in sorted(combined_samples_dict.items()):

        for gene_name, paralog_dict in genes_dict.items():
            # if sample == 'Flindersia_brayleyana_MELUD105881A' and gene_name == '4471':
            #     print(f'paralog_dict: {paralog_dict}')

            # Initialise counts at zero:
            chimera_count_dict[gene_name][sample]['non_chimera_count'] = 0
            chimera_count_dict[gene_name][sample]['chimera_count'] = 0
            chimera_count_dict[gene_name][sample]['not_tested'] = 0

            for paralog_name, chimera_text in paralog_dict.items():
                all_seqs_count += 1

                # Recover the sequence regardless of whether it's chimeric or not
                seqrecord = sample_to_gene_seqrecord_dict[sample][gene_name][paralog_name]
                unfiltered_seqs_dict[gene_name][sample].append(seqrecord)

                if chimera_text in ['no_chimera_from_multi_hit', 'no_chimera_as_single_BLAT_hit']:

                    seqrecord = sample_to_gene_seqrecord_dict[sample][gene_name][paralog_name]
                    non_chimeric_seqs_dict[gene_name][sample].append(seqrecord)
                    chimera_count_dict[gene_name][sample]['non_chimera_count'] += 1
                    non_chimeric_seqs_count += 1

                elif chimera_text in ['no_detectable_chimera_as_single_BLAT_hit_after_length_filtering',
                                      'no_detectable_chimera_as_no_seqs_left_after_length_filtering',
                                      'no_chimera_test_performed_as_too_little_seq_length_after_length_filtering',
                                      'no_chimera_test_performed_as_too_few_contigs_after_length_filtering',
                                      'unknown_repeated_subject']:

                    chimera_count_dict[gene_name][sample]['not_tested'] += 1
                    seqs_filtered_out_count += 1
                    seqs_not_tested_count += 1
                else:
                    chimeric_seqs_count += 1
                    seqs_filtered_out_count += 1
                    chimera_count_dict[gene_name][sample]['chimera_count'] += 1

    logger.info(f'{"[INFO]:":10} Number of genes: {len(chimera_count_dict)}')
    logger.info(f'{"[INFO]:":10} Count of all sequences: {all_seqs_count}')
    logger.info(f'{"[INFO]:":10} Count of non-chimeric sequences: {non_chimeric_seqs_count}')
    logger.info(f'{"[INFO]:":10} Count of chimeric sequences: {chimeric_seqs_count}')
    logger.info(f'{"[INFO]:":10} Count of sequences not tested: {seqs_not_tested_count}')
    logger.info(f'{"[INFO]:":10} Count of filtered out sequences (chimeras and not tested): {seqs_filtered_out_count}')

    ####################################################################################################################
    # Write a *.tsv report for sample vs gene showing non-chimera vs chimera counts, and another report showing
    # percentage of chimeras (the latter for a heatmap):
    ####################################################################################################################
    sample_names = set()

    outfile_chimera_report = f'{outdir_reports}/chimera_report_min_samples_threshold_{min_samples_threshold}.tsv'
    outfile_chimera_report_heatmap = (f'{outdir_reports}/chimera_report_heatmap_data_min_samples_threshold'
                                      f'_{min_samples_threshold}.tsv')

    outfile_chimera_report_heatmap_showing_not_tested = \
        f'{outdir_reports}/chimera_report_heatmap_data_min_samples_threshold_{min_samples_threshold}.not_tested.tsv'

    logger.info(f'{"[INFO]:":10} Writing report: {outfile_chimera_report}')
    logger.info(f'{"[INFO]:":10} Writing report: {outfile_chimera_report_heatmap}')
    logger.info(f'{"[INFO]:":10} Writing report: {outfile_chimera_report_heatmap_showing_not_tested}')

    with (open(outfile_chimera_report, 'w') as chimera_report_handle,
          open(outfile_chimera_report_heatmap, 'w') as heatmap_data_handle,
          open(outfile_chimera_report_heatmap_showing_not_tested, 'w') as heatmap_data_showing_not_tested_handle):

        # Get sorted sample names for first row of report:
        for gene, sample_dict in sorted(chimera_count_dict.items()):
            for sample, count_dict in sorted(sample_dict.items()):
                sample_names.add(sample)

        sample_names_sorted = sorted(list(sample_names))
        sample_names_joined = '\t'.join(sample_names_sorted)

        chimera_report_handle.write(f'gene\t{sample_names_joined}\n')
        heatmap_data_handle.write(f'gene\t{sample_names_joined}\n')
        heatmap_data_showing_not_tested_handle.write(f'gene\t{sample_names_joined}\n')

        # Get gene stats:
        for gene, sample_dict in sorted(chimera_count_dict.items()):
            gene_row = [gene]
            gene_row_heatmap_data = [gene]
            gene_row_heatmap_data_showing_not_tested = [gene]

            for sample in sample_names_sorted:
                sample_count_dict = sample_dict[sample]

                if len(sample_count_dict) == 0:
                    logger.debug(f'No captus sequences for gene {gene} sample {sample}!')
                    count_values = 'None'
                    gene_row.append(count_values)
                    gene_row_heatmap_data.append('NaN')
                    gene_row_heatmap_data_showing_not_tested.append('NaN')
                else:
                    non_chimera_count = sample_count_dict['non_chimera_count']
                    chimera_count = sample_count_dict['chimera_count']

                    count_values = f'{non_chimera_count}/{chimera_count}'
                    gene_row.append(count_values)

                    total_count = non_chimera_count + chimera_count

                    try:
                        chimera_ratio = chimera_count / total_count
                        # print(f'chimera_ratio: {chimera_ratio}')
                        gene_row_heatmap_data.append(str(chimera_ratio))
                        gene_row_heatmap_data_showing_not_tested.append(str(chimera_ratio))
                    except ZeroDivisionError:
                        gene_row_heatmap_data.append('NaN')  # it will show in heatmap same as 'no seqs from captus'
                        gene_row_heatmap_data_showing_not_tested.append('2')

            gene_row_to_write = '\t'.join(gene_row)
            chimera_report_handle.write(f'{gene_row_to_write}\n')

            gene_row_heatmap_data_to_write = '\t'.join(gene_row_heatmap_data)
            heatmap_data_handle.write(f'{gene_row_heatmap_data_to_write}\n')

            gene_row_heatmap_data_to_write_showing_not_tested = '\t'.join(gene_row_heatmap_data_showing_not_tested)
            heatmap_data_showing_not_tested_handle.write(f'{gene_row_heatmap_data_to_write_showing_not_tested}\n')

    ####################################################################################################################
    # Write a heatmap for sample vs gene showing proportion of chimeras:
    ####################################################################################################################
    heatmap_name = f'{outdir_reports}/chimera_proportion_heatmap_captus_{min_samples_threshold}.png'
    logger.info(f'{"[INFO]:":10} Writing heatmap: {heatmap_name}')

    # Heatmap - read in the table created above:
    df = pd.read_csv(outfile_chimera_report_heatmap, delimiter='\t', )

    # Melt the dataframe to make it suitable for the df.pivot method:
    df = df.melt(id_vars=['gene'], var_name='sample_id', value_name='chimera_proportion')

    # Pivot the dataframe for input into the seaborn heatmap function:
    df = df.pivot(index='gene', columns='sample_id', values='chimera_proportion')

    # sns.set(rc={'figure.figsize': (10, 60)})
    sns.set(rc={'figure.figsize': (30, 90)})
    sns.set_style('ticks')
    cmap = 'plasma'

    heatmap = sns.heatmap(df,
                          vmin=0,
                          cmap=cmap,
                          xticklabels=1,
                          yticklabels=1,
                          square=True,
                          linewidth=0.5,
                          cbar=True,
                          cbar_kws={"orientation": "vertical",
                                    "pad": 0.01,
                                    "shrink": 0.25},
                          )

    heatmap.set_facecolor('lightgrey')  # Makes transparent NaN values to show up as grey squares

    heatmap.set_yticklabels(heatmap.get_yticklabels(), rotation=0)

    # Save the figure:
    plt.savefig(heatmap_name, dpi=300, bbox_inches='tight')

    ####################################################################################################################
    # Write a heatmap for sample vs gene showing proportion of chimeras, and also showing 'not tested':
    ####################################################################################################################
    plt.close()
    heatmap_name = f'{outdir_reports}/chimera_proportion_heatmap_captus_showing_not_tested_{min_samples_threshold}.png'
    logger.info(f'{"[INFO]:":10} Writing heatmap: {heatmap_name}')

    # Heatmap - read in the table created above:
    df = pd.read_csv(outfile_chimera_report_heatmap_showing_not_tested, delimiter='\t', )

    # Melt the dataframe to make it suitable for the df.pivot method:
    df = df.melt(id_vars=['gene'], var_name='sample_id', value_name='chimera_proportion')

    # Pivot the dataframe for input into the seaborn heatmap function:
    df = df.pivot(index='gene', columns='sample_id', values='chimera_proportion')

    # sns.set(rc={'figure.figsize': (10, 60)})
    sns.set(rc={'figure.figsize': (30, 90)})
    sns.set_style('ticks')
    # cmap = 'plasma'
    cmap = matplotlib.colormaps["plasma"]
    cmap.set_over("chartreuse")
    cmap.set_under("chartreuse")

    # Create a copy of the dataframe for color scaling
    df_for_scaling = df.copy()

    # Replace values of 2 with NaN so they don't affect the color scale
    df_for_scaling = df_for_scaling.replace(2, np.nan)

    # Calculate vmin and vmax based on the data excluding values of 2
    vmin = df_for_scaling.min().min()
    vmax = df_for_scaling.max().max()

    # If all values are 2, set default vmin and vmax
    if pd.isna(vmin) or pd.isna(vmax):
        vmin = 0
        vmax = 1

    heatmap = sns.heatmap(df,
                          vmin=vmin,
                          vmax=vmax,
                          cmap=cmap,
                          xticklabels=1,
                          yticklabels=1,
                          square=True,
                          linewidth=0.5,
                          cbar=True,
                          cbar_kws={"orientation": "vertical",
                                    "pad": 0.01,
                                    "shrink": 0.25},
                          )

    heatmap.set_facecolor('lightgrey')  # Makes transparent NaN values to show up as grey squares

    heatmap.set_yticklabels(heatmap.get_yticklabels(), rotation=0)

    # Save the figure:
    plt.savefig(heatmap_name, dpi=300, bbox_inches='tight')

    ####################################################################################################################
    # For each sample that contains 2 or more non-chimeric captus paralog sequences, recover the seqrecords. Note
    # that this is done for Harvey's Flindersia dataset as we _want_ paralogous genes in this case. 28 Feb 2025 -
    # changed from 2 to a default of 0, with user changeable option <min_num_paralogs_per_sample>
    ####################################################################################################################
    sample_to_gene_paralog_seqrecord_min_num_paralogs_for_sample_dict = defaultdict(lambda: defaultdict())
    single_seq_genes_dict = defaultdict(lambda: defaultdict())

    for gene_name, sample_dict in sorted(non_chimeric_seqs_dict.items()):

        # 28 November 2024:
        gene_number_of_seqs = []

        for sample, seqrecord_list in sample_dict.items():
            gene_number_of_seqs.append(len(seqrecord_list))

            if len(seqrecord_list) >= min_num_paralogs_per_sample:
                sample_to_gene_paralog_seqrecord_min_num_paralogs_for_sample_dict[gene_name][sample] = seqrecord_list

        number_of_samples = len(gene_number_of_seqs)
        samples_with_only_one_seq = 0
        for value in gene_number_of_seqs:
            if value == 1:
                samples_with_only_one_seq += 1

        ratio_of_single_seqs = samples_with_only_one_seq / number_of_samples
        if ratio_of_single_seqs >= 0.90:
            logger.info(
                f'{"[INFO]:":10} For gene {gene_name}, the number of samples with a single sequence, which is also '
                f'non-chimeric: {samples_with_only_one_seq}. The total number of samples is: {number_of_samples}. '
                f'The ratio of single sequences is: {ratio_of_single_seqs}')

            single_seq_genes_dict[gene_name]['samples_with_only_one_seq'] = samples_with_only_one_seq
            single_seq_genes_dict[gene_name]['ratio_of_single_seqs'] = ratio_of_single_seqs
            single_seq_genes_dict[gene_name]['number_of_samples'] = number_of_samples

    outfile_genes_with_one_seq_per_sample_in_90_percent_samples = \
        (f'{outdir_reports}/genes_with_one_seq_per_sample_in_90_percent_samples_min_samples_threshold'
         f'_{min_samples_threshold}.tsv')
    logger.info(f'{"[INFO]:":10} Writing report: {outfile_genes_with_one_seq_per_sample_in_90_percent_samples}')

    with open(outfile_genes_with_one_seq_per_sample_in_90_percent_samples,
              'w') as single_seq_handle:
        single_seq_handle.write(f'gene_name\t'
                                f'number_samples_with_single_non_chimeric_seq\t'
                                f'total_samples\t'
                                f'ratio_of_single_seqs\n')

        for gene_name, stats_dict in single_seq_genes_dict.items():
            single_seq_handle.write(f'{gene_name}\t'
                                    f'{stats_dict["samples_with_only_one_seq"]}\t'
                                    f'{stats_dict["number_of_samples"]}\t'
                                    f'{stats_dict["ratio_of_single_seqs"]}\n')

    ####################################################################################################################
    # For a given gene, if >= min_threshold samples have <min_num_paralogs_per_sample> or more non-chimeric captus
    # paralog sequences, recover the sequences:
    ####################################################################################################################
    non_chimeric_genes_and_paralogs_after_sample_threshold_filtering = defaultdict(lambda: defaultdict(int))
    genes_not_recovered_dict = dict()

    for gene, sample_dict in sample_to_gene_paralog_seqrecord_min_num_paralogs_for_sample_dict.items():

        gene_seqs_to_write = []

        if len(sample_dict) >= min_threshold:
            for sample, seqrecord_list in sample_dict.items():
                for seq in seqrecord_list:
                    gene_seqs_to_write.append(seq)
                    non_chimeric_genes_and_paralogs_after_sample_threshold_filtering[gene][sample] += 1

            gene_outfile = f'{outdir_non_chimeras}/{gene}.fasta'
            with open(gene_outfile, 'w') as gene_outfile_handle:
                SeqIO.write(gene_seqs_to_write, gene_outfile_handle, 'fasta')
        else:
            logger.debug(f'Not recovering sequences from any samples for gene {gene}, as there are fewer than'
                         f' {min_threshold} of {total_number_of_samples} samples containing '
                         f'{min_num_paralogs_per_sample} or more sequences!')

            genes_not_recovered_dict[gene] = len(sample_dict)

            for sample, seqrecord_list in sample_dict.items():
                non_chimeric_genes_and_paralogs_after_sample_threshold_filtering[gene][sample] = 0

        if len(genes_not_recovered_dict) != 0:
            with open(f'genes_not_recovered_{min_samples_threshold}.tsv', 'w') as genes_not_recovered_handle:
                genes_not_recovered_handle.write(f'min_samples_threshold_fraction\t'
                                                 f'total_number_of_samples\t'
                                                 f'min_threshold_number\n'
                                                 f'{min_samples_threshold}\t'
                                                 f'{total_number_of_samples}\t'
                                                 f'{min_threshold}\n')

                genes_not_recovered_handle.write(f'\n'
                                                 f'gene_name\t'
                                                 f'number_samples_with_paralogs'
                                                 f'\t'
                                                 f'\n')

                for gene_name, number_of_samples_with_paralogs in genes_not_recovered_dict.items():
                    genes_not_recovered_handle.write(f'{gene_name}\t'
                                                     f'{number_of_samples_with_paralogs}'
                                                     f'\t'
                                                     f'\n')

    ####################################################################################################################
    # Write a report for sample vs gene showing number of non-chimeric paralogs after sample threshold filtering:
    ####################################################################################################################
    outfile_non_chimeric_paralogs_report = \
        (f'{outdir_reports}/non_chimeric_paralogs_after_sample_threshold_filtering_report_min_samples_threshold'
         f'_{min_samples_threshold}.tsv')

    logger.info(f'{"[INFO]:":10} Writing report: {outfile_non_chimeric_paralogs_report}')

    with open(outfile_non_chimeric_paralogs_report, 'w') as non_chimeric_paralogs_report_handle:

        # Get sorted sample names for first row of report:
        for gene, sample_dict in sorted(chimera_count_dict.items()):  # contains all samples
            for sample, count_dict in sorted(sample_dict.items()):
                sample_names.add(sample)

        sample_names_sorted = sorted(list(sample_names))

        sample_names_joined = '\t'.join(sample_names_sorted)

        non_chimeric_paralogs_report_handle.write(f'gene\t{sample_names_joined}\n')

        # Get gene stats:
        for gene, sample_dict in sorted(non_chimeric_genes_and_paralogs_after_sample_threshold_filtering.items()):
            gene_row = [gene]

            for sample_name in sample_names_sorted:
                if sample_name not in sample_dict:
                    logger.warning(f'{"[WARNING]:":10} sample_name: {sample_name} not found for gene {gene}!')
                    gene_row.append('NaN')
                else:
                    paralog_count = sample_dict[sample_name]
                    if paralog_count == 0:
                        gene_row.append('NaN')
                    else:
                        gene_row.append(str(paralog_count))

            gene_row_to_write = '\t'.join(gene_row)
            non_chimeric_paralogs_report_handle.write(f'{gene_row_to_write}\n')

    ####################################################################################################################
    # Write a heatmap for sample vs gene showing number of non-chimeric paralogs after sample threshold filtering:
    ####################################################################################################################
    plt.close()

    outfile_heatmap_name = (f'{outdir_reports}/'
                            f'non_chimeric_paralog_count_heatmap_captus_after_filtering_min_samples_threshold'
                            f'_{min_samples_threshold}.png')

    logger.info(f'{"[INFO]:":10} Writing heatmap: {outfile_heatmap_name}')

    # Heatmap - read in the table created above:
    df = pd.read_csv(outfile_non_chimeric_paralogs_report, delimiter='\t', )

    # Melt the dataframe to make it suitable for the df.pivot method:
    df = df.melt(id_vars=['gene'], var_name='sample_id', value_name='non_chimeric_paralog_count')

    # Pivot the dataframe for input into the seaborn heatmap function:
    df = df.pivot(index='gene', columns='sample_id', values='non_chimeric_paralog_count')

    sns.set(rc={'figure.figsize': (30, 90)})
    sns.set_style('ticks')
    cmap = 'plasma'

    heatmap_paralogs = sns.heatmap(df,
                                   vmin=0,
                                   cmap=cmap,
                                   xticklabels=1,
                                   yticklabels=1,
                                   square=True,
                                   linewidth=0.5,
                                   cbar=True,
                                   cbar_kws={"orientation": "vertical",
                                             "pad": 0.01,
                                             "shrink": 0.25},
                                   )

    heatmap_paralogs.set_facecolor('lightgrey')  # Makes transparent NaN values to show up as grey squares

    heatmap_paralogs.set_yticklabels(heatmap_paralogs.get_yticklabels(), rotation=0)

    # Save the figure:
    plt.savefig(outfile_heatmap_name, dpi=300, bbox_inches='tight')

    ####################################################################################################################
    # Write a heatmap for sample vs gene showing number of non-chimeric paralogs after sample threshold filtering,
    # with any paralog count above 10 standardised to 10 [STEPPED COLOURS]:
    ####################################################################################################################
    plt.close()

    outfile_heatmap_name = \
        (f'{outdir_reports}/non_chimeric_paralog_count_heatmap_captus_after_filtering_max_count_10_stepped'
         f'_{min_samples_threshold}.png')

    logger.info(f'{"[INFO]:":10} Writing heatmap: {outfile_heatmap_name}')

    # Heatmap - read in the table created above:
    df = pd.read_csv(outfile_non_chimeric_paralogs_report, delimiter='\t', )
    not_null_mask = df.loc[:, df.columns[1]:].notna().any(axis=1)
    not_null_rows = df[not_null_mask]

    null_mask = df.loc[:, df.columns[1]:].isna().all(axis=1)
    null_rows = df[null_mask]
    df = pd.concat([null_rows, not_null_rows], axis=0)

    # For each paralog count, if the value is greater than 10, assign it to 10:
    df.where((df.loc[:, df.columns[1]:] < 10) | (df.loc[:, df.columns[1]:].isna()), 10, inplace=True)

    # Melt the dataframe to make it suitable for the df.pivot method:
    df = df.melt(id_vars=['gene'], var_name='sample_id', value_name='non_chimeric_paralog_count')

    # Pivot the dataframe for input into the seaborn heatmap function:
    df = df.pivot_table(index='gene',
                        columns='sample_id',
                        values='non_chimeric_paralog_count',
                        sort=False,
                        dropna=False)

    sns.set(rc={'figure.figsize': (30, 90)})
    sns.set_style('ticks')
    cmap = matplotlib.colors.ListedColormap(matplotlib.colormaps["Set3"].colors[:10])

    heatmap_paralogs = sns.heatmap(df,
                                   vmin=0,
                                   cmap=cmap,
                                   xticklabels=1,
                                   yticklabels=1,
                                   square=True,
                                   linewidth=0.5,
                                   cbar=True,
                                   cbar_kws={"orientation": "vertical",
                                             "pad": 0.01,
                                             "shrink": 0.25},
                                   )

    heatmap_paralogs.set_facecolor('lightgrey')  # Makes transparent NaN values to show up as grey squares

    heatmap_paralogs.set_yticklabels(heatmap_paralogs.get_yticklabels(), rotation=0)

    # Save the figure:
    plt.savefig(outfile_heatmap_name, dpi=300, bbox_inches='tight')

    ####################################################################################################################
    # Write a heatmap for sample vs gene showing number of non-chimeric paralogs after sample threshold filtering,
    # with any paralog count above 10 standardised to 10 [CONTINUOUS COLOURS]:
    ####################################################################################################################
    plt.close()

    outfile_heatmap_name = (f'{outdir_reports}/non_chimeric_paralog_count_heatmap_captus_after_filtering_max_'
                            f'count_10_continuous_{min_samples_threshold}.png')

    logger.info(f'{"[INFO]:":10} Writing heatmap: {outfile_heatmap_name}')

    # Heatmap - read in the table created above:
    df = pd.read_csv(outfile_non_chimeric_paralogs_report, delimiter='\t', )
    not_null_mask = df.loc[:, df.columns[1]:].notna().any(axis=1)
    not_null_rows = df[not_null_mask]

    null_mask = df.loc[:, df.columns[1]:].isna().all(axis=1)
    null_rows = df[null_mask]
    df = pd.concat([null_rows, not_null_rows], axis=0)

    # For each paralog count, if the value is greater than 10, assign it to 10:
    df.where((df.loc[:, df.columns[1]:] < 10) | (df.loc[:, df.columns[1]:].isna()), 10, inplace=True)

    # Melt the dataframe to make it suitable for the df.pivot method:
    df = df.melt(id_vars=['gene'], var_name='sample_id', value_name='non_chimeric_paralog_count')

    # Pivot the dataframe for input into the seaborn heatmap function:
    df = df.pivot_table(index='gene',
                        columns='sample_id',
                        values='non_chimeric_paralog_count',
                        sort=False,
                        dropna=False)

    sns.set(rc={'figure.figsize': (30, 90)})
    sns.set_style('ticks')
    cmap = 'summer'

    heatmap_paralogs = sns.heatmap(df,
                                   vmin=0,
                                   cmap=cmap,
                                   xticklabels=1,
                                   yticklabels=1,
                                   square=True,
                                   linewidth=0.5,
                                   cbar=True,
                                   cbar_kws={"orientation": "vertical", "pad": 0.01, "shrink": 0.25},
                                   )

    heatmap_paralogs.set_facecolor('lightgrey')  # Makes transparent NaN values to show up as grey squares

    heatmap_paralogs.set_yticklabels(heatmap_paralogs.get_yticklabels(), rotation=0)

    # Save the figure:
    plt.savefig(outfile_heatmap_name, dpi=300, bbox_inches='tight')

    ####################################################################################################################
    # Write fasta sequences for all captus seqs, regardless of whether they are chimeric or not, for comparisons:
    ####################################################################################################################
    for gene, sample_dict in unfiltered_seqs_dict.items():
        gene_seqs_to_write = []

        for sample, seqrecord_list in sample_dict.items():
            for seq in seqrecord_list:
                gene_seqs_to_write.append(seq)

        gene_outfile = f'{outdir_all_seqs}/{gene}.all.fasta'
        with open(gene_outfile, 'w') as gene_outfile_handle:
            SeqIO.write(gene_seqs_to_write, gene_outfile_handle, 'fasta')


def map_captus_stitched_contigs(sample_folder,
                                gene_list,
                                min_seq_length,
                                min_length_percentage,
                                min_contig_number_percentage,
                                outdir,
                                threads=1):
    """
    Process a single sample folder to map Captus stitched contigs to a reference genome.

    This function processes a single sample folder containing Captus stitched contigs,
    maps them to a reference genome, and identifies chimeric sequences based on mapping patterns.
    It extracts coding regions from contigs, filters sequences based on length criteria,
    and determines if sequences are chimeric based on their mapping behavior.

    Args:
        sample_folder (str): Path to the sample folder containing Captus assembly data
        gene_list (list): List of gene IDs to process
        min_seq_length (int): Minimum sequence length threshold for filtering
        min_length_percentage (float): Minimum percentage of the paralog sequence length remaining after
                                      filtering contig hits via <min_seq_length> for chimera detection to be performed.
                                      Otherwise, no sequence is retained for this sample/gene. Defaults to 0.80.
        min_contig_number_percentage (float): Minimum percentage of the total number of contig hits remaining for a
                                             given paralog after filtering contig hits via <min_seq_length> for
                                             chimera detection to be performed. Otherwise, no sequence is retained
                                             for this sample/gene. Defaults to 0.80.
        outdir (str): Output directory for storing mapping results
        threads (int, optional): Number of threads to use for mapping. Defaults to 1.

    Returns:
        tuple: A tuple containing:
            - str: Sample name extracted from the folder path
            - dict: Dictionary mapping gene names to dictionaries of paralog results
                  Each paralog result indicates whether it's chimeric or not
    """

    sample_basename = os.path.basename(sample_folder)
    sample_name = sample_basename.split('__')[0]
    logger.info(f'{"[INFO]:":10} Processing Captus sample: {sample_name}')

    # Create sample dict for reports:
    sample_dict = dict()

    ####################################################################################################################
    # Parse <sample_name>_recovery_stats.tsv to get hit contigs names, strands and coordinates (recovers paralogs too):
    ####################################################################################################################
    gene_stats_dict = defaultdict(list)
    with open(f'{sample_folder}/06_assembly_annotated/{sample_name}_recovery_stats.tsv') as recovery_stats_handle:
        lines = iter(recovery_stats_handle.readlines())
        header = next(lines)
        for line in lines:
            gene = line.split('\t')[2]

            if gene not in gene_list:  # DEBUG !!
                continue

            megahit_contigs = [item.strip() for item in line.split('\t')[22].split(';')]
            megahit_contigs_strands = [item.strip() for item in line.split('\t')[23].split(';')]
            megahit_contigs_coords = [item.strip() for item in line.split('\t')[24].split(';')]

            zipped_stats = list(zip(megahit_contigs,
                                    megahit_contigs_strands,
                                    megahit_contigs_coords))
            gene_stats_dict[gene].append(zipped_stats)

    ####################################################################################################################
    # Get coding regions from megahit contigs and add to dict:
    ####################################################################################################################
    all_raw_megahit_contigs_dict = SeqIO.to_dict(SeqIO.parse(f'{sample_folder}/06_assembly_annotated/'
                                                             f'{sample_name}_hit_contigs.fasta', 'fasta'))

    gene_seq_records_dict = defaultdict(lambda: defaultdict(list))

    for gene, paralog_stats_list in gene_stats_dict.items():
        for count, paralog_stats in enumerate(paralog_stats_list, 1):  # iterate over each paralog

            for node_stats_tuple in paralog_stats:
                node_name, strand, coords = node_stats_tuple
                coords_list = coords.split(',')
                raw_megahit_contig = all_raw_megahit_contigs_dict[node_name]

                if strand == '+':
                    exon_regions = ''

                    for coord in coords_list:
                        start, end = coord.split('-')
                        exon_seq = raw_megahit_contig.seq[int(start) - 1: int(end)]
                        exon_regions = exon_regions + exon_seq

                    seqrecord = SeqRecord(seq=exon_regions,
                                          name=raw_megahit_contig.name,
                                          id=raw_megahit_contig.id,
                                          description=raw_megahit_contig.description)

                    gene_seq_records_dict[gene][f'paralog_{count}'].append(seqrecord)  # name the paralog here

                else:  # hit on negative strand
                    exon_regions = ''

                    for coord in coords_list:
                        start, end = coord.split('-')
                        exon_seq = raw_megahit_contig.seq[int(start) - 1: int(end)]
                        exon_seq = exon_seq.reverse_complement()
                        exon_regions = exon_regions + exon_seq

                    seqrecord = SeqRecord(seq=exon_regions,
                                          name=raw_megahit_contig.name,
                                          id=raw_megahit_contig.id,
                                          description=raw_megahit_contig.description)

                    gene_seq_records_dict[gene][f'paralog_{count}'].append(seqrecord)  # name the paralog here

    ####################################################################################################################
    # Filter the seqrecords in the dict:
    ####################################################################################################################
    for gene, paralogs_dict in gene_seq_records_dict.items():
        for paralog, node_list in paralogs_dict.items():  # process a single paralog

            node_lengths = [len(node) for node in node_list]
            total_paralog_length = sum(node_lengths)
            number_of_nodes_greater_than_min_seq_length = 0
            length_of_nodes_greater_than_min_seq_length = 0

            ############################################################################################################
            # Determine some stats for filtering and thresholds for performing a chimera test or not. More than
            # <min_length_percentage> of the length of the paralog and <min_contig_number_percentage> of the number of
            # contigs for the paralog must be remaining after filtering out short contigs below <min_seq_length>:
            ############################################################################################################
            for node_length in node_lengths:
                if node_length >= min_seq_length:
                    number_of_nodes_greater_than_min_seq_length += 1
                    length_of_nodes_greater_than_min_seq_length += node_length

            if length_of_nodes_greater_than_min_seq_length / total_paralog_length >= min_length_percentage:
                passes_length_cutoff_after_min_seq_length_filtering = True
            else:
                passes_length_cutoff_after_min_seq_length_filtering = False

            if number_of_nodes_greater_than_min_seq_length / len(node_list) >= min_contig_number_percentage:
                passes_contig_number_cutoff_after_min_seq_length_filtering = True
            else:
                passes_contig_number_cutoff_after_min_seq_length_filtering = False

            if len(node_list) == 1:
                logger.debug(f'Only one seq for gene {gene} paralog {paralog} sample {sample_name}, so no chimera from '
                             f'stitching.')
                if gene not in sample_dict:
                    sample_dict[gene] = {}

                sample_dict[gene][paralog] = 'no_chimera_as_single_BLAT_hit'

                outdir_single_hit = createfolder(f'{outdir}/{sample_name}/single_hit')

                # Write to file:
                with open(f'{outdir_single_hit}/{gene}_{paralog}_single_scipio_hit.fasta', 'w') as scipio_single_handle:
                    SeqIO.write(node_list[0], scipio_single_handle, 'fasta')

            elif len(node_list) > 1 and number_of_nodes_greater_than_min_seq_length == 1:
                logger.debug(f'Only one seq for gene {gene} paralog {paralog} sample {sample_name} after filtering '
                             f'out seqs < {min_seq_length}, so no _detectable_ chimera from stitching.')

                if gene not in sample_dict:
                    sample_dict[gene] = {}

                sample_dict[gene][paralog] = 'no_detectable_chimera_as_single_BLAT_hit_after_length_filtering'

            elif len(node_list) > 1 and number_of_nodes_greater_than_min_seq_length == 0:
                logger.debug(f'No seqs left for gene {gene} paralog {paralog} sample {sample_name} after filtering out '
                             f'seqs < {min_seq_length}, so no _detectable_ chimera from stitching.')

                if gene not in sample_dict:
                    sample_dict[gene] = {}

                sample_dict[gene][paralog] = 'no_detectable_chimera_as_no_seqs_left_after_length_filtering'

            elif not passes_length_cutoff_after_min_seq_length_filtering:
                logger.debug(f'Too little of the paralog sequence remaining for gene {gene} paralog {paralog} sample'
                             f' {sample_name} after filtering out seqs < {min_seq_length}, so no chimera detection '
                             f'performed.')

                if gene not in sample_dict:
                    sample_dict[gene] = {}

                sample_dict[gene][paralog] = 'no_chimera_test_performed_as_too_little_seq_length_after_length_filtering'

            elif not passes_contig_number_cutoff_after_min_seq_length_filtering:
                logger.debug(f'Too few contigs remaining for gene {gene} paralog {paralog} sample {sample_name} after '
                             f'filtering out seqs < {min_seq_length}, so no chimera detection performed.')

                if gene not in sample_dict:
                    sample_dict[gene] = {}

                sample_dict[gene][paralog] = 'no_chimera_test_performed_as_too_few_contigs_after_length_filtering'

            else:
                seqs_to_write = [node for node in node_list]

                # Create an output dir for mapping results for this sample/gene:
                outdir_gene = createfolder(f'{outdir}/{sample_name}/{gene}')
                outfile_gene = f'{outdir_gene}/{paralog}_scipio_hits.fasta'

                # Strip any N's from the seqs:
                for seq in seqs_to_write:
                    seq.seq = seq.seq.replace('N', '')

                # Filter any seqs less than <min_seq_length>:
                seqs_to_write_filtered = []
                for seq in seqs_to_write:
                    if len(seq) < min_seq_length:
                        logger.debug(f'Sample {sample_name} gene {gene} paralog {paralog} sequence {seq.name} is less '
                                     f'than min_seq_length = {min_seq_length} bases long, sequence will not be written '
                                     f'or mapped.')
                    else:
                        seqs_to_write_filtered.append(seq)

                # Write seqrecords to file:
                if len(seqs_to_write) > 0:
                    with open(outfile_gene, 'w') as outfile_gene_handle:
                        SeqIO.write(seqs_to_write_filtered, outfile_gene_handle, 'fasta')
                else:
                    raise ValueError(f'No seqs to write for sample {sample_name} gene {gene}')

                ########################################################################################################
                # Map the sequences to the genome ref with BBmap.sh (mapPacBio.sh):
                ########################################################################################################
                expected_alignment = f'{outdir}/{sample_name}/{gene}/{paralog}_mapping.sam'

                try:
                    assert file_exists_and_not_empty(expected_alignment)
                    logger.debug(f'Expected file {expected_alignment} exists, skipping...')

                except AssertionError:

                    try:
                        command = (f'mapPacBio.sh '
                                   f'Xmx=100g '
                                   f'in={outfile_gene} '
                                   f'k=13 '
                                   f'fastareadlen=6000 '
                                   f'maxindel=100000 '
                                   f'minid=0.76 '
                                   f'local=f '
                                   f'ambiguous=all '
                                   f'noheader=t '
                                   f'threads={threads} '
                                   f'out={expected_alignment}')

                        result = subprocess.run(command,
                                                universal_newlines=True,
                                                stdout=subprocess.PIPE,
                                                stderr=subprocess.PIPE,
                                                check=True,
                                                shell=True)

                        logger.debug(f'mapPacBio.sh check_returncode() is: {result.check_returncode()}')
                        logger.debug(f'mapPacBio.sh stdout is: {result.stdout}')
                        logger.debug(f'mapPacBio.sh stderr is: {result.stderr}')

                    except subprocess.CalledProcessError as exc:
                        logger.error(f'mapPacBio.sh FAILED. Output is: {exc}')
                        logger.error(f'mapPacBio.sh stdout is: {exc.stdout}')
                        logger.error(f'mapPacBio.sh stderr is: {exc.stderr}')

                        raise ValueError('There was an issue running mapPacBio.sh. Check input files!')

                assert file_exists_and_not_empty(expected_alignment)

                ########################################################################################################
                # Parse the mapping results:
                ########################################################################################################
                gene_dict = defaultdict(list)

                with open(expected_alignment) as mapping_sam_handle:
                    lines = mapping_sam_handle.readlines()
                    for line in lines:
                        stats = line.split('\t')
                        query_name = stats[0].split(',')[0]
                        sam_flag = stats[1]
                        subject_name = stats[2]
                        mapping_pos = stats[3]

                        gene_dict[query_name].append(subject_name)

                ########################################################################################################
                # Get the number of times each query occurs in the "{paralog}_scipio_hits.fasta" file:
                ########################################################################################################
                scipio_hits_names = [seq.name for seq in SeqIO.parse(outfile_gene, 'fasta')]
                scipio_hits_names_counter = Counter(scipio_hits_names)

                ########################################################################################################
                # Check that there is a mapping in the *.sam file for each sequence in the input paralog contig file:
                ########################################################################################################
                for name in scipio_hits_names:
                    try:
                        assert name in gene_dict
                    except AssertionError:
                        logger.warning(f'{"[WARNING]:":10} No mapping result in file {expected_alignment} for Scipio '
                                       f'hit {name}.')
                        raise

                ########################################################################################################
                # Sort dictionary by number of subject hits:
                gene_dict_sorted = {key: value for key, value in sorted(gene_dict.items(), key=lambda x: len(x[1]))}
                ########################################################################################################

                first_dict_entry = dict(itertools.islice(gene_dict_sorted.items(), 1))
                remaining_dict_entries = dict(itertools.islice(gene_dict_sorted.items(), 1, len(gene_dict_sorted) + 1))

                # Set some defaults:
                chimera_found = False
                unknown_repeated_subject = False

                ########################################################################################################
                # Iterate over the subject hits for the first query (i.e. smallest number of subject hits):
                ########################################################################################################
                for query, subject_list in first_dict_entry.items():
                    query_count = scipio_hits_names_counter[query]
                    assert query_count != 0
                    if query_count > 1:
                        print(f'query_count: {query_count}')

                    ####################################################################################################
                    # Check if there are too many subject hits (e.g. due to multiple copies of the gene on the
                    # genomic contig):
                    ####################################################################################################
                    first_dict_entry_counter = Counter(subject_list)
                    for subject_contig_name, count in first_dict_entry_counter.items():
                        if count > query_count:
                            logger.debug(
                                f'{"[WARNING]:":10} Sample {sample_name}, gene {gene}: subject {subject_contig_name} '
                                f'occurs more ({count}) than the number of times the query {query} occurs in the '
                                f'"{paralog}_scipio_hits.fasta" file ({query_count}). Sequence will be flagged as '
                                f'"unknown_repeated_subject".')

                            unknown_repeated_subject = True

                    for subject in subject_list:
                        for remaining_query, remaining_subject_list in remaining_dict_entries.items():
                            remaining_query_count = scipio_hits_names_counter[remaining_query]
                            assert remaining_query_count != 0
                            if remaining_query_count > 1:
                                print(f'remaining_query_count: {remaining_query_count}')
                            remaining_dict_entries_counter = Counter(remaining_subject_list)

                            if subject not in remaining_subject_list:
                                chimera_found = True
                            else:
                                remaining_dict_entries_subject_count = remaining_dict_entries_counter[subject]
                                assert remaining_dict_entries_subject_count != 0

                                if remaining_dict_entries_subject_count > remaining_query_count:
                                    logger.debug(
                                        f'{"[WARNING]:":10} Sample {sample_name}, gene {gene}: subject {subject} '
                                        f'occurs more ({remaining_dict_entries_subject_count}) than the number of '
                                        f'times the query {remaining_query} occurs in the "'
                                        f'{paralog}_scipio_hits.fasta" file  ({remaining_query_count}). Sequence '
                                        f'will be flagged as "unknown_repeated_subject".')

                                    unknown_repeated_subject = True

                if gene not in sample_dict:
                    sample_dict[gene] = {}

                if not unknown_repeated_subject:
                    if chimera_found:
                        sample_dict[gene][paralog] = 'chimera_from_multi_hit'
                    else:
                        # INSERT CHECK THAT HIT SAM FLAGS ARE THE SAME FOR A GIVEN NODE, AND NODE GENOMIC LOCATIONS
                        # ARE IN SAME ORDER AS gene_dict().
                        sample_dict[gene][paralog] = 'no_chimera_from_multi_hit'
                else:
                    sample_dict[gene][paralog] = 'unknown_repeated_subject'

    return sample_name, sample_dict


########################################################################################################################
# argparse options:
########################################################################################################################
def parse_arguments():
    """
    Parse command line arguments for the script.

    This function sets up the command line argument parser and defines all the
    required and optional arguments for the script.

    Returns:
        argparse.Namespace: An object containing the parsed command line arguments
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    group_1 = parser.add_mutually_exclusive_group(required=False)
    group_1.add_argument('--version', '-v',
                         dest='version',
                         action='version',
                         version=f'filter_chimeras_with_genome_ref_captus.py version {__version__}',
                         help='Print the script version number.')
    parser.add_argument('reference_genome_fasta',
                        type=str,
                        help='High quality reference genome for mapping target capture contigs from Captus.')
    parser.add_argument('target_file_fasta',
                        type=str,
                        help='The target file used for the Captus runs.')
    parser.add_argument('captus_sample_parent_folder',
                        type=str,
                        help='Parent folder containing Captus sample folders.')
    parser.add_argument('--min_seq_length',
                        type=int,
                        default=75,
                        help='The minimum sequence length to map to the genome. Default is: %(default)s')
    parser.add_argument('--min_length_percentage',
                        type=float,
                        default=0.80,
                        help='The minimum percentage of the paralog sequence length remaining after '
                             'filtering contig hits via <min_seq_length> for chimera detection to be performed.'
                             'Otherwise, no sequence is retained for this sample/gene.Default is: %(default)s')
    parser.add_argument('--min_contig_number_percentage',
                        type=float,
                        default=0.80,
                        help='The minimum percentage of the total number of contig hits remaining for a given paralog '
                             'after filtering contig hits via <min_seq_length> for chimera detection to be performed.'
                             'Otherwise, no sequence is retained for this sample/gene. Default is: %(default)s')
    parser.add_argument('--min_samples_threshold',
                        type=float,
                        default=0.75,
                        help='For a given gene, the minimum percentage of total samples to have '
                             '>= <min_num_paralogs_per_sample> non-chimeric sequences for sequences to be written to '
                             'file. Useful to increase stringency of alignment sample occupancy, as the cost of fewer '
                             'gene alignments. Default is: %(default)s')
    parser.add_argument('--min_num_paralogs_per_sample',
                        type=int,
                        default=0,
                        help='For a given gene for a given sample, the minimum number of non-chimeric paralog '
                             'sequences recovered for the sequences to be written to file. This can be useful when you '
                             'expect paralogs to be present (e.g. due to polyploidy), and want to skip genes that '
                             'might have hidden paralogy due to sequence or assembly issues. Default is: %(default)s')
    parser.add_argument('--pool',
                        type=int,
                        default=1,
                        help='Number of samples to run concurrently. Default is: %(default)s')
    parser.add_argument('--threads',
                        type=int,
                        default=1,
                        help='Number of threads to use for mapping for a given sample. Default is: %(default)s')

    results = parser.parse_args()

    return results


########################################################################################################################
# Function 'main':
########################################################################################################################
def main():
    """
    Main function that orchestrates the execution of the script.

    This function parses command line arguments, sets up output directories,
    and calls the main processing functions in the correct order.

    Returns:
        int: Exit code (0 for success, non-zero for failure)
    """
    args = parse_arguments()

    ####################################################################################################################
    # Create a dictionary from the argparse Namespace and print to screen:
    ####################################################################################################################
    parameters = vars(args)

    logger.info(f'{"[INFO]:":10} Script {__file__} version {__version__} was called with these arguments:\n')

    for parameter, value in parameters.items():
        logger.info(f'{" " * 10} {parameter}: {value}')
    logger.info('')

    ####################################################################################################################
    # Set base output folder names:
    ####################################################################################################################
    folder_00 = f'{cwd}/00_reports'
    folder_01 = f'{cwd}/01_captus_seqs_mapped'
    folder_02 = f'{cwd}/02_non_chimeric_captus_seqs_{args.min_samples_threshold}'
    folder_03 = f'{cwd}/03_all_captus_seqs'

    ####################################################################################################################
    # main functions:
    ####################################################################################################################
    create_genome_mapping_reference(args.reference_genome_fasta)

    gene_list = parse_target_file(args.target_file_fasta)

    map_captus_stitch_contigs_to_genome_mp(gene_list,
                                           args.reference_genome_fasta,
                                           args.captus_sample_parent_folder,
                                           min_seq_length=args.min_seq_length,
                                           min_samples_threshold=args.min_samples_threshold,
                                           min_num_paralogs_per_sample=args.min_num_paralogs_per_sample,
                                           outdir_reports=folder_00,
                                           outdir_captus_seqs_mapped=folder_01,
                                           outdir_non_chimeras=folder_02,
                                           outdir_all_seqs=folder_03,
                                           min_length_percentage=args.min_length_percentage,
                                           min_contig_number_percentage=args.min_contig_number_percentage,
                                           pool=args.pool,
                                           threads=args.threads)

    logger.info(f'\n{"*" * 90}\n')
    logger.info(f'{"[INFO]:":10} Script finished!')

    return 0


########################################################################################################################
# Run script:
########################################################################################################################

if __name__ == '__main__':

    if not len(sys.argv) >= 1:
        print(__doc__)
        sys.exit()

    sys.exit(main())

########################################################################################################################
########################################################################################################################
