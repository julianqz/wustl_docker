"""
External application wrappers
"""
# Info
__author__    = 'Jason Anthony Vander Heiden, Namita Gupta'
from presto import __version__, __date__

# Imports
import csv
import os
import re
import tempfile
import pandas as pd
from io import StringIO
from itertools import groupby
from subprocess import CalledProcessError, check_output, PIPE, Popen, STDOUT
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Presto imports
from presto.Defaults import default_muscle_exec, default_usearch_exec, \
                            default_blastn_exec, default_blastdb_exec, \
                            default_cdhit_exec
from presto.IO import readReferenceFile, printError, printWarning

# Defaults
default_cluster_ident = 0.9
default_length_ratio = 0.0
default_align_ident = 0.5
default_evalue = 1e-5
default_max_hits = 100
default_max_memory = 6000 #*JZ increased from 3000 to 6000


def runMuscle(seq_list, aligner_exec=default_muscle_exec):
    """
    Multiple aligns a set of sequences using MUSCLE

    Arguments:
      seq_list : a list of SeqRecord objects to align
      aligner_exec : the MUSCLE executable

    Returns:
      Bio.Align.MultipleSeqAlignment : Multiple alignment results.
    """
    # Return sequence if only one sequence in seq_list
    if len(seq_list) < 2:
        align = MultipleSeqAlignment(seq_list)
        return align

    # Set MUSCLE command
    cmd = [aligner_exec, '-diags', '-maxiters', '2']

    # Convert sequences to FASTA and write to string
    stdin_handle = StringIO()
    SeqIO.write(seq_list, stdin_handle, 'fasta-2line')
    stdin_str = stdin_handle.getvalue()
    stdin_handle.close()

    # Open MUSCLE process
    child = Popen(cmd, bufsize=-1, stdin=PIPE, stdout=PIPE, stderr=PIPE,
                  universal_newlines=True)

    # Send sequences to MUSCLE stdin and retrieve stdout, stderr
    stdout_str, __ = child.communicate(stdin_str)

    # Capture sequences from MUSCLE stdout
    stdout_handle = StringIO(stdout_str)
    align = AlignIO.read(stdout_handle, 'fasta')
    stdout_handle.close()

    return align


def runUClust(seq_list, ident=default_cluster_ident, length_ratio=default_length_ratio,
              seq_start=0, seq_end=None,
              threads=1, cluster_exec=default_usearch_exec):
    """
    Cluster a set of sequences using the UCLUST algorithm from USEARCH

    Arguments:
      seq_list (list): a list of SeqRecord objects to align.
      ident (float): the sequence identity cutoff to be passed to usearch.
      length_ratio (float): usearch parameter defining the minimum short/long length
                            ratio allowed within a cluster.
      seq_start (int): the start position to trim sequences at before clustering.
      seq_end (int): the end position to trim sequences at before clustering.
      threads (int): number of threads for usearch.
      cluster_exec (str): the path to the usearch executable.

    Returns:
      dict: {cluster id: list of sequence ids}.
    """
    # Function to trim and mask sequences
    gap_trans = str.maketrans({'-': 'N', '.': 'N'})
    def _clean(rec, i, j):
        seq = str(rec.seq[i:j])
        seq = seq.translate(gap_trans)
        return SeqRecord(Seq(seq), id=rec.id, name=rec.name, description=rec.description)

    # Make a trimmed and masked copy of each sequence so we don't mess up originals
    seq_trimmed = [_clean(x, seq_start, seq_end) for x in seq_list]

    # Return sequence if only one sequence in seq_iter
    if len(seq_trimmed) < 2:
        return {1:[seq_trimmed[0].id]}

    # If there are any empty sequences after trimming return None
    if any([len(x.seq) == 0 for x in seq_trimmed]):
        return None

    # Open temporary files
    in_handle = tempfile.NamedTemporaryFile(mode='w+t', encoding='utf-8')
    out_handle = tempfile.NamedTemporaryFile(mode='w+t', encoding='utf-8')

    # Define usearch command
    cmd = [cluster_exec,
           '-cluster_fast', in_handle.name,
           '-uc', out_handle.name,
           '-id', str(ident),
           '-minsl', str(length_ratio),
           '-qmask', 'none',
           '-minseqlength', '1',
           '-threads', str(threads)]

    # Write usearch input fasta file
    SeqIO.write(seq_trimmed, in_handle, 'fasta-2line')
    in_handle.seek(0)

    # Run usearch uclust algorithm
    try:
        stdout_str = check_output(cmd, stderr=STDOUT, shell=False,
                                  universal_newlines=True)
    except CalledProcessError as e:
        printError('Running command: %s\n%s' % (' '.join(cmd), e.output))

    # Parse the results of usearch
    # Output columns for the usearch 'uc' output format
    #   0 = entry type -- S: centroid seq, H: hit, C: cluster record (redundant with S)
    #   1 = group the sequence is assigned to
    #   8 = the id of the sequence
    #   9 = id of the centroid for cluster
    cluster_dict = {}
    for row in csv.reader(out_handle, delimiter='\t'):
        if row[0] in ('H', 'S'):
            # Trim sequence label to portion before space for usearch v9 compatibility
            key = int(row[1]) + 1
            # Trim sequence label to portion before space for usearch v9 compatibility
            hit = row[8].split()[0]
            # Update cluster dictionary
            cluster = cluster_dict.setdefault(key, [])
            cluster.append(hit)

    return cluster_dict if cluster_dict else None


def runCDHit(seq_list, ident=default_cluster_ident, length_ratio=default_length_ratio,
             seq_start=0, seq_end=None, max_memory=default_max_memory,
             threads=1, cluster_exec=default_cdhit_exec):
    """
    Cluster a set of sequences using CD-HIT

    Arguments:
      seq_list (list): a list of SeqRecord objects to align.
      ident (float): the sequence identity cutoff to be passed to cd-hit-est.
      length_ratio (float): cd-hit-est parameter defining the minimum short/long length
                            ratio allowed within a cluster.
      seq_start (int): the start position to trim sequences at before clustering.
      seq_end (int): the end position to trim sequences at before clustering.
      max_memory (int): cd-hit-est max memory limit (Mb)
      threads (int): number of threads for cd-hit-est.
      cluster_exec (str): the path to the cd-hit-est executable.

    Returns:
      dict: {cluster id: list of sequence ids}.
    """
    # Function to trim and mask sequences
    gap_trans = str.maketrans({'-': 'N', '.': 'N'})
    def _clean(rec, i, j):
        seq = str(rec.seq[i:j])
        seq = seq.translate(gap_trans)
        return SeqRecord(Seq(seq), id=rec.id, name=rec.name, description=rec.description)

    # Make a trimmed and masked copy of each sequence so we don't mess up originals
    seq_trimmed = [_clean(x, seq_start, seq_end) for x in seq_list]

    # Return sequence if only one sequence in seq_iter
    if len(seq_trimmed) < 2:
        return {1:[seq_trimmed[0].id]}

    # If there are any empty sequences after trimming return None
    if any([len(x.seq) == 0 for x in seq_trimmed]):
        return None

    # Open temporary files
    in_handle = tempfile.NamedTemporaryFile(mode='w+t', encoding='utf-8')
    out_handle = tempfile.NamedTemporaryFile(mode='w+t', encoding='utf-8')

    # Define usearch command
    cmd = [cluster_exec,
           '-i', in_handle.name,
           '-o', out_handle.name,
           '-c', str(ident),
           '-s', str(length_ratio),
           '-n', '3',
           '-d', '0',
           '-M', str(max_memory),
           '-T', str(threads)]

    # Write usearch input fasta file
    SeqIO.write(seq_trimmed, in_handle, 'fasta-2line')
    in_handle.seek(0)

    # Run CD-HIT
    try:
        stdout_str = check_output(cmd, stderr=STDOUT, shell=False,
                                  universal_newlines=True)
    except CalledProcessError as e:
        printError('Running command: %s\n%s' % (' '.join(cmd), e.output))

    # Parse the results of CD-HIT
    # Output of the .clstr file
    #   >Cluster 0
    #   0	17nt, >S01|BARCODE=CTAAGTGACTGGAGTTC... *
    #   1	17nt, >S02|BARCODE=CTAAGTGACTGGAGTTC... at +/100.00%
    #   2	17nt, >S07|BARCODE=CTAAGTGACTGGACTTC... at +/94.12%
    #   >Cluster 1
    #   0	17nt, >S12|BARCODE=TTTTTTTTTTTTTTTTT... *
    # Parsing regex
    block_regex = re.compile('>Cluster [0-9]+')
    id_regex = re.compile('([0-9]+\t[0-9]+nt, \>)(.+)(\.\.\.)')

    # Parse .clstr file
    cluster_dict = {}
    cluster_file = '%s.clstr' % out_handle.name
    with open(cluster_file, 'r') as cluster_handle:
        # Define parsing blocks
        clusters = groupby(cluster_handle, key=lambda x: block_regex.match(x))
        # Iterate over clusters and update return dict
        count = 1
        for key, group in clusters:
            if key is not None:
                __, block = next(clusters)
                cluster_dict[count] = [id_regex.match(x).group(2) for x in block]
                count += 1

    # Delete temp file
    os.remove(cluster_file)

    return cluster_dict if cluster_dict else None


def makeUBlastDb(ref_file, db_exec=default_usearch_exec):
    """
    Makes a ublast database file

    Arguments:
      ref_file : path to the reference database file.
      db_exec : path to the usearch executable.

    Returns:
      tuple : (location of the database, handle of the tempfile.NamedTemporaryFile)
    """
    # Open temporary files
    seq_handle = tempfile.NamedTemporaryFile(suffix='.fasta', mode='w+t', encoding='utf-8')
    db_handle = tempfile.NamedTemporaryFile(suffix='.udb')

    # Write temporary ungapped reference file
    ref_dict = readReferenceFile(ref_file)
    SeqIO.write(ref_dict.values(), seq_handle, format='fasta-2line')
    seq_handle.seek(0)

    # Define usearch command
    cmd = [db_exec,
           '-makeudb_ublast', seq_handle.name,
           '-wordlength', '9',
           '-output', db_handle.name,
           '-dbmask', 'none']
    try:
        stdout_str = check_output(cmd, stderr=STDOUT, shell=False,
                                  universal_newlines=True)
    except:
        seq_handle.close()
        printError('Failed to make usearch database.')

    # Close temporary sequence file
    seq_handle.close()

    return (db_handle.name, db_handle)


def makeBlastnDb(ref_file, db_exec=default_blastdb_exec):
    """
    Makes a blastn database file

    Arguments:
      ref_file : the path to the reference database file
      db_exec : the path to the makeblastdb executable

    Returns:
      tuple : (name and location of the database, handle of the tempfile.TemporaryDirectory)
    """
    # Open temporary files
    seq_handle = tempfile.NamedTemporaryFile(suffix='.fasta', mode='w+t', encoding='utf-8')
    db_handle = tempfile.TemporaryDirectory()

    # Write temporary ungapped reference file
    ref_dict = readReferenceFile(ref_file)
    SeqIO.write(ref_dict.values(), seq_handle, format='fasta-2line')
    seq_handle.seek(0)

    # Define usearch command
    cmd = [db_exec,
           '-in', seq_handle.name,
           '-out', os.path.join(db_handle.name, 'reference'),
           '-dbtype', 'nucl',
           '-title', 'reference',
           '-parse_seqids']
    try:
        stdout_str = check_output(cmd, stderr=STDOUT, shell=False,
                                  universal_newlines=True)
    except:
        seq_handle.close()
        printError('Failed to make blastn database.')

    # Close temporary sequence file
    seq_handle.close()

    return (os.path.join(db_handle.name, 'reference'), db_handle)


def runUBlast(seq, database, evalue=default_evalue, max_hits=default_max_hits,
              aligner_exec=default_usearch_exec):
    """
    Aligns a sequence against a reference database using the usearch_local algorithm of USEARCH

    Arguments:
      seq : a list of SeqRecord objects to align.
      database : the path to the ublast database or a fasta file.
      evalue : the E-value cut-off.
      maxhits : the maximum number of hits returned.
      aligner_exec : the path to the usearch executable.

    Returns:
      pandas.DataFrame : Alignment results.
    """
    # Open temporary files
    in_handle = tempfile.NamedTemporaryFile(mode='w+t', encoding='utf-8')
    out_handle = tempfile.NamedTemporaryFile(mode='w+t', encoding='utf-8')

    # Define usearch command
    cmd = [aligner_exec,
           '-ublast', in_handle.name,
           '-db', database,
           '-strand', 'plus',
           '-evalue', str(evalue),
           '-maxhits', str(max_hits),
           '-wordlength', '9',
           '-maxaccepts', '0',
           '-maxrejects', '0',
           '-userout', out_handle.name,
           '-userfields', 'query+target+qlo+qhi+tlo+thi+alnlen+evalue+id',
           '-qmask', 'none',
           '-dbmask', 'none',
           '-threads', '1']

    # Write usearch input fasta file
    SeqIO.write(seq, in_handle, format='fasta-2line')
    in_handle.seek(0)

    # Run ublast algorithm
    try:
        stdout_str = check_output(cmd, stderr=STDOUT, shell=False, universal_newlines=True)

        # child = Popen(cmd, bufsize=1, stdout=PIPE, stderr=STDOUT, shell=False,
        #              universal_newlines=True)
        # while child.poll() is None:
        #    out = child.stdout.readline()
        #    sys.stdout.write(out)
        #    sys.stdout.flush()
        # child.wait()
    except CalledProcessError:
        return None

    # Parse usearch output
    field_names = ['query', 'target', 'query_start', 'query_end',
                   'target_start', 'target_end',
                   'length', 'evalue', 'identity']
    align_df = pd.read_csv(out_handle, header=None, names=field_names, encoding='utf-8', sep='\t')
    # Convert to base-zero indices
    align_df[['query_start', 'query_end', 'target_start', 'target_end']] -= 1

    # Close temp file handles
    in_handle.close()
    out_handle.close()

    return align_df


def runBlastn(seq, database, evalue=default_evalue, max_hits=default_max_hits,
              aligner_exec=default_blastn_exec):
    """
    Aligns a sequence against a reference database using BLASTN

    Arguments:
      seq : a list of SeqRecord objects to align.
      database : the path and name of the blastn database.
      evalue : the E-value cut-off.
      maxhits : the maximum number of hits returned.
      aligner_exec : the path to the blastn executable.

    Returns:
      pandas.DataFrame : Alignment results.
    """
    seq_fasta = seq.format('fasta')

    # Define blastn command
    cmd = [aligner_exec,
           '-query',  '-',
           '-db', database,
           '-strand', 'plus',
           '-evalue', str(evalue),
           '-max_target_seqs', str(max_hits),
           '-word_size', '9',
           '-dust', 'no',
           # '-reward', '2',
           # '-penalty', '-1',
           #'-num_descriptions', str(max_hits),
           #'-num_alignments', str(max_hits),
           #'-max_hsps', str(max_hits),
           '-outfmt', '6 qseqid sseqid qstart qend sstart send length evalue pident',
           '-num_threads', '1']

    # Run blastn
    child = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE,
                 shell=False, universal_newlines=True)
    stdout_str, stderr_str = child.communicate(seq_fasta)
    out_handle = StringIO(stdout_str)

    # child = Popen(cmd, bufsize=1, stdout=PIPE, stderr=STDOUT, shell=False,
    #              universal_newlines=True)
    # while child.poll() is None:
    #    out = child.stdout.readline()
    #    sys.stdout.write(out)
    #    sys.stdout.flush()
    # child.wait()

    # Parse blastn output
    field_names = ['query', 'target', 'query_start', 'query_end', 'target_start', 'target_end',
                   'length', 'evalue', 'identity']
    align_df = pd.read_csv(out_handle, header=None, names=field_names, sep='\t')
    # Convert to base-zero indices
    align_df[['query_start', 'query_end', 'target_start', 'target_end']] -= 1

    return align_df
