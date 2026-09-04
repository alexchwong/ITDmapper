__version__ = '1.0.0'


import datetime
import multiprocessing
import argparse

import decimal as dc
dc.getcontext().prec = 5

import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)

import copy
import subprocess
import os
import gzip
import timeit
import shutil

from tqdm import trange
from collections import Counter
from Bio import Align

def save_config(config, filename):
    """
    Write timestamp and commandline arguments to file.

    Args:
        config (dict): Config parameters and values to write.
        filename (str): Name of the file to write to.
    """
    with open(filename, "w") as f:
        f.write("Commandline_argument\tValue\n")
        f.write(f'Time\t{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%d")}\n')
        f.write(f'mergeITD_version\t{__version__}\n')
        for param in sorted(config.keys()):
            if param not in ["ANNO", "DOMAINS"]:
                f.write(f'{param}\t{config[param]}\n')

def load_config(filename):
    """
    Load config parameters from file.

    Args:
        filename (str): Name of the file to read config from.

    Returns:
        Dictionary with config parameter - value pairs.
    """
    config = {}
    with open(filename, "r") as f:
        for line in f:
            key, val = line.strip("\n").split("\t")
            if key not in ["Time", "Commandline_argument"]:
                try:
                    config[key] = int(val)
                except:
                    try:
                        config[key] = float(val)
                    except:
                        config[key] = val

    # recognize string as dict
    # if "COST_ALIGNED" in config:
        # config["COST_ALIGNED"] = eval(config["COST_ALIGNED"])
    return config



# child processes spawned on Windows by multiprocessing do not
# receive variables set in __main__ of parent process
# -->  they cannot access config {} values set in __main__
# --> to circumvent this, __main__ saves config and children
#     spawned by multiprocessing load it from file
if __name__ in ['__mp_main__', 'getitd']:
    try:
        current_dir = os.getcwd()
        config = load_config(os.path.join(current_dir, "config.txt"))
    except OSError:
        print("NO CONFIG FOUND")
        config = {}
else:
    # mimic Windows style process spawning on Linux:
    # multiprocessing.set_start_method("spawn")
    try:
        config
    except NameError:
        config = {}


def parallelize(function, args, cores):
    """
    Parallelize a given function across a given number of cores.

    Args:
        function (function): Function or method to parallelize.

        args (tuple): Tuple of function's arguments.

        cores (int): Number of cores to utilize.

    Returns:
        List of function's outputs.
    """
    with multiprocessing.Pool(cores) as p:
        return p.map(function, args)

def bbmap_process(config):
    fastq1, fastq2 = config["R1"], config["R2"]
    bbmap_path, temp_path = config["BBMAP_PATH"], config["TMP_DIR"]

    assert os.path.isdir(bbmap_path)
    assert os.path.isfile(os.path.join(bbmap_path, 'bbmerge.sh'))
    assert os.path.isfile(os.path.join(bbmap_path,'bbduk.sh'))
    
    start_time = timeit.default_timer()
    
    if os.path.isdir(temp_path):
        shutil.rmtree(temp_path)
        
    bbmap_log = ""
    
    # Phase 1 merging
    ret = subprocess.run([
        os.path.join(bbmap_path, 'bbmerge.sh'), 
        f'in1={fastq1}', f'in2={fastq2}',
        f'out={os.path.join(temp_path, "merged.fastq")}', 
        f'outu={os.path.join(temp_path, "unmerged.fastq")}',
        f'ihist={os.path.join(config["OUT_DIR"], "hist.tsv")}'
    ], capture_output=True, text=True)    
    
    blog = ret.stderr
    logLines = blog.split('\n')
    numPairs = 0
    numMerged = 0
    for i in range(len(logLines)):
        if logLines[i].find("Pairs:") > -1:
            numPairs = int(logLines[i].split('\t')[1])
        elif logLines[i].find("Joined:") > -1:
            numMerged = int(logLines[i].split('\t')[1])
    
    assert numPairs > 0
    pcMerged = round(numMerged * 100 / numPairs, 2)
    tt = round(timeit.default_timer() - start_time, 2)
    if config["BBMAP_TRIMQ"] > -1:
        save_stats(f'Phase 1 merging {numMerged} of {numPairs} pairs merged ({pcMerged} %) - {tt} sec', config["STATS_FILE"])
    else:
        save_stats(f'Merging {numMerged} of {numPairs} pairs merged ({pcMerged} %) - {tt} sec', config["STATS_FILE"])
        
    bbmap_log += blog + '\n'
    
    # Phase 2 merging
    numPairs2 = 0
    numMerged2 = 0
    if config["BBMAP_TRIMQ"] > -1:
        start_time2 = timeit.default_timer()
        ret = subprocess.run([
            f'{os.path.join(bbmap_path, "bbduk.sh")}', 
            f'in={os.path.join(temp_path, "unmerged.fastq")}', 
            f'out={os.path.join(temp_path, "qtrimmed.fastq")}',
            "qtrim=r", f'trimq={config["BBMAP_TRIMQ"]}'
        ], capture_output=True, text=True)
        bbmap_log += ret.stderr + '\n'
        
        ret = subprocess.run([
            f'{os.path.join(bbmap_path, "bbmerge.sh")}', 
            f'in={os.path.join(temp_path, "qtrimmed.fastq")}', 
            f'out={os.path.join(temp_path, "merged2.fastq")}',
            f'ihist={os.path.join(config["OUT_DIR"], "hist2.tsv")}'
        ], capture_output=True, text=True)
        
        # Concatenate into first file
        f1 = open(os.path.join(temp_path, "merged.fastq"), 'a+')
        f2 = open(os.path.join(temp_path, "merged2.fastq"), 'r')
        f1.write(f2.read())
        f1.close()
        f2.close()

        blog = ret.stderr
        logLines = blog.split('\n')
        for i in range(len(logLines)):
            if logLines[i].find("Pairs:") > -1:
                numPairs2 = int(logLines[i].split('\t')[1])
            elif logLines[i].find("Joined:") > -1:
                numMerged2 = int(logLines[i].split('\t')[1])
        
        if numPairs2 > 0:
            tt = round(timeit.default_timer() - start_time2, 2)
            pcMerged2 = round(numMerged2 * 100 / numPairs2, 2)
            pcMergedT2 = round(numMerged2 * 100 / numPairs, 2)
            save_stats(f'Phase 2 merging {numMerged2} of {numPairs2} pairs merged ({pcMerged2} % of unmerged, {pcMergedT2} % of total) - {tt} sec', config["STATS_FILE"])
        else:
            save_stats("Phase 2 merging not performed due to zero reads available", config["STATS_FILE"])
            
        bbmap_log += blog + '\n'


    # Phase 3 - average bqs filtering
    if config["BBMAP_BQS"] > -1:
        start_time3 = timeit.default_timer()
        ret = subprocess.run([
            f'{os.path.join(bbmap_path, "bbduk.sh")}', 
            f'in={os.path.join(temp_path, "merged.fastq")}', 
            f'out={os.path.join(temp_path, "cleaned.fastq")}',
            f'maq={config["BBMAP_BQS"]}'
        ], capture_output=True, text=True)

        blog = ret.stderr
        logLines = blog.split('\n')
        numInput = 0
        numRetained = 0
        for i in range(len(logLines)):
            if logLines[i].find("Input:") > -1:
                txtInput = logLines[i].split('\t')[1]
                txtInput = txtInput.split(' ')[0]
                numInput = int(txtInput)
            elif logLines[i].find("Result:") > -1:
                txtRetained = logLines[i].split('\t')[1]
                txtRetained = txtRetained.split(' ')[0]
                numRetained = int(txtRetained)
        
        assert numInput > 0
        tt = round(timeit.default_timer() - start_time3, 2)

        pcRetained = round(numRetained * 100 / numInput, 2)
        save_stats(f'BQS filtering retained {numRetained} of {numInput} merged reads ({pcRetained} %) - {tt} sec', config["STATS_FILE"])

        pcTotal = round(numRetained * 100 / numPairs, 2)
        save_stats(f'BBmap merging final result: {numRetained} of {numPairs} merged reads ({pcTotal} %)', config["STATS_FILE"])
            
        bbmap_log += blog + '\n'
    else:
        os.rename(
            os.path.join(temp_path, "merged.fastq"), 
            os.path.join(temp_path, "cleaned.fastq")
        )
        if numMerged2 > 0:
            numMerged3 = numMerged + numMerged2
            pcTotal = round(numMerged3 * 100 / numPairs, 2)
            save_stats(f'BBmap merging final result: {numMerged3} of {numPairs} merged reads ({pcTotal} %)', config["STATS_FILE"])
    
    save_stats(f'BBmap total time taken - {round(timeit.default_timer() - start_time, 2)} sec',
        config["STATS_FILE"])

    with open(config["BBLOG"], 'w') as f:
        f.write(bbmap_log)

    assert os.path.isfile(os.path.join(temp_path, "cleaned.fastq"))
    return os.path.join(temp_path, "cleaned.fastq")

def is_gz_file(filename):
    """
    Check whether a given file is gzipped or not,
    using its magic number.

    Args:
        filename: Name of the file to read.

    Returns:
        bool, True when gzipped, False otherwise.
    """
    with open(filename, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'

def read_fastq(fastq_file):
    """
    Read sequence fastq file and extract sequences and BQS.

    Args:
        fastq_file: Name of the fastq file to read, R1 or R2.

    Returns:
        List of Read() objects.
    """
    reads = [] # simple list of fastq sequences
    readbqs = []
    try:
        if is_gz_file(fastq_file):
            open_fct = gzip.open
        else:
            open_fct = open
        
        with open_fct(fastq_file, 'rt') as f:
            line = f.readline()
            while line:
                _ = line
                read_seq = f.readline().rstrip(os.linesep)
                _ = f.readline()
                _ = f.readline().rstrip(os.linesep)
                reads.append(read_seq)
                # readbqs.append(average_bqs(read_bqs))
                
                line = f.readline()
    except IOError as e:
        print(f'---\nCould not read fastq file {fastq_file}!\n---')
    
    return reads, readbqs

def average_bqs(bqs):
    """
    Calculate the mean BQS of a given string of quality scores.
    Assumes BQS are in Sanger format, encoded as Phred +33.

    Args:
        bqs (str): String of base quality scores.

    Returns:
        Mean BQS of that string.
    """
    return sum([ord(x) - 33 for x in bqs]) / len(bqs)


def read_reference(filename):
    """
    Read in WT reference sequence.

    Args:
        filename (str): Name of the file to be read.

    Returns:
        Reference sequence, stripped of trailing newlines.
    """
    with open(filename, 'r') as f:
        ref = f.read()
    ref = ref.splitlines()
    assert len(ref) == 1
    return ref[0]

# add column names!
def read_annotation(filename):
    """
    Read in WT reference annotation file.

    For each bp of the WT reference, provides genomic, transcriptomic
    and proteomic coordinate, exon/intron annotation and the respective
    reference bp.

    Args:
        filename (str): Name of the file to be read.

    Returns:
        pd.DataFrame of the annotation.
    """
    try:
        return pd.read_csv(filename, sep='\t')
    except IOError as e:
        print("\nAnnotation file was not provided or cannot be accessed!\n")
        return None

def annotateCoords(anno_df):
    """
    Annotate HGVS coordinates to given reference file
    """
    df = anno_df.copy(deep=True)
    
    # find first exon
    firstExonCoord = 0
    for i in range(len(df)):
        if df.iloc[i]["region"].find("exon") > -1:
            firstExonCoord = i
            break

    assert firstExonCoord > -1
    assert int(df.iloc[firstExonCoord]["transcript_bp"]) > 0
    df["HGVScoord"] = ""
    
    # Annotate upstream intron if required
    if firstExonCoord > 0:
        cdot = int(df.iloc[i]["transcript_bp"])
        for i in range(firstExonCoord, -1, -1):
            df.loc[i, "HGVScoord"] = f'{int(cdot)}{int(i)-int(firstExonCoord)}'
    
    # Now annotate first exon
    i = firstExonCoord
    inIntron = False
    while i < len(anno_df):
        if df.iloc[i]["region"].find("exon") > -1:
            if inIntron:
                inIntron = False
            df.loc[i, "HGVScoord"] = f'{int(df.iloc[i]["transcript_bp"])}'
        elif not inIntron:
            inIntron = True
            cdot = int(df.iloc[i-1]["transcript_bp"])
            # find next exon coord
            nextExonCoord = 0
            lastExonCoord = i-1
            for j in range(i, len(df)):
                if df.iloc[j]["region"].find("exon") > -1:
                    nextExonCoord = j
                    break
            # annotate first base of intron
            df.loc[i, "HGVScoord"] = f'{int(cdot)}+{1}'
        elif nextExonCoord > 0:
            # need to find whether position is closer to donor or acceptor
            if i - lastExonCoord < nextExonCoord - i:
                df.loc[i, "HGVScoord"] = f'{int(cdot)}+{i - lastExonCoord}'
            else:
                df.loc[i, "HGVScoord"] = f'{int(cdot)+1}-{nextExonCoord - i}'
        else:
            assert inIntron
            df.loc[i, "HGVScoord"] = f'{int(cdot)}+{i - lastExonCoord}'
        i += 1

    return df

def ar_to_vaf(ar):
    """
    Convert AR to VAF.

    VAF (variant allele frequency) = V-AF
    AR (allele ratio) = V-AF / WT-AF
    V-AF + WT-AF = 100 (%)

    Args:
        ar (float): AR to convert.

    Returns:
        VAF (float)

    """
    return ar/(ar + 1) * 100 # * 100 because VAF is in %

def vaf_to_ar(vaf):
    """
    Convert VAF to AR.

    VAF (variant allele frequency) = V-AF
    AR (allele ratio) = V-AF / WT-AF
    V-AF + WT-AF = 100 (%)

    Note:
        if VAF == 100:
            AR = -1
            (instead of 100 / 0)

    Args:
        vaf (dc.Decimal): VAF to convert.

    Returns:
        AR (dc.Decimal)
    """
    if vaf == 100:
        return -1
    return vaf/(100 - vaf)

def save_stats(stat, filename):
    """
    Write statistics to file.

    Args:
        stat (str): Statistic to save.
        filename (str): Name of the file to write to.
    """
    print(stat)
    with open(filename, "a") as f:
        f.write(stat + "\n")

def str_to_bool(string):
    """
    Convert literal str to bool, such that
        'True' ->  True
        'False' -> False
    This is required to parse boolean command line
    arguments, which are all initially of type str.
    Calling `bool(str)` evaluates to `True` whenever
    the str is not empty. Thus, without this function,
    'False' would also evaluate to `True`.

    Args:
        string (str): To be converted.

    Returns:
        Literal boolean (bool) of string.
    """
    if string.lower() in ('true', 'True'):
        return True
    elif string.lower() in ('false', 'False'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value (True or False) expected.')

class Mutation(object):
    def __init__(
        self,
        mutType = "ins",
        pos_str = None,
        ins_str = "",
        counts = 0):
        
        self.mutType = mutType
        assert pos_str is not None
        self.pos = [int(i) for i in pos_str.split("_")]
        if len(self.pos) == 1:
            self.pos.append(self.pos[0])
        self.ins_str = ins_str
        self.counts = counts
        self.comutations = {}
        
    def add(self, count):
        self.counts += count
        return(self)
    
    def nameMut(self):
        if self.mutType == "snp":
            return(f'{self.pos[0]}{self.ins_str}')
        elif self.pos[0] == self.pos[1]:
            return(f'{self.pos[0]}{self.mutType}{self.ins_str}')
        return(f'{self.pos[0]}_{self.pos[1]}{self.mutType}{self.ins_str}')

    def nameActualMut(self, config):
        pos0 = config["ANNO"].iloc[self.pos[0]]["HGVScoord"]
        if self.mutType == "snp":
            return(f'c.{pos0}{self.ins_str}')
        elif self.pos[0] == self.pos[1]:
            return(f'c.{pos0}{self.mutType}{self.ins_str}')
        pos1 = config["ANNO"].iloc[self.pos[1]]["HGVScoord"]
        return(f'c.{pos0}_{pos1}{self.mutType}{self.ins_str}')
    
    def addComutation(self, mutName, counts):
        if mutName in self.comutations.keys():
            self.comutations[mutName] += counts
        else:
            self.comutations[mutName] = counts
        return(self)
    
    def netInsert(self):
        if self.mutType[0:3] == "ins":
            return len(self.ins_str)
        elif self.mutType == "del":
            return -self.pos[1] + self.pos[0] - 1
        elif self.mutType[0:6] == "delins":
            return len(self.ins_str) - self.pos[1] + self.pos[0] - 1
        elif self.mutType == "dup":
            return self.pos[1] - self.pos[0] + 1
        else:
            return 0
        
    def getInsertPos(self):
        if self.mutType[0:3] == "ins":
            return self.pos[0]
        elif self.mutType == "del":
            return self.pos[0]
        elif self.mutType[0:6] == "delins":
            return self.pos[0]
        elif self.mutType == "dup":
            return self.pos[1] + 1
        else:
            return self.pos[0]

def getHGVS(seq, ref, config, verbose = False):
    """
    Aligns a sequence with respect to the reference
    - returns 3 lists: coordinates, operations, insert sequence
    - Also returns start / end coordinates of alignment (-1 if not aligned)
    
    Hardcoded parameters (for now)
    - minAlignLen (6): number of consecutive nucleotides required for a block alignment
    - minRefAlignFraction (0.4): min fraction of reference aligned by read
    - maxFracIsIndel (0.7): max fraction of merged read that is part of in/del (i.e. not aligned to reference)
    """
    
    q_alns = []
    r_alns = []

    aligner2 = config["ALIGNER"]

    aln = aligner2.align(seq, ref)
    if not aln:
        return [],[],[], -1, -1
    if verbose:
        print(aln[-1])
        
    coords = aln[-1].coordinates
    # print(coords)
    minAlignLen = config["MIN_ALIGN_LEN"]
    for j in range(len(coords[0]) - 1):
        if coords[0][j+1] - coords[0][j] >= minAlignLen and coords[1][j+1] - coords[1][j] >= minAlignLen:
            # if match must be a block of at least minAlignLen
            q_alns.append(coords[0][j:j+2])
            r_alns.append(coords[1][j:j+2])
            
    q_alns = np.array(q_alns)
    r_alns = np.array(r_alns)

    # No alignment
    if q_alns.shape[0] == 0:
        return [],[],[], -1, -1

    # Boundaries of alignment
    ref_start, ref_end = r_alns[0][0], r_alns[-1][1]

    # Insufficient alignment to reference
    minRefAlignFraction = config["MIN_REF_ALN_FRACTION"]
    if (ref_end - ref_start) < minRefAlignFraction * len(ref):
        return [],[],[], -1, -1                
    
    # Single alignment - likely WT
    if q_alns.shape[0] == 1:
        return [],[],[], ref_start, ref_end
   
    if verbose:
        print(q_alns)
        print(r_alns)

    # calculate the sum inserted / deleted. Abort if this exceeds 70% of length of sequence
    nTotIndel = 0    
    ops, rC, rSeq = [],[],[]
    for i in range(len(q_alns) - 1):
        if q_alns[i+1][0] == q_alns[i][1]:
            # no gaps in query sequence, i.e. no insertion
            if r_alns[i+1][0] > r_alns[i][1]:
                # gap in reference alignment deletion
                ops.append("del")
                rSeq.append("")
                rC.append(f'{str(r_alns[i][1])}_{str(r_alns[i+1][0]-1)}')
                nTotIndel += (r_alns[i][1] - r_alns[i+1][0])
            else:
                pass
        else:
            if r_alns[i+1][0] == r_alns[i][1]:
                # simple insertion
                insSeq = seq[q_alns[i][1]:q_alns[i+1][0]]
                nTotIndel += len(insSeq)
                
                # check if duplication
                isDup = False
                if len(insSeq) <= q_alns[i][1]:
                    dupSeq = seq[(q_alns[i][1] - len(insSeq)):(q_alns[i+1][0] - len(insSeq))]
                    if insSeq == dupSeq:
                        isDup = True
                        ops.append("dup")
                        rSeq.append("")
                        rC.append(f'{str(r_alns[i][1] - len(insSeq))}_{str(r_alns[i][1] - 1)}')
                if not isDup:
                    rC.append(f'{str(r_alns[i][1])}_{str(r_alns[i+1][0]+1)}')
                    ops.append(f'ins[{len(insSeq)}]')
                    rSeq.append(insSeq)

            elif r_alns[i+1][0] > r_alns[i][1]:
                # delins
                insSeq = seq[q_alns[i][1]:q_alns[i+1][0]]
                rC.append(f'{str(r_alns[i][1])}_{str(r_alns[i+1][0]-1)}')
                ops.append(f'delins[{r_alns[i+1][0] - r_alns[i][1]},{len(insSeq)}]')
                rSeq.append(insSeq)
                nTotIndel += len(insSeq) + (r_alns[i+1][0] - r_alns[i][1]) # del + ins
            else:
                # overlapping alignment preceeded by novel insert, treat as insertion
                # true insert length is longer than mapped
                # attach duplicated alignment to end of novel insert
                dupLen = r_alns[i][1] - r_alns[i+1][0]
                insSeq = seq[q_alns[i][1]:(q_alns[i+1][0] + dupLen)]
                rC.append(f'{str(r_alns[i][1])}_{str(r_alns[i][1]+1)}')
                ops.append(f'ins[{len(insSeq)}]')
                rSeq.append(insSeq)

    maxFracIsIndel = config["MAX_FRAC_INDEL"]
    if (nTotIndel + ref_end - ref_start) * maxFracIsIndel < nTotIndel:
        return [],[],[], -1, -1

    if verbose:
        print([f'{c}{o}{s}' for c, o, s in zip(rC, ops, rSeq)])
    
    return rC, ops, rSeq, ref_start, ref_end

def generateMutSeq(ref, hgvs, returnComplex = False):
    """
    Generates a sequence which is a mutated sequence from the given reference,
    using the given hgvs instructions (inserts / deletions only, not SNPs)
    """
    
    refStarts = []
    refEnds = []
    insSeqs = []
    
    for hg in hgvs:
        posStart = -1
        posEnd = -1
        seq = ""
        if hg.find("delins") > -1:
            mutStart = hg.find("delins")
            pos = hg[0:mutStart]
            seq = hg[mutStart+6:]
            if seq.find("]") > -1:
                seq = seq.split("]")[1]
            if pos.find("_") > 0:
                posStart, posEnd = pos.split("_")
                posStart = int(posStart)
                posEnd = int(posEnd) + 1
            else:
                posStart = int(pos)
                posEnd = int(pos) + 1
        elif hg.find("ins") > -1:
            mutStart = hg.find("ins")
            pos = hg[0:mutStart]
            seq = hg[mutStart+3:]
            if seq.find("]") > -1:
                seq = seq.split("]")[1]
            if pos.find("_") > 0:
                posStart, posEnd = pos.split("_")
                posStart = int(posStart)
                posEnd = int(posEnd) - 1
            else:
                posStart = int(pos)
                posEnd = int(pos) - 1
        elif hg.find("del") > -1:
            mutStart = hg.find("del")
            pos = hg[0:mutStart]
            if pos.find("_") > 0:
                posStart, posEnd = pos.split("_")
                posStart = int(posStart)
                posEnd = int(posEnd) + 1
            else:
                posStart = int(pos)
                posEnd = int(pos) + 1
        elif hg.find("dup") > -1:
            mutStart = hg.find("dup")
            pos = hg[0:mutStart]
            if pos.find("_") > 0:
                posStart, posEnd = pos.split("_")
                posStart = int(posStart)
                posEnd = int(posEnd) + 1
            else:
                posStart = int(pos)
                posEnd = int(pos) + 1
            seq = ref[posStart:posEnd]
            posEnd = posStart # as this is duplication
        else:
            pass
        
        if posStart > -1:
            refStarts.append(posStart)
            refEnds.append(posEnd)
            insSeqs.append(seq)
        
    refStarts = np.array(refStarts)
    refEnds = np.array(refEnds)
    insSeqs = np.array(insSeqs)
    
    arg_sort = np.argsort(refStarts)
    refStarts = refStarts[arg_sort]
    refEnds = refEnds[arg_sort]
    insSeqs = insSeqs[arg_sort]
    
    refPos = 0
    newRef = ""
    for rS, rE, iS in zip(refStarts, refEnds, insSeqs):
        if rS > refPos:
            newRef += ref[refPos:rS]
            refPos = rE
            newRef += iS
    
    newRef += ref[refPos:]
    
    if returnComplex:
        return newRef, refStarts, refEnds, insSeqs
    
    return newRef

def getRefLoc(pos, refStarts, refEnds, insSeqs):
    # gets the true reference position given the outputs of generateMutSeq
    
    fudgeF = 0
    for rS, rE, iS in zip(refStarts, refEnds, insSeqs):
        if pos + fudgeF < rS:
            return pos + fudgeF
        altF = rE - rS + len(iS)
        if altF > 0:
            # net insertion
            if pos + fudgeF - rS < altF:
                # inside insert
                return -1
        fudgeF -= altF
    return pos + fudgeF

def findSNP(seq, ref, delins_hgvs, config):
    """
    Finds and describes any SNPs of sequence, with respect to reference
    mutated with deletion / insertion hgvs instructions
    """
    
    synRef, refStarts, refEnds, insSeqs = generateMutSeq(ref, delins_hgvs, returnComplex = True)
    
    aligner2 = config["ALIGNER"]
    aln = aligner2.align(seq, synRef)
    
    q_alns, r_alns = [], []
    coords = aln[-1].coordinates
    minAlignLen = 6
    for j in range(len(coords[0]) - 1):
        if coords[0][j+1] - coords[0][j] >= minAlignLen and coords[1][j+1] - coords[1][j] >= minAlignLen:
            # if match must be a block of at least 6
            q_alns.append(coords[0][j:j+2])
            r_alns.append(coords[1][j:j+2])
            
    q_alns = np.array(q_alns)
    r_alns = np.array(r_alns)

    snpCoords = []
    snpRefs = []
    snpSubs = []
    
    for i in range(len(q_alns)):
        qlen = q_alns[i][1] - q_alns[i][0]
        rlen = r_alns[i][1] - r_alns[i][0]

        qseq = seq[q_alns[i][0]:q_alns[i][1]]
        rseq = synRef[r_alns[i][0]:r_alns[i][1]]
        assert len(qseq) == len(rseq)
            
        r_start = r_alns[i][0]
        for j in range(qlen):
            if qseq[j] != rseq[j]:
                snpCoords.append(r_start + j)
                snpRefs.append(rseq[j])
                snpSubs.append(qseq[j])
    
    # Correct snpCoords
    realCoords, realRefs, realSubs = [], [], []
    for i in range(len(snpCoords)):
        pos = getRefLoc(snpCoords[i], refStarts, refEnds, insSeqs)
        if pos > 0:
            realCoords.append(pos)
            realRefs.append(snpRefs[i])
            realSubs.append(snpSubs[i])

    return realCoords, realRefs, realSubs

def get_amplicons(config):
    ampliconome = config["OME_FILE"]
    amp_names, amp_lens = [], []
    header, length = None, 0
    with open(ampliconome) as fasta:
        for line in fasta:
            # Trim newline
            line = line.rstrip() # remove return carriage and any trailing spaces
            if line.startswith('>'):
                # If we captured one before, print it now
                if header is not None:
                    amp_names.append(header)
                    amp_lens.append(length)
                    length = 0
                header = line[1:]
            else:
                line.replace(" ", "") # remove spaces
                length += len(line)
    # final seq
    amp_names.append(header)
    amp_lens.append(length)        
    return amp_names, amp_lens

def generate_bam(config, mode = "PairedEnd"):
    outPrefix = mode
    
    sampleName = config["SAMPLE"]
    samplePath = config["OUT_DIR"]
    initSam = os.path.join(samplePath, f"{outPrefix}.sam")
    tmpSam = os.path.join(samplePath, "tmp.sam")
    cleanedBam = os.path.join(samplePath, f"{outPrefix}_cleaned.bam")
    sortedCleanedBam = os.path.join(samplePath, f"{outPrefix}_cleaned_sorted.bam")
    bbmap_path = config["BBMAP_PATH"]

    assert os.path.isdir(bbmap_path)
    assert os.path.isfile(os.path.join(bbmap_path, "bbmap.sh"))

    start_time = timeit.default_timer()
    save_stats(f'\nUsing BBMap to align {outPrefix} reads', config["STATS_FILE"])
    
    # Get sequences names and lengths from ampliconome
    amp_names, amp_lens = get_amplicons(config)
    assert len(amp_names) > 1
    
    if mode == "PairedEnd":
        ret = subprocess.run([
            os.path.join(bbmap_path, "bbmap.sh"), 
            f'in={config["R1"]}', f'in2={config["R2"]}',
            f'ref={config["OME_FILE"]}',
            f'out={initSam}',
            "maxindel=2", "strictmaxindel=t",
            f'minaveragequality={config["BBMAP_BQS"]}',
            "nodisk"
        ], capture_output=True, text=True)
    elif mode == "Merged":
        ret = subprocess.run([
            os.path.join(bbmap_path, "bbmap.sh"), 
            f'in={config["MERGED_READS"]}',
            f'ref={config["OME_FILE"]}',
            f'out={initSam}',
            "maxindel=2", "strictmaxindel=t",
            f'minaveragequality={config["BBMAP_BQS"]}',
            "nodisk"
        ], capture_output=True, text=True)

    tt = round(timeit.default_timer() - start_time, 2)
    save_stats(f'BBMap {outPrefix} initial alignment completed - {tt} sec', config["STATS_FILE"])
    
    start_time2 = timeit.default_timer()
    # filter reads by fragment length > amplicon length minus max unaligned
    alignDiff = config["ALIGN_LEN_DIFF"]
    with open(tmpSam, 'w') as tmp:
        with open(tmpSam, 'a') as tmp:
            ps2 = subprocess.Popen(["awk", 'substr($0,1,1)=="@"', initSam], stdout = tmp)
            p_status = ps2.wait()
            
    printAll = "{print $0}"
    for i in range(len(amp_names)):
        fLen = amp_lens[i] - alignDiff
        fLen2 = amp_lens[i] + alignDiff
        mapqT = config["ALIGN_MAPQ"]
        with open(tmpSam, 'a') as tmp:
            if mode == "PairedEnd":
                ps2 = subprocess.Popen(
                    ["awk", "-F\t", f'($3 == "{amp_names[i]}") && ($5 >= {mapqT}) && (($9 >= {fLen}) || ($9 <= -{fLen-1})) {printAll}', initSam], 
                    stdout = tmp)
                p_status = ps2.wait()
            elif mode == "Merged":
                ps2 = subprocess.Popen(
                    ["awk", "-F\t", f'($3 == "{amp_names[i]}") && ($5 >= 30) && (length($10) >= {fLen}) && (length($10) <= {fLen2}) {printAll}', initSam], 
                    stdout = tmp)
                p_status = ps2.wait()

    os.remove(initSam)
    with open(cleanedBam, 'w') as bam:
        pb = subprocess.Popen(["samtools", "view", "-b", tmpSam], stdout = bam)
        p_status = pb.wait()

    os.remove(tmpSam)    
    with open(sortedCleanedBam, 'w') as bam:
        pb = subprocess.Popen(["samtools", "sort", cleanedBam], stdout = bam)
        p_status = pb.wait()

    os.remove(cleanedBam)
    ret = subprocess.run(["samtools", "index", sortedCleanedBam])

    tt = round(timeit.default_timer() - start_time2, 2)
    save_stats(f'{outPrefix} Alignment fidelity filter complete - {tt} sec', config["STATS_FILE"])
    
    idxStatsFile = os.path.join(config["OUT_DIR"], f'{outPrefix}_idxstats.txt')
    with open(idxStatsFile, 'w') as log:
        p = subprocess.Popen(["samtools", "idxstats", sortedCleanedBam], stdout=log)
        p_status = p.wait()
    
    # clean final idxstats
    idx = pd.read_csv(idxStatsFile, sep = '\t', header = None)
    amplicon_names= idx.iloc[:-1, 0].tolist()
    amplicon_lens = idx.iloc[:-1, 1].tolist()    
    amplicon_aligns = idx.iloc[:-1, 2].tolist()    
    sum_aligned = sum([int(i) for i in amplicon_aligns])
    vafs = [i * 100 / sum_aligned for i in amplicon_aligns]
    
    res_s = pd.read_csv(config["MUTATION_FILE_FILTERED"])
    mut_names = res_s["name"].tolist()
    mut_names.insert(0, "WildType")
    
    dict = {'Amplicon': mut_names, 'Alias': amplicon_names,
        'Length': amplicon_lens, 'Reads Aligned': amplicon_aligns,
        'VAF%': vafs}
    res = pd.DataFrame(dict)
    res.to_csv(os.path.join(samplePath, f'{outPrefix}_aligned_stats.csv'),
        sep = ",", index=False)
    os.remove(idxStatsFile)
    return(0)
    

def alignITD(prealigns_df, config):
    REF = config["REF"]
    aligner2 = config["ALIGNER"]
    
    start_time = timeit.default_timer()
    
    df = prealigns_df.copy(deep=True)

    allow_insert_mismatch = 3

    df["Aligned"] = False
    df["alignRefCoords"] = ""
    df["HGVS"] = ""
    df["SNP"] = ""
    df["absHGVS"] = ""
    df["absSNP"] = ""
    df["net_insertSize"] = 0
    df["Ns"] = ""
    # df["idealSequence"] = ""
    # df["gaps_and_mismatches"] = -1
    # df["Ns"] = -1
    # df["nonN_gaps_and_mismatches"] = -1

    covIncrement = np.zeros((len(REF) + 1,))
    mutList = []
    mutNames = []

    if config["PROGRESSBAR"]:
        opt_range = trange
    else:
        opt_range = range

    for i in opt_range(len(df)):
        seq = df.iloc[i]["Sequence"]
        rC_S, ops_S, rSeq, startC, endC = getHGVS(seq, REF, config)
        if startC == -1 or endC == -1:
            continue

        seqCount = int(df.iloc[i]["Counts"])
            
        covIncrement[startC] += seqCount
        covIncrement[endC] -= seqCount
        
        seqMutList = []
        for rC, ops, rS in zip(rC_S, ops_S, rSeq):
            tmpMut = Mutation(ops, rC, rS, seqCount)
            seqMutList.append(tmpMut)
        
        if any(m.netInsert() > 2 for m in seqMutList):
            for ii in range(len(seqMutList)):
                netInsert = seqMutList[ii].netInsert()
                if netInsert > 2:
                    thisMutName = seqMutList[ii].nameMut()
                    # find highest count mutation with match under threshold
                    mutOpts = [m for m in mutList if (m.netInsert() >= netInsert - 3 and m.netInsert() <= netInsert + 3)]

                    if len(mutOpts) > 0:                                              
                        mutOpts.sort(key=lambda x: x.counts, reverse=True)

                        hasScore = False
                        for mut in mutOpts:
                            if mut.nameMut() == thisMutName:
                                seqMutList[ii] = mut
                                break
                                
                            if not hasScore:
                                # get ideal scoring
                                mutNameList = [m.nameMut() for m in seqMutList]
                                synSeq = generateMutSeq(REF, mutNameList)
                                aln = aligner2.align(seq, synSeq)
                                cts_ideal = aln[-1].counts()
                                hasScore = True
                                
                            tmpHGVS = copy.deepcopy(mutNameList)
                            tmpHGVS[ii] = mut.nameMut()
                            synSeq = generateMutSeq(REF, tmpHGVS)
                            aln = aligner2.align(seq, synSeq)
                            cts = aln[-1].counts()
                            if cts.gaps + cts.mismatches <= cts_ideal.gaps + cts_ideal.mismatches + allow_insert_mismatch:
                                seqMutList[ii] = mut
                                break

        HGVSMutNames = [m.nameMut() for m in seqMutList]
        actualHGVSMutNames = [m.nameActualMut(config) for m in seqMutList]
                                                
        # Infer SNPs
        snpNameList = []
        actualSnpNameList = []
        N_List = []

        snpCoords, snpRefs, snpSubs = findSNP(seq, REF, HGVSMutNames, config)

        for sC, sR, sS in zip(snpCoords, snpRefs, snpSubs):
            if sS == "N":
                N_List.append(sC)
            else:
                snp = f'{sC}{sR}>{sS}'
                snpNameList.append(snp)
                snpMut = Mutation("snp", str(sC), f'{sR}>{sS}', seqCount)
                seqMutList.append(snpMut)
        
        actualSnpNameList = [m.nameActualMut(config) for m in seqMutList]
        
        # add mutations to main list
        for ii in range(len(seqMutList)):
            thisMutName = seqMutList[ii].nameMut()
            mutIdx = [i for i, m in enumerate(mutNames) if m == thisMutName]
            assert len(mutIdx) < 2
            if len(mutIdx) == 1:
                mutList[mutIdx[0]].add(seqCount)
            else:
                mutList.append(seqMutList[ii])
                mutNames.append(seqMutList[ii].nameMut())                    
                    
        # determine comutations
        if len(seqMutList) > 1:
            seqMutNames = [m.nameMut() for m in seqMutList]
            for ii in range(len(seqMutNames)):
                mutIdx = [i for i, m in enumerate(mutNames) if m == seqMutNames[ii]]
                assert len(mutIdx) == 1
                for iii in range(len(seqMutNames)):
                    if iii == ii:
                        continue
                    mutList[mutIdx[0]].addComutation(seqMutNames[iii], seqCount)                    
                    
        df.loc[i,"Aligned"] = True
        df.loc[i, "alignRefCoords"] = f'{startC}-{endC}'
        mutNameList = [m.nameMut() for m in seqMutList]
        df.loc[i, "HGVS"] = ";".join(HGVSMutNames)
        df.loc[i, "absHGVS"] = ";".join(actualHGVSMutNames)
        df.loc[i, "SNP"] = ";".join(snpNameList)
        df.loc[i, "absSNP"] = ";".join(actualSnpNameList)

        # newSeq = generateMutSeq(REF, mutNameList)
        # df.loc[i, "net_insertSize"] = len(newSeq) - len(REF)
        df.loc[i, "net_insertSize"] = sum([int(m.netInsert()) for m in seqMutList])
        
        df.loc[i, "Ns"] = ",".join([str(i) for i in N_List])
        
        # df.loc[i, "idealSequence"] = newSeq
        # aln2 = aligner2.align(df.iloc[i]["Sequence"], newSeq)
        # cts = aln2[-1].counts()
        # df.loc[i, "gaps_and_mismatches"] = cts.gaps + cts.mismatches
        # df.loc[i, "Ns"] = df.iloc[i]["Sequence"].count('N')
        # df.loc[i, "nonN_gaps_and_mismatches"] = cts.gaps + cts.mismatches - df.iloc[i]["Sequence"].count('N')

    ######################################################################
    # Write results
    df.to_csv(config["ALIGN_FILE"], sep = ",", index=False)
    
    anno = config["ANNO"]
    
    mutList.sort(key=lambda x: x.counts, reverse=True)

    mutName = [m.nameMut() for m in mutList]
    absmutName = [m.nameActualMut(config) for m in mutList]
    mutCount = [m.counts for m in mutList]
    netIns = [m.netInsert() for m in mutList]

    totCov = []
    incCov = 0
    iref_coverage_total, iref_coverage_frwd, iref_coverage_rev = {}, {}, {}
    for i in range(len(covIncrement)):
        incCov += covIncrement[i]
        totCov.append(incCov)
        iref_coverage_total[i] = incCov
        iref_coverage_frwd[i] = 0
        iref_coverage_rev[i] = 0
    iref_coverage = {"all_reads": iref_coverage_total, 
        "forward_reads": iref_coverage_frwd, 
        "reverse_reads": iref_coverage_rev}
    save_coverage(iref_coverage, config)
    if config["PLOT"]:
        plot_coverage(iref_coverage, config)
    
    insertPos = []
    insertRegion = []
    mutNorm = []
    mutVaf = []
    coMuts = []
    for i in range(len(mutList)):
        pos = mutList[i].getInsertPos()
        assert pos > -1
        insertPos.append(pos)
        insertRegion.append(anno.iloc[pos]["region"])
        mutNorm.append(totCov[pos])
        mutVaf.append(100 * mutCount[i] / totCov[pos])
        
        cM_dict = {k: v for k, v in sorted(mutList[i].comutations.items(), key=lambda x: x[1], reverse = True)}
        cM_dict_pc = {}
        for k, v in zip(cM_dict.keys(), cM_dict.values()):
            cM_dict_pc[k] = f'{round(v * 100/ mutCount[i], 2)}%'
        coMuts.append(cM_dict_pc)

    dict = {'netInsert': netIns, 'counts': mutCount, 'vaf_percent': mutVaf,
        'coverage': mutNorm, 'insertPos': insertPos, 'insertRegion' : insertRegion,
        'name': mutName, 'HGVS': absmutName, 'co_mutations' : coMuts} 
    res = pd.DataFrame(dict)
    res.to_csv(config["MUTATION_FILE"], sep = ",", index=False)

    save_stats("\nTop inserts", config["STATS_FILE"])    
    
    # Filtered inserts
    res_s = res.copy(deep=True)
    res_s = res_s[res_s["netInsert"] >= config["MIN_INSERT_SEQ_LENGTH"]]
    res_s = res_s[res_s["counts"] >= config["MIN_TOTAL_READS"]]
    res_s = res_s[res_s["vaf_percent"] >= config["MIN_VAF"]]
    
    res_s = res_s[["netInsert", "counts", "vaf_percent", "insertPos", 
        "insertRegion", "coverage", "name", "HGVS", "co_mutations"]]
    res_s.to_csv(config["MUTATION_FILE_FILTERED"], sep = ",", index=False)
    
    res_s2 = res_s[["netInsert", "counts", "vaf_percent", "insertPos",
        "insertRegion", "coverage", "HGVS"]]
    res_txt = res_s2.head(10).to_string(index_names = False, index = False)
    save_stats(res_txt, config["STATS_FILE"])    

    # Amplicon-ome
    with open(config["OME_FILE"], 'w') as f:
        f.write(f'>WildType\n')
        f.write(f'{config["REF"]}\n')
        
    # Write each insert
    with open(config["OME_FILE"], 'a') as f:
        for i, ins_name in enumerate(res_s.name):
            f.write(f'>mutation_{i}\n')
            f.write(f'{generateMutSeq(config["REF"], [ins_name])}\n')
    
    ######################################################################
    # Write results    
    ######################################################################

    summa = res_s.groupby(["netInsert"])[["vaf_percent", "counts"]].sum().reset_index()
    summa = summa.sort_values(by = 'vaf_percent', ascending = False)
    save_stats("\nTop insert lengths", config["STATS_FILE"])    
    res_txt = summa.head(10).to_string(index_names = False, index = False)
    save_stats(res_txt, config["STATS_FILE"])    
    
    ######################################################################
    # Write results    
    summa.to_csv(config["NETINSERT_FILE"], sep = ",", index=False)
    ######################################################################
    
    tt = round(timeit.default_timer() - start_time, 2)
    save_stats(f'mergeITD time taken - {tt} sec', config["STATS_FILE"])    
    return 0

def parse_config_from_cmdline(config):
    """
    Get analysis parameters from commandline.

    Args:
        config (dict): Dict to save parameters to

    Returns:
        Filled config dict
    """
    parser = argparse.ArgumentParser()
    
    # Required parameters
    parser.add_argument("sampleID", help="sample ID used as output folder prefix (REQUIRED)")
    parser.add_argument("fastq1", help="FASTQ file (optionally gzipped) of forward reads (REQUIRED)")
    parser.add_argument("fastq2", help="FASTQ file (optionally gzipped) of reverse reads (REQUIRED)")
    
    # BBmap path
    parser.add_argument("-bbmap", help="Path to bbmap directory (default ~/bin/bbmap)", default="~/bin/bbmap", type=str)

    parser.add_argument("-reference", help="WT amplicon sequence as reference for read alignment (default ./anno/amplicon.txt)", default="./anno/amplicon.txt", type=str)
    parser.add_argument("-anno", help="WT amplicon sequence annotation (default ./anno/amplicon_kayser.tsv)", default="./anno/amplicon_kayser.tsv", type=str)
    
    parser.add_argument("-plot_coverage", help="If True, plot read coverage across the reference to 'coverage.png' in the respective output folder (default False)", default=False, type=str_to_bool)

    parser.add_argument("-progress_bar", help="If True, displays progress bar when aligning unique fragment sequences", default=True, type=str_to_bool)

    # BBMap alignment BAM
    parser.add_argument("-bam_from_reads", help="If True, uses BBmap to generate a sorted, cleaned BAM file from paired-end reads", default=False, type=str_to_bool)
    parser.add_argument("-bam_from_merged", help="If True, uses BBmap to generate a sorted, cleaned BAM file from merged reads", default=False, type=str_to_bool)
    parser.add_argument("-keep_merged_reads", help="If True, preserve a gzipped copy of merged reads after analysis", default=False, type=str_to_bool)
    parser.add_argument("-keep_temp_files", help="If True, preserve the temporary folder storing fastq intermediates", default=False, type=str_to_bool)
    
    parser.add_argument('-alignment_length_difference', help="filter out alignments with length difference from reference amplicon above this value, lower is more strict (default 20)", default="20", type=int)
    parser.add_argument('-alignment_mapq_filter', help="max number of bases of reference amplicon unaligned, lower is more strict (default 30)", default="30", type=int)
    
    # Not used
    parser.add_argument('-nkern', help="number of cores to use for parallel tasks (default 12)", default="12", type=int)
    
    # Alignment parameters
    parser.add_argument('-match', help="alignment cost of base match (default 5)", default="5", type=int)
    parser.add_argument('-mismatch', help="alignment cost of base mismatch (default -15)", default="-15", type=int)
    parser.add_argument('-gap_open', help="alignment cost of gap opening (default -36)", default="-36", type=int)
    parser.add_argument('-gap_extend', help="alignment cost of gap extension (default -0.5)", default="-0.5", type=float)
    
    # getHGVS() parameters
    parser.add_argument('-minAlignLen', help="minimum number of nucleotides that must be aligned in a block alignment (default 6)", default="6", type=int)
    parser.add_argument('-minRefAlignFraction', help="min fraction of reference that must be aligned by the read (default 0.4)", default="0.4", type=float)
    parser.add_argument('-maxFracIsIndel', help="max fraction of merged read that is allowed to be part of in/del -i.e. not aligned to reference. (default 0.7)", default="0.7", type=float)

    # bbmerge parameters
    parser.add_argument('-trimq_merging', help="Whether to attempt to merge reads a second time by first 3'-trimming reads by quality score cutoff prior (default = 20). -1 to disable", default="20", type=int)
    parser.add_argument("-min_bqs", help="minimum average base quality score (BQS) required by each read (default 25). -1 to disable", type=int, default=25)

    parser.add_argument('-min_read_copies', help="minimum number of copies of each read required for processing (1 to turn filter off, 2 (default) to discard unique reads)", default="2", type=int)
    parser.add_argument('-min_insert_seq_length', help="minimum number of insert basepairs which must be sequenced of each insert for it to be considered (default 6).", default="6", type=int)
    # parser.add_argument("-max_seq_Ns", help="maximum number of N's before these are filtered prior to alignment", type=int, default=-1)

    parser.add_argument('-filter_ins_total_reads', help="minimum number of total reads required to support an insertion for it to be considered (default 1)", default="1", type=int)
    parser.add_argument('-filter_ins_vaf', help="minimum variant allele frequency (VAF percent) required for an insertion to be considered 'high confidence' (default 0.006)", default="0.006", type=float)
    cmd_args = parser.parse_args()

    config["R1"] = cmd_args.fastq1
    config["R2"] = cmd_args.fastq2
    config["SAMPLE"] = cmd_args.sampleID
    config["NKERN"] = cmd_args.nkern

    config["REF_FILE"] = cmd_args.reference
    config["ANNO_FILE"] = cmd_args.anno
    
    config["BBMAP_PATH"] = cmd_args.bbmap
    config["BBMAP_TRIMQ"] = cmd_args.trimq_merging
    config["BBMAP_BQS"] = cmd_args.min_bqs
        
    # config["TECH"] = cmd_args.technology
    # if config["TECH"] == "454":
        # config["INFER_SENSE_FROM_ALIGNMENT"] = True
    # else:
        # config["INFER_SENSE_FROM_ALIGNMENT"] = cmd_args.infer_sense_from_alignment
    config["PLOT"] = cmd_args.plot_coverage
    config["PROGRESSBAR"] = cmd_args.progress_bar

    config["BAM_FROM_READS"] = cmd_args.bam_from_reads
    config["BAM_FROM_MERGED"] = cmd_args.bam_from_merged
    config["KEEP_MERGED"] = cmd_args.keep_merged_reads
    config["KEEP_TEMP"] = cmd_args.keep_temp_files

    config["ALIGN_LEN_DIFF"] = cmd_args.alignment_length_difference
    config["ALIGN_MAPQ"] = cmd_args.alignment_mapq_filter

    # R2 reads are reverse-complemented prior to alignment to the WT reference sequence
    # --> reverse-complement any sequence later to be found within reverse-complemented R2 reads
    # config["FORWARD_PRIMERS"] = [primer.upper() for primer in cmd_args.forward_primer]
    # config["REVERSE_PRIMERS"] = [primer.upper().translate(str.maketrans('ATCGatcg','TAGCtagc'))[::-1] for primer in cmd_args.reverse_primer]
    # config["FORWARD_ADAPTER"] = cmd_args.forward_adapter.upper()
    # config["REVERSE_ADAPTER"] = cmd_args.reverse_adapter.upper().translate(str.maketrans('ATCGatcg','TAGCtagc'))[::-1]

    config["COST_MATCH"] = cmd_args.match
    config["COST_MISMATCH"] = -abs(cmd_args.mismatch)
    config["COST_GAPOPEN"] = -abs(cmd_args.gap_open)
    config["COST_GAPEXTEND"] = -abs(cmd_args.gap_extend)

    config["MIN_ALIGN_LEN"] = cmd_args.minAlignLen
    config["MIN_REF_ALN_FRACTION"] = cmd_args.minRefAlignFraction
    config["MAX_FRAC_INDEL"] = cmd_args.maxFracIsIndel

    # config["MIN_SCORE_INSERTS"] = cmd_args.minscore_inserts
    # config["MIN_SCORE_ALIGNMENTS"] = cmd_args.minscore_alignments

    # config["MIN_BQS"] = cmd_args.min_bqs
    # config["MAX_NS"] = cmd_args.max_seq_Ns
    # config["MIN_READ_LENGTH"] = cmd_args.min_read_length
    config["MIN_INSERT_SEQ_LENGTH"] = cmd_args.min_insert_seq_length
    # if config["TECH"] == "454":
        # config["MIN_READ_COPIES"] = 1
    # else:
    config["MIN_READ_COPIES"] = cmd_args.min_read_copies
    # config["REQUIRE_INDEL_FREE_PRIMERS"] = cmd_args.require_indel_free_primers
    # config["MAX_TRAILING_BP"] = cmd_args.max_trailing_bp

    config["MIN_TOTAL_READS"] = cmd_args.filter_ins_total_reads
    config["MIN_VAF"] = cmd_args.filter_ins_vaf

    return config

def make_file_path_absolute(file_):
    if not os.path.isabs(file_):
        file_ = os.path.join(os.getcwd(), file_)
    return file_


def save_coverage(iref_coverage, config):
    """
    Write coverage distribution per inter-bp space
    to file `config["OUT_COV_FILE"]` in the `config["OUT_DIR"]` folder.

    Args:
        iref_coverage ([dict]): List oft three dictionaries which each contain
                the inter-bp coverage of the reference for i) forward reads only,
                ii) reverse reads only and iii) all reads, merged at the fragment
                level so that paired reads of the same DNA fragments are not counted
                twice at any given position.
        config (dict): Dictionary containing analysis parameters.
    """
    cov = pd.DataFrame(iref_coverage)
    cov.to_csv(config["OUT_COV_FILE"], sep="\t")


def plot_coverage(iref_coverage, config):
    """
    Plot coverage distribution per inter-bp space
    to file `config["OUT_COV_PLOT"]` in the `config["OUT_DIR"]` folder.

    Args:
        iref_coverage ([dict]): List oft three dictionaries which each contain
                the inter-bp coverage of the reference for i) forward reads only,
                ii) reverse reads only and iii) all reads, merged at the fragment
                level so that paired reads of the same DNA fragments are not counted
                twice at any given position.
        config (dict): Dictionary containing analysis parameters.
    """
    # import only when plotting is desired to avoid depending on matplotlib install?
    import matplotlib.pyplot as plt
    plt.switch_backend('Agg')

    fig, axs = plt.subplots(3, figsize=(20, 8), sharex=True, sharey=True)
    fig.suptitle("Final coverage achieved for " + config["SAMPLE"], fontsize=20)

    forward_plot = axs[0].bar(
            iref_coverage["all_reads"].keys(),
            iref_coverage["all_reads"].values(),
            label="total fragments",
            linewidth=0,
            width=1,
            color="dimgray")
    forward_plot = axs[1].bar(
            iref_coverage["forward_reads"].keys(),
            iref_coverage["forward_reads"].values(),
            label="forward reads",
            linewidth=0,
            width=1,
            color="tab:blue")
    forward_plot = axs[2].bar(
            iref_coverage["reverse_reads"].keys(),
            iref_coverage["reverse_reads"].values(),
            label="reverse reads",
            linewidth=0,
            width=1,
            color="tab:orange")

    for ax in axs:
        # Add some text for labels, title and custom x-axis tick labels, etc.
        ax.legend()

    axs[2].set_xlabel('reference bp', fontsize=18)
    axs[1].set_ylabel('# of reads aligned', fontsize=18)

    plt.tight_layout()
    plt.savefig(config["OUT_COV_PLOT"], dpi=300)


def main(config):

    # PROCESS INPUTS
    config["OUT_DIR"] = '_'.join([config["SAMPLE"], "mergeitd"])
    config["TMP_DIR"] = '_'.join([config["SAMPLE"], "mergeitd/temp_fastq"])

    config["OUT_COV_PLOT"] = os.path.join(config["OUT_DIR"], "coverage.png")
    config["OUT_COV_FILE"] = os.path.join(config["OUT_DIR"], "coverage.txt")
    config["STATS_FILE"] = os.path.join(config["OUT_DIR"], "stats.txt")
    config["CONFIG_FILE"] = os.path.join(config["OUT_DIR"], "config.txt")

    config["BBLOG"] = os.path.join(config["OUT_DIR"], "bbmap.log")
    config["ALIGN_FILE"] = os.path.join(config["OUT_DIR"], "alignClasses.csv")
    config["MUTATION_FILE"] = os.path.join(config["OUT_DIR"], "mutation_vaf.csv")
    config["MUTATION_FILE_FILTERED"] = os.path.join(config["OUT_DIR"], "filtered_mut_vaf.csv")
    config["NETINSERT_FILE"] = os.path.join(config["OUT_DIR"], "netInserts_vaf.csv")
    config["OME_FILE"] = os.path.join(config["OUT_DIR"], f'{config["SAMPLE"]}_ampliconome.fa')
    
    # make all input & output file / folder names absolute paths
    for file_ in ["R1", "R2", "REF_FILE", "ANNO_FILE", 
        "OUT_DIR", "TMP_DIR", "BBMAP_PATH",
        "OUT_COV_PLOT", "OUT_COV_FILE", "STATS_FILE", "CONFIG_FILE",
        "BBLOG", "ALIGN_FILE", "MUTATION_FILE", "MUTATION_FILE_FILTERED", "NETINSERT_FILE"
    ]:
        if config[file_]:
            config[file_] = make_file_path_absolute(config[file_])

    config["ANNO"] = read_annotation(config["ANNO_FILE"])
    config["ANNO"] = annotateCoords(config["ANNO"])
    
    # config["DOMAINS"] = get_domains(config["ANNO"])
    config["REF"] = read_reference(config["REF_FILE"]).upper()

    ## CREATE OUTPUT FOLDER
    if not os.path.exists(config["OUT_DIR"]):
        os.makedirs(config["OUT_DIR"])

    ## CREATE TEMP DIRECTORY FOR MERGED FASTQ
    if not os.path.exists(config["TMP_DIR"]):
        os.makedirs(config["TMP_DIR"])
        
    ## CHANGE TO OUTPUT FOLDER
    #  this is required for parallel child processes to retrieve
    #  the correct config.txt file later on despite static / constant filename
    # os.chdir(config["OUT_DIR"])
    save_config(config, config["CONFIG_FILE"])

    ## REMOVE OLD STATS & LOG FILE & START CREATING A NEW ONE
    try:
        os.remove(config["STATS_FILE"])
        # os.remove(os.path.join(config["OUT_DIR"], "incomplete-wt-tandem.log"))
    except OSError:
        pass
    save_stats(f'\n==== PROCESSING SAMPLE {config["SAMPLE"]} ====', config["STATS_FILE"])

    ### NEW MERGEITD PIPELINE

    config["ALIGNER"] = Align.PairwiseAligner()

    config["ALIGNER"].mode = 'global'
    config["ALIGNER"].match_score = config["COST_MATCH"]
    config["ALIGNER"].mismatch_score = config["COST_MISMATCH"]
    config["ALIGNER"].open_gap_score = config["COST_GAPOPEN"]
    config["ALIGNER"].extend_gap_score = config["COST_GAPEXTEND"]
    config["ALIGNER"].target_end_gap_score = 0.0
    config["ALIGNER"].query_end_gap_score = 0.0

    cleaned_fastq = bbmap_process(config)
    
    ### READS MERGED & CLEANED FASTQ READS
    readseq, _ = read_fastq(cleaned_fastq)
    # read_dict = {'Sequence' : readseq, 'BQS' : readbqs}
    # reads = pd.DataFrame(read_dict)
    
    ### GET UNIQUE READS
    unique_reads = Counter(readseq)
    # gb = reads.groupby(['Sequence'])
    # prealigns = gb.size().to_frame(name='Counts')
    # prealigns = prealigns.join(gb.agg({'BQS': 'mean'}).rename(columns={'BQS': 'avgBQS'}))
    
    ### MAKE PANDAS DF OF UNIQUE READS AND COUNTS
    prealigns = pd.DataFrame({
        "Sequence": list(Counter(unique_reads).keys()),
        "Counts" : list(Counter(unique_reads).values())
    })
    totalReads = prealigns["Counts"].sum()
    
    prealigns = prealigns.sort_values(by = "Counts", ascending = False).reset_index(drop = True)
    prealigns = prealigns[prealigns["Counts"] >= config["MIN_READ_COPIES"]]
    filteredReads = prealigns["Counts"].sum()
    pcReads = round(filteredReads * 100 / totalReads, 2)

    ### MEASURE SEQUENCE LENGTH
    prealigns["SeqLength"] = 0
    for i in range(len(prealigns)):
        prealigns.loc[i, "SeqLength"] = len(prealigns.iloc[i]["Sequence"])

    save_stats(f'\nAligning - {len(prealigns)} unique fragment sequences - {filteredReads} of {totalReads} total fragments ({pcReads} %)', config["STATS_FILE"])        
    alignITD(prealigns, config)

    # export cleaned fastq to main folder if required
    if config["KEEP_MERGED"]:            
        tmpFile = os.path.join(config["OUT_DIR"], "cleaned.fastq.gz")
        ret = subprocess.run(["gzip", cleaned_fastq])
        shutil.move(f'{cleaned_fastq}.gz', tmpFile)
        config["MERGED_READS"] = tmpFile
    elif config["BAM_FROM_MERGED"]:
        tmpFile = os.path.join(config["OUT_DIR"], "cleaned.fastq")
        shutil.move(cleaned_fastq, tmpFile)
        config["MERGED_READS"] = tmpFile         
    if not config["KEEP_TEMP"]:
        shutil.rmtree(config["TMP_DIR"])
    
    amp_names, amp_lens = get_amplicons(config)
    if len(amp_names) > 1:
        if config["BAM_FROM_READS"]:
            generate_bam(config, "PairedEnd")
        if config["BAM_FROM_MERGED"]:
            generate_bam(config, "Merged")
    elif config["BAM_FROM_READS"] or config["BAM_FROM_MERGED"]:
        save_stats(f'No inserts or ITDs detected, skipping BAM alignment steps', config["STATS_FILE"])    

    if not config["KEEP_MERGED"] and os.path.isfile(cleaned_fastq):
        os.remove(config["MERGED_READS"])

    ### END MERGEITD PIPELINE

    ########################################
    # CHANGE BACK TO ORIGINAL / PARENT DIRECTORY
    # os.chdir("..")

########## MAIN ####################
if __name__ == '__main__':

    config = parse_config_from_cmdline(config)
    main(config)
