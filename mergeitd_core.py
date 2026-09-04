"""Verbatim mergeITD alignment/calling core extracted from mergeitd.py.

Do not edit these functions independently of mergeitd.py.  New BAM orchestration
lives in itdmapper.py / mergeitd_runner.py.
"""

import copy
import timeit

import numpy as np
import pandas as pd
from tqdm import trange


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
