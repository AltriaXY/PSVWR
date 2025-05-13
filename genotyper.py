import argparse
import os
import subprocess
import shutil
import sys
from multiprocessing import Manager
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed, wait, FIRST_COMPLETED
import time
import gzip
import re
import math
import mmap
import networkx as nx

def check_path(path):
    if os.path.exists(path):
        return True
    else:
        print(f"Path {path} does not exist.")
        return False

def index_gfa(gfa_file):
    '''
    Indexing gfa file for quick reading.
    :param gfa_file: pan-genome in gfa format.
    :param segments_in_bubble: segments belonging to large bubbles.
    :return: offset of lines in a list(offset_list).
    :return: line numbers of segments and links start from these segments in a dictionary(index_dict gfaID:[coverage, index, length]).
    '''
    offset_list = []
    index_dict = {}

    ref_seg2chr = {}
    chr2ref_seg = {}

    offset_byte = 0
    line_number = 0      #0 corresponding to the 1st line in gfa file, 1 corresponding to 2rd line in gfa file...

    genome_size = 0

    for eachline in gfa_file:
        offset_list.append(offset_byte)

        text = eachline.decode("utf-8").strip()
        if text[0] == "S":      #Segment line
            text_split = text.split("\t")
            gfaID = text_split[1]
            index_segment = line_number
            LN = len(text_split[2])
            index_dict[gfaID] = [{"seg":0}, [index_segment], LN]      #[coverage, index, length]
            if len(text_split) > 3:
                SR, SO, chr = -1, -1, ""
                i = 3
                for i in range(3, len(text_split)):
                    element = text_split[i]
                    if element[0:2] == "SR":
                        SR = int(element.split(":")[-1])
                    elif element[0:2] == "SO":
                        SO = int(element.split(":")[-1])
                    elif element[0:2] == "SN":
                        chr = element.split(":")[-1]
                if SR == 0:
                    if chr not in chr2ref_seg:
                        chr2ref_seg[chr] = {}
                    chr2ref_seg[chr][gfaID] = SO
                    ref_seg2chr[gfaID] = chr
                    genome_size += LN

        elif text[0] == "L":
            text_split = text.split("\t")
            gfaID1, gfaID2 = text_split[1], text_split[3]
            index_link = line_number
            index_dict[gfaID1][1].append(index_link)
            index_dict[gfaID1][0][gfaID2] = 0
            index_dict[gfaID2][1].append(index_link)
            index_dict[gfaID2][0][gfaID1] = 0

        offset_byte += len(eachline)
        line_number += 1

    offset_list.append(offset_byte)

    for chro in chr2ref_seg:
        chr2ref_seg[chro] = dict(sorted(chr2ref_seg[chro].items(), key = lambda item: item[1]))
    return offset_list, index_dict, ref_seg2chr, chr2ref_seg, genome_size

def process_alignment(alignment_list, index_dict):
    '''
    Processing result of alignment to pan-genome.
    :param alignment_list: an alignment record in list(splited with "\t")
    :param index_dict: output of index_gfa() ({gfaID: [coverage, index, length], ...})
    :return: coverage in dict(coverage_dict: {s1111:1, S2222: 0.458})
    :return: breakpoint in dict(breakpoint_list: [read_id, ref_path, start_on_ref_path, end_on_ref_path, deletion/insertion, None/[start_on_query, end_on_query, strand]])
    breakpoint的形式还需要思考，从CIGAR信息出发，两端的segment如何确定？如何确定断点在两端segment上位置？
    直接返回断点所在路径已经在路径（path的截取，包含<>）上的起始和终止？
    返回断点类型（del额外返回None，ins额外返回reads信息和插入序列在reads上的起始和终止）？
    '''
    start_on_read, end_on_read, path, path_length, start_on_path, end_on_path = \
        (int(alignment_list[2]), int(alignment_list[3]), alignment_list[5], int(alignment_list[6]),
         int(alignment_list[7]), int(alignment_list[8]))
    ##coverage
    coverage_dict = dict()
    seg = list()
    if ((">" in path) or ("<" in path)) and ((":" not in path) and ("-" not in path)):      #graphaligner
        seg = re.split(r"[<>]+", path)
        del seg[0]
        i = 0
        for i in range(len(seg)):
            seg_id ,coverage = seg[i], 0
            if i == 0:
                coverage = (index_dict[seg_id][2] - start_on_path)/index_dict[seg_id][2]
            elif i == len(seg) - 1:
                coverage = (index_dict[seg_id][2] - (path_length - end_on_path))/index_dict[seg_id][2]
            else:
                coverage = 1
            coverage_dict[seg_id] = coverage

    return coverage_dict, seg, end_on_read - start_on_read

def process_batch_size_alignment(several_alignment, index_dict):
    return [process_alignment(alignment_list, index_dict) for alignment_list in several_alignment]

def read_gaf(gaf_file, index_dict, genome_size, manager, thread):
    '''
    Reading gaf file, calculate coverage of each segment and find possible breakpoints.
    :param gaf_file: result of alignment to pan-genome in gaf format.
    :param thread: thread(int).
    :param index_dict: output of index_gfa().
    :return: write coverage into previous index_dict.
    :return: a new dict storing breakpoint(breakpoint).
    '''
    total_coverage_dict = {}
    total_coverage_link_dict = {}
    total_reads_length = 0

    share_index_dict = manager.dict(index_dict)
    batch_size = 10      #10 reads input into function process_batch_size_alignment() at a time
    batch = []
    futures = set()      #store tasks in a set
    MAX_PENDING = thread * 2      #max size of futures

    name = ""
    length = -1
    text_list = []

    with ProcessPoolExecutor(max_workers=thread) as executor:
        for eachline in gaf_file:
            text_split = eachline.strip().split("\t")
            na = text_split[0]
            le = int(text_split[3]) - int(text_split[2])
            if name == "":
                name, length, text_list = na, le, text_split
            elif name == na:
                if le > length:
                    length, text_list = le, text_split
            elif name != "" and name != na:
                batch.append(text_list)
                name, length, text_list = na, le, text_split

            if len(batch) == batch_size:
                future = executor.submit(process_batch_size_alignment, batch, share_index_dict)
                futures.add(future)
                batch = []

                if len(futures) >= MAX_PENDING:
                    done, futures = wait(futures, return_when = FIRST_COMPLETED)
                    for f in done:
                        results = f.result()
                        for result in results:
                            coverage_dict, segs, read_length = result
                            for key in coverage_dict:
                                if key not in total_coverage_dict:
                                    total_coverage_dict[key] = coverage_dict[key]
                                else:
                                    total_coverage_dict[key] += coverage_dict[key]
                            i = 0
                            for i in range(len(segs) - 1):
                                link = (segs[i], segs[i + 1])
                                if link not in total_coverage_link_dict:
                                    total_coverage_link_dict[link] = 1
                                else:
                                    total_coverage_link_dict[link] += 1
                            total_reads_length += read_length
        if text_list:
            batch.append(text_list)
        if batch:
            future = executor.submit(process_batch_size_alignment, batch, share_index_dict)
            futures.add(future)

        for f in wait(futures).done:
            results = f.result()
            for result in results:
                coverage_dict, segs, read_length = result
                for key in coverage_dict:
                    if key not in total_coverage_dict:
                        total_coverage_dict[key] = coverage_dict[key]
                    else:
                        total_coverage_dict[key] += coverage_dict[key]
                i = 0
                for i in range(len(segs) - 1):
                    link = (segs[i], segs[i + 1])
                    if link not in total_coverage_link_dict:
                        total_coverage_link_dict[link] = 1
                    else:
                        total_coverage_link_dict[link] += 1
                total_reads_length += read_length

    for key in total_coverage_dict:
        index_dict[key][0]["seg"] += total_coverage_dict[key]
    for key in total_coverage_link_dict:
        index_dict[key[0]][0][key[1]] += total_coverage_link_dict[key]
        index_dict[key[1]][0][key[0]] += total_coverage_link_dict[key]

    return total_reads_length / genome_size

def function_weight(coverage, x):
    '''
    Calcualting weight of edge from its coverage.
    :param coverage: from index_gfa
    :param x: sequencing depth, user input
    :return: weight(float)
    '''
    return math.exp(x - coverage)

def calculate_path_coverage(path, ref_d, alt_d):
    reads_length, path_length = 0, 0
    for seg in path:
        seg_coverage, seg_length = 0, 0
        if seg in ref_d:
            seg_coverage, seg_length = ref_d[seg][0]["seg"], ref_d[seg][2]
        elif seg in alt_d:
            seg_coverage, seg_length = alt_d[seg][0]["seg"], alt_d[seg][2]
        reads_length += seg_coverage * seg_length
        path_length += seg_length
    return reads_length/path_length, path_length

def extract_consecutive_groups(sorted_list):
    groups = []
    start = prev = sorted_list[0]

    for num in sorted_list[1:]:
        if num == prev + 1:
            prev = num
        else:
            groups.append(list(range(start, prev + 1)))
            start = prev = num
    groups.append(list(range(start, prev + 1)))
    return groups

def reverse_bases(bases):
    reverse_bases = ""
    length = len(bases)
    i = 0
    for i in range(length):
        base = bases[length-i-1]
        if base == "A":
            reverse_bases += "T"
        elif base == "T":
            reverse_bases += "A"
        elif base == "G":
            reverse_bases += "C"
        elif base == "C":
            reverse_bases += "G"
        else:
            reverse_bases += "N"
    return reverse_bases

def extract_path_beses_from_graph(path, graph):
    bases = ""
    if len(path) == 0:
        pass
    elif len(path) == 1:
        if graph.nodes[path[0]]["strand"] == "+":
            bases += graph.nodes[path[0]]["base"]
        else:
            bases += reverse_bases(graph.nodes[path[0]]["base"])
    else:
        for seg in path:
            if graph.nodes[seg]["strand"] == "+":
                bases += graph.nodes[seg]["base"]
            else:
                bases += reverse_bases(graph.nodes[seg]["base"])
    return bases

def extract_path_format_from_graph(path, graph):
    path_format = ""
    for seg in path:
        if graph.nodes[seg]["strand"] == "+":
            path_format += ">"
        else:
            path_format += "<"
        path_format += seg
    return path_format

def process_one_path(graph, path, bubble_chrom, shared_chr2ref_seg, start_seg, end_seg, ref_seg_dict, alt_seg_dict):
    vcf = list()
    tuple_ = tuple(shared_chr2ref_seg[bubble_chrom].keys())
    ref_path = tuple_[tuple_.index(start_seg): tuple_.index(end_seg) + 1]
    if ref_path == tuple(path):
        pass
    else:
        index_non_ref = list()
        index_ref = list()
        i = 0
        for i in range(len(path)):
            if path[i] in ref_path:
                index_ref.append(ref_path.index(path[i]))
            else:
                index_non_ref.append(i)
        if index_non_ref == list():
            index_ref_groups = extract_consecutive_groups(index_ref)
            group1, group2 = index_ref_groups[0], index_ref_groups[-1]
            i1, i2 = group1[-1], group2[0]
            s1, s2 = ref_path[i1], ref_path[i2]
            ref_p = ref_path[i1 + 1 : i2]
            alt_p = path[path.index(s1) + 1 : path.index(s2)]
            alt_p_complete = path[path.index(s1) : path.index(s2) + 1]
            POS = shared_chr2ref_seg[bubble_chrom][ref_p[0]]
            REF = graph.nodes[s1]["base"][-1] + extract_path_beses_from_graph(ref_p, graph)
            ALT = graph.nodes[s1]["base"][-1] + extract_path_beses_from_graph(alt_p, graph)
            SVPATH = ">{}>{}".format(s1, s2)
            SVTYPE = "DEL"
            SVLEN = len(ALT) - len(REF)
            if SVLEN <= -30:
                END = shared_chr2ref_seg[bubble_chrom][s2]
                COVER, _ = calculate_path_coverage(alt_p_complete, ref_seg_dict, alt_seg_dict)
                vcf.append((bubble_chrom, POS, "", REF, ALT, 255, ".", (SVPATH, SVTYPE, SVLEN, END, COVER), "GT", "1/."))
        else:
            i1, i2 = index_non_ref[0] - 1, index_non_ref[-1] + 1
            tmp1, tmp2 = index_ref[ : i1 + 1], index_ref[len(index_ref) - (len(path) - i2) : ]
            gs1, gs2 = extract_consecutive_groups(tmp1), extract_consecutive_groups(tmp2)
            if gs1 == []:
                pass
            elif len(gs1) == 1:
                pass
            else:
                i1 = path.index(ref_path[gs1[0][-1]])
            if gs2 == []:
                pass
            elif len(gs2) == 1:
                pass
            else:
                i2 = path.index(ref_path[gs2[-1][0]])
            s1, s2 = path[i1], path[i2]
            alt_p = path[i1 + 1 : i2]
            alt_p_complete = path[i1 : i2 + 1]
            ref_p = ref_path[ref_path.index(s1) + 1 : ref_path.index(s2)]
            if ref_p == tuple():
                POS = shared_chr2ref_seg[bubble_chrom][s2]
            else:
                POS = shared_chr2ref_seg[bubble_chrom][ref_p[0]]
            REF, ALT = (graph.nodes[s1]["base"][-1] + extract_path_beses_from_graph(ref_p, graph),
                        graph.nodes[s1]["base"][-1] + extract_path_beses_from_graph(alt_p, graph))
            SVPATH = extract_path_format_from_graph(alt_p_complete, graph)
            if len(REF) >= 30 or len(ALT) >= 30:
                if len(REF) > len(ALT):
                    SVTYPE = "DEL"
                elif len(REF) < len(ALT):
                    SVTYPE = "INS"
                else:
                    SVTYPE = "REPLACE"
                SVLEN = len(ALT) - len(REF)
                END = shared_chr2ref_seg[bubble_chrom][s2]
                COVER, _ = calculate_path_coverage(alt_p_complete, ref_seg_dict, alt_seg_dict)
                vcf.append((bubble_chrom, POS, "", REF, ALT, 255, ".", (SVPATH, SVTYPE, SVLEN, END, COVER), "GT", "1/."))
    return set(vcf)

'''
def process_one_path(graph, path, bubble_chrom, shared_chr2ref_seg, start_seg, end_seg, ref_seg_dict, alt_seg_dict):
    vcf = list()
    index_non_ref = list()
    index_ref = list()
    i = 0
    for i in range(len(path)):
        if path[i] in ref_seg_dict:
            index_ref.append(i)
        else:
            index_non_ref.append(i)

    tuple_ = tuple(shared_chr2ref_seg[bubble_chrom].keys())
    ref_path = tuple_[tuple_.index(start_seg): tuple_.index(end_seg) + 1]
    #check continuity of reference parts in this path
    index_ref_groups = extract_consecutive_groups(index_ref)
    for group in index_ref_groups:
        if len(group) == 1:
            pass
        else:
            j = 0
            for j in range(len(group) - 1):
                s1, s2 = path[group[j]], path[group[j+1]]
                i1, i2 = ref_path.index(s1) + 1, ref_path.index(s2)
                if i1 == i2:
                    pass
                else:
                    bases = graph.nodes[s1]["base"][-1]
                    bases += extract_path_beses_from_graph(ref_path[i1 : i2], graph)
                    SVLEN = 1 - len(bases)
                    if SVLEN <= -30:
                        pos = shared_chr2ref_seg[bubble_chrom][ref_path[i1]]
                        SVPATH = ">{}>{}".format(s1, s2)
                        SVTYPE = "DEL"
                        END = shared_chr2ref_seg[bubble_chrom][s2]
                        COVER, _ = calculate_path_coverage([s1, s2], ref_seg_dict, alt_seg_dict)
                        vcf.append(
                            (bubble_chrom, pos, "", bases, bases[0], 255, ".", (SVPATH, SVTYPE, SVLEN, END, COVER), "GT", "1/."))
        #other svs on branches
    if index_non_ref == list():
        pass
    else:
        index_groups= extract_consecutive_groups(index_non_ref)
        for group in index_groups:
            ref_seg_s, ref_seg_e = path[group[0] - 1], path[group[-1] + 1]
            alt_p = path[group[0] : group[-1] + 1]
            alt_p_complete = path[group[0] - 1 : group[-1] + 2]
            ref_p = ref_path[ref_path.index(ref_seg_s) + 1 : ref_path.index(ref_seg_e)]
            if ref_p == tuple():
                pos = shared_chr2ref_seg[bubble_chrom][ref_seg_e]
            else:
                pos = shared_chr2ref_seg[bubble_chrom][ref_p[0]]
            ref, alt = (graph.nodes[ref_seg_s]["base"][-1] + extract_path_beses_from_graph(ref_p, graph),
                        graph.nodes[ref_seg_s]["base"][-1] + extract_path_beses_from_graph(alt_p, graph))
            SVPATH = extract_path_format_from_graph(alt_p_complete, graph)
            if len(ref) >= 30 or len(alt) >= 30:
                if len(ref) > len(alt):
                    SVTYPE = "DEL"
                elif len(ref) < len(alt):
                    SVTYPE = "INS"
                else:
                    SVTYPE = "REPLACE"
                SVLEN = len(alt) - len(ref)
                END = shared_chr2ref_seg[bubble_chrom][ref_seg_e]
                COVER, _ = calculate_path_coverage(alt_p_complete, ref_seg_dict, alt_seg_dict)
                vcf.append(
                    (bubble_chrom, pos, "", ref, alt, 255, ".", (SVPATH, SVTYPE, SVLEN, END, COVER), "GT", "1/."))
    return set(vcf)
'''

def genotyping_one_bubble(bubble_chrom, segments, shared_ref_seg2chr, shared_chr2ref_seg,
                          shared_gfa_index_dict, shared_gfa_offset_list, gfa_dir, average_cover):
    cutoff = max(0.3 * average_cover, 3)
    vcf_list = list()
    ref_seg, alt_seg = dict(), dict()
    t = False
    for seg in segments:
        coverage_d, index_l, LN = shared_gfa_index_dict[seg][0], shared_gfa_index_dict[seg][1], shared_gfa_index_dict[seg][2]
        if seg in shared_ref_seg2chr:
            ref_seg[seg] = [coverage_d, index_l, LN]
            if coverage_d["seg"] < cutoff:
                t = True
        else:
            alt_seg[seg] = [coverage_d, index_l, LN]
            if coverage_d["seg"] >= cutoff:
                t = True
    if t == False:
        pass
    elif t == True:
        segments_remained = set(segments)  # set of segments whose edges haven't been added to graph
        segments_remained_2 = set(segments)  # set of segments whose strand haven't been determined
        G = nx.DiGraph()  # subgraph construction
        for seg in segments:
            G.add_node(seg, strand = "", base = "")
        with open(gfa_dir, "rb") as fi:
            mm = mmap.mmap(fi.fileno(), 0, access=mmap.ACCESS_READ)

            now_segments = []
            for ss in ref_seg:
                G.nodes[ss]["strand"] = "+"
                segments_remained_2.discard(ss)
            for ss in ref_seg:
                segments_remained.discard(ss)
                seg_text_start = mm[shared_gfa_offset_list[ref_seg[ss][1][0]] : shared_gfa_offset_list[ref_seg[ss][1][0] + 1]].decode("utf-8")
                seg_text_start_split = seg_text_start.strip().split("\t")
                G.nodes[ss]["base"] = seg_text_start_split[2]
                link_index_start = ref_seg[ss][1][1:]
                for line in link_index_start:
                    text = mm[shared_gfa_offset_list[line]: shared_gfa_offset_list[line + 1]].decode("utf-8")
                    text_split = text.strip().split("\t")
                    gfaID1, gfaID2, strand1, strand2 = text_split[1], text_split[3], text_split[2], text_split[4]
                    if gfaID1 == ss:
                        if strand1 == "+":
                            if gfaID2 in segments:
                                if gfaID2 in segments_remained_2:
                                    G.nodes[gfaID2]["strand"] = strand2
                                    segments_remained_2.discard(gfaID2)
                                    now_segments.append(gfaID2)
                                if not G.has_edge(gfaID1, gfaID2):
                                    G.add_edge(gfaID1, gfaID2,
                                               weight=function_weight(ref_seg[gfaID1][0][gfaID2], average_cover))
                        elif strand1 == "-":
                            if gfaID2 in segments:
                                if gfaID2 in segments_remained_2:
                                    if strand2 == "+":
                                        G.nodes[gfaID2]["strand"] = "-"
                                    elif strand2 == "-":
                                        G.nodes[gfaID2]["strand"] = "+"
                                    segments_remained_2.discard(gfaID2)
                                    now_segments.append(gfaID2)
                                if not G.has_edge(gfaID2, gfaID1):
                                    G.add_edge(gfaID2, gfaID1,
                                               weight=function_weight(ref_seg[gfaID1][0][gfaID2], average_cover))
                    elif gfaID2 == ss:
                        if strand2 == "+":
                            if gfaID1 in segments:
                                if gfaID1 in segments_remained_2:
                                    G.nodes[gfaID1]["strand"] = strand1
                                    segments_remained_2.discard(gfaID1)
                                    now_segments.append(gfaID1)
                                if not G.has_edge(gfaID1, gfaID2):
                                    G.add_edge(gfaID1, gfaID2,
                                               weight=function_weight(ref_seg[gfaID2][0][gfaID1], average_cover))
                        elif strand2 == "-":
                            if gfaID1 in segments:
                                if gfaID1 in segments_remained_2:
                                    if strand1 == "+":
                                        G.nodes[gfaID1]["strand"] = "-"
                                    elif strand1 == "-":
                                        G.nodes[gfaID1]["strand"] = "+"
                                    segments_remained_2.discard(gfaID1)
                                    now_segments.append(gfaID1)
                                if not G.has_edge(gfaID2, gfaID1):
                                    G.add_edge(gfaID2, gfaID1,
                                               weight=function_weight(ref_seg[gfaID2][0][gfaID1], average_cover))
            while segments_remained != set():
                now_segments_ts = []
                for seg in now_segments:
                    segments_remained.discard(seg)
                    seg_text = mm[shared_gfa_offset_list[alt_seg[seg][1][0]]: shared_gfa_offset_list[alt_seg[seg][1][0] + 1]].decode("utf-8")
                    seg_text_split = seg_text.strip().split("\t")
                    G.nodes[seg]["base"] = seg_text_split[2]
                    link_index = alt_seg[seg][1][1:]
                    for line in link_index:
                        text = mm[shared_gfa_offset_list[line]: shared_gfa_offset_list[line + 1]].decode("utf-8")
                        text_split = text.strip().split("\t")
                        gfaID1, gfaID2, strand1, strand2 = text_split[1], text_split[3], text_split[2], text_split[4]
                        if gfaID1 == seg:
                            if strand1 == G.nodes[gfaID1]["strand"]:
                                if gfaID2 in segments:
                                    if gfaID2 in segments_remained_2:
                                        G.nodes[gfaID2]["strand"] = strand2
                                        segments_remained_2.discard(gfaID2)
                                        now_segments_ts.append(gfaID2)
                                    if not G.has_edge(gfaID1, gfaID2):
                                        G.add_edge(gfaID1, gfaID2,
                                                   weight=function_weight(alt_seg[gfaID1][0][gfaID2], average_cover))
                            elif strand1 != G.nodes[gfaID1]["strand"]:
                                if gfaID2 in segments:
                                    if gfaID2 in segments_remained_2:
                                        if strand2 == "+":
                                            G.nodes[gfaID2]["strand"] = "-"
                                        elif strand2 == "-":
                                            G.nodes[gfaID2]["strand"] = "+"
                                        segments_remained_2.discard(gfaID2)
                                        now_segments_ts.append(gfaID2)
                                    if not G.has_edge(gfaID2, gfaID1):
                                        G.add_edge(gfaID2, gfaID1,
                                                   weight=function_weight(alt_seg[gfaID1][0][gfaID2], average_cover))
                        elif gfaID2 == seg:
                            if strand2 == G.nodes[gfaID2]["strand"]:
                                if gfaID1 in segments:
                                    if gfaID1 in segments_remained_2:
                                        G.nodes[gfaID1]["strand"] = strand1
                                        segments_remained_2.discard(gfaID1)
                                        now_segments_ts.append(gfaID1)
                                    if not G.has_edge(gfaID1, gfaID2):
                                        G.add_edge(gfaID1, gfaID2,
                                                   weight=function_weight(alt_seg[gfaID2][0][gfaID1], average_cover))
                            elif strand2 != G.nodes[gfaID2]["strand"]:
                                if gfaID1 in segments:
                                    if gfaID1 in segments_remained_2:
                                        if strand1 == "+":
                                            G.nodes[gfaID1]["strand"] = "-"
                                        elif strand1 == "-":
                                            G.nodes[gfaID1]["strand"] = "+"
                                        segments_remained_2.discard(gfaID1)
                                        now_segments_ts.append(gfaID1)
                                    if not G.has_edge(gfaID2, gfaID1):
                                        G.add_edge(gfaID2, gfaID1,
                                                   weight=function_weight(alt_seg[gfaID2][0][gfaID1], average_cover))
                now_segments = list(now_segments_ts)
            mm.close()

        start_seg, end_seg = segments[0], segments[-1]
        if start_seg not in G or end_seg not in G:
            pass
        else:
            paths = nx.shortest_simple_paths(G, start_seg, end_seg, weight="weight")
            first_path, first_path_coverage = next(paths, None), -1
            second_path, second_path_coverage = list(), -1
            first_path_coverage, first_path_length = calculate_path_coverage(first_path, ref_seg, alt_seg)
            genetic_constitution = -1  # 0 = Homo, 1 = Hetero, -1 = reads coverage is not enough
            if first_path_coverage < max(0.3 * average_cover, 3):
                pass
            elif first_path_coverage >= 0.7 * average_cover:
                genetic_constitution = 0
            else:
                second_path, second_path_coverage = next(paths, None), -1
                if second_path == None:
                    genetic_constitution = 0
                else:
                    second_path_coverage, second_path_length = calculate_path_coverage(second_path, ref_seg, alt_seg)
                    two_paths_intersection = set(first_path)
                    two_paths_intersection.intersection(set(second_path))
                    _, two_paths_intersection_length = calculate_path_coverage(two_paths_intersection, ref_seg, alt_seg)
                    if second_path_coverage < max(0.3 * average_cover, 3):
                        genetic_constitution = 0
                    else:
                        if (two_paths_intersection_length / first_path_length >= 0.9) or (two_paths_intersection_length / second_path_length >= 0.9):
                            genetic_constitution = 0
                        else:
                            genetic_constitution = 1
            if genetic_constitution == -1:
                pass
            elif genetic_constitution == 0:
                first_path = tuple(first_path)
                l1 = process_one_path(G, first_path, bubble_chrom, shared_chr2ref_seg, start_seg, end_seg,
                                            ref_seg, alt_seg)
                for vcf in l1:
                    vcf_list.append((vcf[0], vcf[1], vcf[2], vcf[3], vcf[4], vcf[5], vcf[6],
                                     vcf[7], vcf[8], "1/1"))
            elif genetic_constitution == 1:
                first_path, second_path = tuple(first_path), tuple(second_path)
                l1 = process_one_path(G, first_path, bubble_chrom, shared_chr2ref_seg, start_seg, end_seg,
                                      ref_seg, alt_seg)
                l2 = process_one_path(G, second_path, bubble_chrom, shared_chr2ref_seg, start_seg, end_seg,
                                      ref_seg, alt_seg)
                for vcf in l2:
                    if vcf in l1:
                        vcf_list.append((vcf[0], vcf[1], vcf[2], vcf[3], vcf[4], vcf[5], vcf[6],
                                         vcf[7], vcf[8], "1/1"))
                        l1.discard(vcf)
                    else:
                        vcf_list.append((vcf[0], vcf[1], vcf[2], vcf[3], vcf[4], vcf[5], vcf[6],
                                         vcf[7], vcf[8], "0/1"))
                for vcf in l1:
                    vcf_list.append((vcf[0], vcf[1], vcf[2], vcf[3], vcf[4], vcf[5], vcf[6],
                                     vcf[7], vcf[8], "0/1"))
    return vcf_list

def genotyping_several_bubbles(several_bubbles_list, shared_ref_seg2chr, shared_chr2ref_seg,
                               shared_gfa_index_dict, shared_gfa_offset_list, gfa_dir, average_cover):
    return [genotyping_one_bubble(*args, shared_ref_seg2chr, shared_chr2ref_seg,
                                  shared_gfa_index_dict, shared_gfa_offset_list,
                                  gfa_dir, average_cover) for args in several_bubbles_list]
def genotyping(bubble_dict, shared_ref_seg2chr, shared_chr2ref_seg, shared_gfa_index_dict, shared_gfa_offset_list,
               gfa_dir, average_cover, thread):
    variants_in_large_bubbles = list()
    batch_size = 10
    batch = list([] for i in range(thread))
    futures = set()
    max_pending = thread * 2
    with ProcessPoolExecutor(max_workers=thread) as executor:
        count = 0
        order = 0
        for chrom in bubble_dict:
            for start in bubble_dict[chrom]:
                arguments = list()
                arguments.append(chrom)
                arguments.append(bubble_dict[chrom][start][1])
                batch[order].append(arguments)
                order = (order + 1) % thread
                if len(batch[thread - 1]) == batch_size:
                    for ttt in batch:
                        future = executor.submit(genotyping_several_bubbles, ttt, shared_ref_seg2chr,
                                             shared_chr2ref_seg, shared_gfa_index_dict, shared_gfa_offset_list,
                                             gfa_dir, average_cover)
                        futures.add(future)
                        if len(futures) >= max_pending:
                            done, futures = wait(futures, return_when=FIRST_COMPLETED)
                            for f in done:
                                results = f.result()
                                for result in results:
                                    variants_in_large_bubbles.extend(result)
                                count += 1
                                print(count)
                    batch = list([] for i in range(thread))
        print("final")
        for f in wait(futures).done:
            results = f.result()
            for result in results:
                variants_in_large_bubbles.extend(result)
            count += 1
            print(count)
    sorted_variants_in_large_bubbles = sorted(variants_in_large_bubbles, key = lambda x: (x[0], x[1]))
    return sorted_variants_in_large_bubbles

def main():
    parser = argparse.ArgumentParser(description="This program is writen for genotyping large bubbles" \
                                                 "based on long reads alignment result and pan-genome graph." \
                                                 "")
    parser.add_argument("-i", "--input", type=str, help="alignment to pan-genome in gaf format")
    parser.add_argument("-g", "--gfa_in", type=str, help="input pan-genome graph in gfa format")
    parser.add_argument("-o", "--output", type=str, help="output_SV")
    parser.add_argument("-t", "--thread", type=int, default=1, help="thread(default=1)")
    args = parser.parse_args()

    path_to_gfatools = shutil.which("gfatools")
    if path_to_gfatools == None:
        print("Can not find gfatools.")
        sys.exit(1)
    os.environ["PATH"] += (
            os.pathsep + path_to_gfatools
    )

    if args.input:
        pass
    else:
        print("Need gaf input, see -i/--input.")
    check_path(args.input)

    if args.gfa_in:
        pass
    else:
        print("Need graph(gfa format) input, see -g/--gfa_in.")
    check_path(args.gfa_in)

    if args.output:
        pass
    else:
        print("Need directory for output, see -o/--output.")
    check_path(args.output)

    gaf = args.input  # input gaf
    gfa = args.gfa_in  # input gfa
    out_dir = args.output  # output directory
    thread = args.thread  # thread

    #bubble
    time1 = time.time()
    print("Calling bubbles")

    bubble, gfaID2bubble = dict(), dict()  #

    with subprocess.Popen(["gfatools", "bubble", gfa], stderr=subprocess.PIPE, stdout=subprocess.PIPE,
                          universal_newlines=True) as proc:
        for line in proc.stdout:
            info = line.strip().split("\t")
            chr, start, end, segments = info[0], int(info[1]), int(info[2]), tuple(info[11].split(","))
            ref, alt = info[12], info[13]
            if ((alt != "*") and (len(alt) >= 50)) or ((ref != "*") and (len(ref) >= 50)):
                if chr in bubble:
                    pass
                elif chr not in bubble:
                    bubble[chr] = dict()
                bubble[chr][start] = [end, segments]
                for segment in segments:
                    if segment not in gfaID2bubble:
                        gfaID2bubble[segment] = list()
                    if (chr, start) not in gfaID2bubble[segment]:
                        gfaID2bubble[segment].append((chr, start))  # In case that some segments may belong to two bubbles, because they are the linkages between two bubbles.
            else:
                pass
    for chrom in bubble:
        bubble[chrom] = dict(sorted(bubble[chrom].items(), key = lambda x: -len(x[1][1])))
        print("{} has {} large bubbles".format(chrom, len(bubble[chrom])))

    print("Finish calling bubbles in {}s".format(time.time() - time1))

    #gfa
    time2 = time.time()
    print("indexing gfa.")

    gfa_offset_list, gfa_index_dict, reference_seg2chr, chr2reference_seg, genome_size = list(), dict(), dict(), dict(), 0
    if gfa[len(gfa) - 7:] == ".gfa.gz":
        with gzip.open(gfa, "rb") as gfa_file:
            gfa_offset_list, gfa_index_dict, reference_seg2chr, chr2reference_seg, genome_size = index_gfa(gfa_file)
    elif gfa[len(gfa) - 4:] == ".gfa":
        with open(gfa, "rb") as gfa_file:
            gfa_offset_list, gfa_index_dict, reference_seg2chr, chr2reference_seg, genome_size = index_gfa(gfa_file)
    else:
        print("wrong gfa format")

    print("Finish indexing gfa in {}s.".format(time.time() - time2))
    print("Genome size: {} bp".format(genome_size))
    average_cover = 0
    with Manager() as manager:
        #gaf
        time3 = time.time()
        print("reading gaf.")

        if gaf[len(gaf) - 7:] == ".gaf.gz":
            with gzip.open(gaf, "rt") as gaf_file:
                average_cover = read_gaf(gaf_file, gfa_index_dict, genome_size, manager, thread)
        elif gaf[len(gaf) - 4:] == ".gaf":
            with open(gaf, "rt") as gaf_file:
                average_cover = read_gaf(gaf_file, gfa_index_dict, genome_size, manager, thread)

        print("Finish reading gaf in {}s.".format(time.time() - time3))
        print("average_cover: {}X".format(average_cover))

        #genotyping
        time4 = time.time()
        print("Genotyping.")

        vcf_info = list()

        shared_reference_seg2chr = manager.dict(reference_seg2chr)
        del reference_seg2chr
        shared_chr2reference_seg = manager.dict()
        for k, v in chr2reference_seg.items():
            if isinstance(v, dict):
                shared_chr2reference_seg[k] = manager.dict(v)
            else:
                shared_chr2reference_seg[k] = v
        del chr2reference_seg
        shared_gfa_index_dict = manager.dict(gfa_index_dict)
        del gfa_index_dict
        shared_gfa_offset_list = manager.list(gfa_offset_list)
        del gfa_offset_list

        vcf_info = genotyping(bubble, shared_reference_seg2chr, shared_chr2reference_seg, shared_gfa_index_dict,
                              shared_gfa_offset_list, gfa, average_cover, thread)
        del shared_gfa_offset_list, shared_reference_seg2chr

        print("Finish genotyping in {}s.".format(time.time() - time4))

        #out
        time5 = time.time()
        print("Outputing vcf.")

        out_vcf_dir = os.path.join(out_dir, "NEW_SV.vcf")
        with open(out_vcf_dir, "w") as fo_vcf:
            fo_vcf.write("##fileformat=VCFv4.2\n")
            fo_vcf.write("##source=LargeBubbleGenotyper\n")
            fo_vcf.write("##reference={}\n".format(gfa))
            fo_vcf.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
            fo_vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            fo_vcf.write(
                '##INFO=<ID=SVPATH,Number=1,Type=String,Description="ALT path in pan-genome graph.">\n')
            fo_vcf.write(
                '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant, INS/DEL/REPLACE.">\n')
            fo_vcf.write(
                '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="ALT length minus REF length.">\n')
            fo_vcf.write(
                '##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate of structural variant on reference.">\n')
            fo_vcf.write(
                '##INFO=<ID=COVER,Number=1,Type=Float,Description="Coverage of structural variant.">\n')
            for chromsome in shared_chr2reference_seg:
                final_seg = next(reversed(list(shared_chr2reference_seg[chromsome].keys())))
                so, length_final_seg = shared_chr2reference_seg[chromsome][final_seg], shared_gfa_index_dict[final_seg][2]
                fo_vcf.write('##contig=<ID={},length={}>\n'.format(chromsome, so + length_final_seg))
                del shared_gfa_index_dict, shared_chr2reference_seg
                fo_vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
            for element in vcf_info:
                CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO_TUPLE, FORMAT, SAMPLE = element
                ID = "{}_{}_{}_{}".format(CHROM, POS, INFO_TUPLE[3], INFO_TUPLE[1])
                fo_vcf.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\tSVPATH={};SVTYPE={};SVLEN={};END={};COVER={}\t{}\t{}\n".format(
                    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO_TUPLE[0], INFO_TUPLE[1], INFO_TUPLE[2], INFO_TUPLE[3],
                    INFO_TUPLE[4], FORMAT, SAMPLE
                ))

        print("Finish outputing vcf in {}s.".format(time.time() - time5))
        print("Finished.")

if __name__ == "__main__":
    main()






