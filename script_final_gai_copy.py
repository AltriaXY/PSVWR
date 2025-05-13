import argparse
import copy
import os
import sys
import tempfile
import shutil
import subprocess
from multiprocessing import Manager
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed, wait, FIRST_COMPLETED
import time
import gzip
import re
import numpy as np
import math
import mmap
import networkx as nx
from scipy.stats import binom
import logging



######read gfa file
def index_gfa(gfa_file):
    '''
    Indexing gfa file for quick reading.
    :param gfa_file: pan-genome in gfa format.
    :return: offset of lines in a list(offset_list).
    :return: line numbers of segments and links start from these segments in a dictionary(index_dict gfaID:[coverage, index, length]).
    '''
    offset_list = []
    index_dict = {}

    ref_seg2chr = {}
    chr2ref_seg = {}

    offset_byte = 0
    line_number = 0      #0 corresponding to the 1st line in gfa file, 1 corresponding to 2rd line in gfa file...

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
    return offset_list, index_dict, ref_seg2chr, chr2ref_seg



######call bubbles



######Reading gaf
def count_query_and_ref_length_from_cigar(op, num, start, end, query_or_ref):
    '''
    Counting query length and reference length depending on CIGAR
    :param op:
    :param num:
    :param start: close
    :param end: open
    :param query_or_ref:
    :return:
    '''
    if query_or_ref == "query":
        query_length = 0
        i = 0
        for i in range(start, end):
            if op[i] in {"M", "X", "="}:
                query_length += int(num[i])
            elif op[i] in {"S", "I"}:
                query_length += int(num[i])
        return query_length
    elif query_or_ref == "ref":
        ref_length = 0
        i = 0
        for i in range(start, end):
            if op[i] in {"M", "X", "="}:
                ref_length += int(num[i])
            elif op[i] in {"N", "D"}:
                ref_length += int(num[i])
        return ref_length

def parse_cigar(cigar):
    num_list = []
    op_list = []
    num = ''
    for ch in cigar:
        if ch.isdigit():
            num += ch
        else:
            num_list.append(int(num))
            op_list.append(ch)
            num = ''
    return num_list, op_list

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
    readID, start_on_read, end_on_read, strand, path, path_length, start_on_path, end_on_path, cigar = \
        alignment_list[0], int(alignment_list[2]), int(alignment_list[3]), alignment_list[4], alignment_list[5], int(alignment_list[6]), int(alignment_list[7]), int(alignment_list[8]), ""
    for info in alignment_list:
        if info.startswith("cg:Z:"):
            cigar = info.split(":")[2]
            break
    ##coverage
    coverage_dict = {}
    seg = []
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
    ##breakpoint
    breakpoint_list_del, breakpoint_list_ins = [], []

    cigar_num_list, cigar_op_list = parse_cigar(cigar)
    possible_sv_del, possible_sv_ins = list(), list()      #posible_sv_del(>=10D) and their index in cigar
    k = 0
    for k in range(len(cigar_op_list)):
        if cigar_op_list[k] == "D":
            if int(cigar_num_list[k]) >= 10:
                possible_sv_del.append(k)
        elif cigar_op_list[k] == "I":
            if int(cigar_num_list[k]) >= 10:
                possible_sv_ins.append(k)
    #deletion
    i = 0
    while i <= len(possible_sv_del) - 1:
        k_start = possible_sv_del[i]
        length_del = int(cigar_num_list[possible_sv_del[i]])
        while (i < len(possible_sv_del) - 1) and (possible_sv_del[i] + 2 == possible_sv_del[i+1]) and (int(cigar_num_list[possible_sv_del[i] + 1]) == 1):
            i += 1
            length_del += int(cigar_num_list[possible_sv_del[i]])
        k_end = possible_sv_del[i]
        if length_del < 30:
            pass
        elif length_del >= 30:
            ref_length1 = count_query_and_ref_length_from_cigar(cigar_op_list, cigar_num_list, 0, k_start, "ref")
            ref_length2 = count_query_and_ref_length_from_cigar(cigar_op_list, cigar_num_list, k_start, k_end + 1, "ref")
            direction = re.findall(r"[<>]+", path)
            length = 0
            ref_path, start_on_ref_path, end_on_ref_path = "", -1, -1
            j = 0
            j_start = j_end = -1
            for j in range(len(seg)):
                le_seg = index_dict[seg[j]][2]
                if j == 0:
                    if ref_length1 >= length and ref_length1 < length + le_seg - start_on_path:
                        j_start = j
                        start_on_ref_path = ref_length1 - length
                    if ref_length1 + ref_length2 >= length and ref_length1 + ref_length2 < length + le_seg - start_on_path:
                        j_end = j
                        end_on_ref_path = start_on_ref_path + ref_length2
                        break
                    length += le_seg - start_on_path
                elif j > 0 and j < len(seg) - 1:
                    if ref_length1 >= length and ref_length1 < length + le_seg:
                        j_start = j
                        start_on_ref_path = ref_length1 - length
                    if ref_length1 + ref_length2 >= length and ref_length1 + ref_length2 < length + le_seg:
                        j_end = j
                        end_on_ref_path = start_on_ref_path + ref_length2
                        break
                    length += le_seg
                elif j == len(seg) - 1:
                    if ref_length1 >= length and ref_length1 <= length + le_seg - (path_length - end_on_path):
                        j_start = j
                        start_on_ref_path = ref_length1 - length
                    if ref_length1 + ref_length2 >= length and ref_length1 + ref_length2 <= length + le_seg - (path_length - end_on_path):
                        j_end = j
                        end_on_ref_path = start_on_ref_path + ref_length2
                        break
            m = 0
            for m in range(j_start, j_end + 1):
                ref_path += "{}{}".format(direction[m], seg[m])

            breakpoint_list_del.append([readID, ref_path, start_on_ref_path, end_on_ref_path, "deletion", None])

        i += 1
    #insertion
    i = 0
    while i <= len(possible_sv_ins) - 1:
        k_start = possible_sv_ins[i]
        length_ins = int(cigar_num_list[possible_sv_ins[i]])
        while (i < len(possible_sv_ins) - 1) and (possible_sv_ins[i] + 2 == possible_sv_ins[i+1]) and (int(cigar_num_list[possible_sv_ins[i] + 1]) == 1):
            i += 1
            length_ins += int(cigar_num_list[possible_sv_ins[i]])
        k_end = possible_sv_ins[i]
        if length_ins < 30:
            pass
        elif length_ins >= 30:
            ref_length1 = count_query_and_ref_length_from_cigar(cigar_op_list, cigar_num_list, 0, k_start, "ref")
            ref_length2 = count_query_and_ref_length_from_cigar(cigar_op_list, cigar_num_list, k_start, k_end + 1,"ref")
            query_length1 = count_query_and_ref_length_from_cigar(cigar_op_list, cigar_num_list, 0, k_start, "query")
            query_length2 = count_query_and_ref_length_from_cigar(cigar_op_list, cigar_num_list, k_start, k_end + 1, "query")
            direction = re.findall(r"[<>]+", path)
            length = 0
            ref_path, start_on_ref_path, end_on_ref_path, start_on_query, end_on_query = "", -1, -1, -1, -1
            if strand == "+":
                start_on_query, end_on_query = start_on_read + query_length1, start_on_read + query_length1 + query_length2
            elif strand == "-":
                start_on_query, end_on_query = end_on_read - query_length1 - query_length2, end_on_read - query_length1
            j = 0
            j_start = j_end = -1
            for j in range(len(seg)):
                le_seg = index_dict[seg[j]][2]
                if j == 0:
                    if ref_length1 >= length and ref_length1 < length + le_seg - start_on_path:
                        j_start = j
                        start_on_ref_path = ref_length1 - length
                    if ref_length1 + ref_length2 >= length and ref_length1 + ref_length2 < length + le_seg - start_on_path:
                        j_end = j
                        end_on_ref_path = start_on_ref_path + ref_length2
                        break
                    length += le_seg - start_on_path
                elif j > 0 and j < len(seg) - 1:
                    if ref_length1 >= length and ref_length1 < length + le_seg:
                        j_start = j
                        start_on_ref_path = ref_length1 - length
                    if ref_length1 + ref_length2 >= length and ref_length1 + ref_length2 < length + le_seg:
                        j_end = j
                        end_on_ref_path = start_on_ref_path + ref_length2
                        break
                    length += le_seg
                elif j == len(seg) - 1:
                    if ref_length1 >= length and ref_length1 <= length + le_seg - (path_length - end_on_path):
                        j_start = j
                        start_on_ref_path = ref_length1 - length
                    if ref_length1 + ref_length2 >= length and ref_length1 + ref_length2 <= length + le_seg - (path_length - end_on_path):
                        j_end = j
                        end_on_ref_path = start_on_ref_path + ref_length2
                        break
            m = 0
            for m in range(j_start, j_end + 1):
                ref_path += "{}{}".format(direction[m], seg[m])

            breakpoint_list_ins.append([readID, ref_path, start_on_ref_path, end_on_ref_path, "insertion", (start_on_query, end_on_query, strand)])

        i += 1

    return coverage_dict, breakpoint_list_del, breakpoint_list_ins, seg, end_on_read - start_on_read

def process_batch_size_alignment(several_alignment, index_dict):
    return [process_alignment(alignment_list, index_dict) for alignment_list in several_alignment]

def read_gaf(gaf_file, thread, index_dict, manager):
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
    total_breakpoint_list_del, total_breakpoint_list_ins = [], []

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
                            coverage_dict, breakpoint_list_del, breakpoint_list_ins, segs, read_length = result
                            for key in coverage_dict:
                                if key not in total_coverage_dict:
                                    total_coverage_dict[key] = coverage_dict[key]
                                else:
                                    total_coverage_dict[key] += coverage_dict[key]
                            if breakpoint_list_del:
                                for bp in breakpoint_list_del:
                                    total_breakpoint_list_del.append(bp)
                            if breakpoint_list_ins:
                                for bp in breakpoint_list_ins:
                                    total_breakpoint_list_ins.append(bp)
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
                coverage_dict, breakpoint_list_del, breakpoint_list_ins, segs, read_length = result
                for key in coverage_dict:
                    if key not in total_coverage_dict:
                        total_coverage_dict[key] = coverage_dict[key]
                    else:
                        total_coverage_dict[key] += coverage_dict[key]
                if breakpoint_list_del:
                    for bp in breakpoint_list_del:
                        total_breakpoint_list_del.append(bp)
                if breakpoint_list_ins:
                    for bp in breakpoint_list_ins:
                        total_breakpoint_list_ins.append(bp)
                i = 0
                for i in range(len(segs) - 1):
                    if (segs[i], segs[i + 1]) not in total_coverage_link_dict:
                        total_coverage_link_dict[(segs[i], segs[i + 1])] = 1
                    else:
                        total_coverage_link_dict[(segs[i], segs[i + 1])] += 1
                total_reads_length += read_length

    for key in total_coverage_dict:
        index_dict[key][0]["seg"] += total_coverage_dict[key]
    for key in total_coverage_link_dict:
        index_dict[key[0]][0][key[1]] += total_coverage_link_dict[key]
        index_dict[key[1]][0][key[0]] += total_coverage_link_dict[key]

    return total_breakpoint_list_del, total_breakpoint_list_ins, total_reads_length / 3.05e9

"""
def read_gaf(gaf_file, thread, index_dict, manager):
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
    total_breakpoint_list_del, total_breakpoint_list_ins = [], []

    total_reads_length = 0

    share_index_dict = manager.dict(index_dict)

    batch_size = thread
    batch = [[] for _ in range(batch_size)]

    name = ""
    length = -1
    text_list = []

    with ProcessPoolExecutor(max_workers=thread) as executor:
        order = 0
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
                batch[order].append(text_list)
                order = (order + 1) % batch_size
                name, length, text_list = na, le, text_split

            if len(batch[batch_size - 1]) == 25:
                futures = [executor.submit(process_batch_size_alignment, several_text_list, share_index_dict) for several_text_list in batch]
                for future in as_completed(futures):
                    results = future.result()
                    for result in results:
                        coverage_dict, breakpoint_list_del, breakpoint_list_ins, segs, read_length = result
                        for key in coverage_dict:
                            if key not in total_coverage_dict:
                                total_coverage_dict[key] = coverage_dict[key]
                            else:
                                total_coverage_dict[key] += coverage_dict[key]
                        if breakpoint_list_del:
                            for bp in breakpoint_list_del:
                                total_breakpoint_list_del.append(bp)
                        if breakpoint_list_ins:
                            for bp in breakpoint_list_ins:
                                total_breakpoint_list_ins.append(bp)
                        i = 0
                        for i in range(len(segs)-1):
                            if (segs[i], segs[i+1]) not in total_coverage_link_dict:
                                total_coverage_link_dict[(segs[i], segs[i+1])] = 1
                            else:
                                total_coverage_link_dict[(segs[i], segs[i + 1])] += 1
                        total_reads_length += read_length
                batch = [[] for _ in range(batch_size)]
        if text_list:
            batch[order].append(text_list)
        if batch:
            futures = [executor.submit(process_batch_size_alignment, several_text_list, share_index_dict) for several_text_list in batch if several_text_list]
            for future in as_completed(futures):
                results = future.result()
                for result in results:
                    coverage_dict, breakpoint_list_del, breakpoint_list_ins, segs, read_length = result[0], result[1], result[2], result[3], result[4]
                    for key in coverage_dict:
                        if key not in total_coverage_dict:
                            total_coverage_dict[key] = coverage_dict[key]
                        else:
                            total_coverage_dict[key] += coverage_dict[key]
                    if breakpoint_list_del:
                        for bp in breakpoint_list_del:
                            total_breakpoint_list_del.append(bp)
                    if breakpoint_list_ins:
                        for bp in breakpoint_list_ins:
                            total_breakpoint_list_ins.append(bp)
                    i = 0
                    for i in range(len(segs) - 1):
                        if (segs[i], segs[i + 1]) not in total_coverage_link_dict:
                            total_coverage_link_dict[(segs[i], segs[i + 1])] = 1
                        else:
                            total_coverage_link_dict[(segs[i], segs[i + 1])] += 1
                    total_reads_length += read_length

    for key in total_coverage_dict:
        index_dict[key][0]["seg"] += total_coverage_dict[key]
    for key in total_coverage_link_dict:
        index_dict[key[0]][0][key[1]] += total_coverage_link_dict[key]
        index_dict[key[1]][0][key[0]] += total_coverage_link_dict[key]

    return total_breakpoint_list_del, total_breakpoint_list_ins, total_reads_length/3.05e9
"""



######Clustering breakpoints into bubbles
def process_breakpoint(breakpoint, gfaID2bubble_dict):
    '''
    Processing each breakpoint to get their bubble path.
    :param breakpoint:
    :param gfaID2bubble_dict:
    :return:
    '''
    ref_path = breakpoint[1]
    bubble_order = dict()
    segments = re.split(r"[<>]+", ref_path)[1:]
    for segment in segments:
        if len(gfaID2bubble_dict[segment]) == 1:
            index0 = gfaID2bubble_dict[segment][0]
            if index0 not in bubble_order:
                bubble_order[index0] = index0[1]
        elif len(gfaID2bubble_dict[segment]) == 2:
            index1, index2 = gfaID2bubble_dict[segment][0], gfaID2bubble_dict[segment][1]
            if index1[1] <= index2[1]:
                if (index1, index2) not in bubble_order:
                    bubble_order[(index1, index2)] = (index1[1] + index2[1])/2
            elif index1[1] > index2[1]:
                if (index2, index1) not in bubble_order:
                    bubble_order[(index2, index1)] = (index2[1] + index1[1]) / 2
    bubble_order = dict(sorted(bubble_order.items(), key=lambda item: item[1]))

    bubble_path = tuple(bubble_order.keys())
    breakpoint_info = copy.deepcopy(breakpoint)

    return bubble_path, breakpoint_info

def breakpoint_cluster(breakpoint_list, share_gfaID2bubble_dict, thread):
    '''
    Clustering breakpoint into bubbles.
    Changing ref_path of each breakpoint to their "bubble path".
    Because for each segment, (1)if it is in a bubble, we can use the start of bubble to represent this bubble;
    (2)if it is the linkage between two bubbles, we can use a list[start of anterior bubble, start of posterior bubble] to represent this link.
    So a ref_path can be changed into a list containing representation described in (1) and (2)---"bubble path".
    Breakpoints indicating the same SV must have same/similar "bubble path";
    breakpoints having same/similar "bubble path" may indicate different SVs.
    :param breakpoint_list: from function read_gaf(), a nested list containing information of breakpoints.
    :param share_gfaID2bubble_dict: from "###calling bubbles", reverse dictionary of bubble_dict.
    :param thread: thread.
    :return: cluster_re: {bubble_path1:[breakpoint1, breakpoint3], bubble_path2:[breakpoint2], ...} (index is the index of breakpoint in breakpoint_list)
    '''
    cluster_re = {}

    batch_size = thread
    batch = []

    count = 0

    with ProcessPoolExecutor(max_workers=thread) as executor:
        while count + thread < len(breakpoint_list):
            batch = copy.deepcopy(breakpoint_list[count:count + thread])
            futures = [executor.submit(process_breakpoint, breakpoint_, share_gfaID2bubble_dict) for breakpoint_ in batch]
            for future in as_completed(futures):
                result = future.result()
                bubble_path, breakpoint_info = result[0], result[1]
                if bubble_path not in cluster_re:
                    cluster_re[bubble_path] = []
                cluster_re[bubble_path].append(breakpoint_info)
            batch = []
            count += thread

        batch = copy.deepcopy(breakpoint_list[count:])
        futures = [executor.submit(process_breakpoint, breakpoint_, share_gfaID2bubble_dict) for breakpoint_ in batch]
        for future in as_completed(futures):
            result = future.result()
            bubble_path, breakpoint_info = result[0], result[1]
            if bubble_path not in cluster_re:
                cluster_re[bubble_path] = []
            cluster_re[bubble_path].append(breakpoint_info)

    return cluster_re



###Decreasing the effect of repeat
def interval_combination(interval_list):
    if not interval_list:
        return []
    else:
        interval_list.sort(key = lambda x: x[0])
        merged_interval_list = []
        for interval in interval_list:
            if not merged_interval_list:
                merged_interval_list.append(interval)
            else:
                if interval[0] > merged_interval_list[-1][1]:
                    merged_interval_list.append(interval)
                else:
                    merged_interval_list[-1][1] = max(interval[1], merged_interval_list[-1][1])
        return merged_interval_list

def read_repeatmasker(repeatmasker_dir, fasta):
    reads_repeat_interval = {}
    if repeatmasker_dir != "":
        if repeatmasker_dir[len(repeatmasker_dir) - 7:] == ".fa.out":
            with open(repeatmasker_dir, "rt") as f_repeat:
                for eachline in f_repeat:
                    text = re.split(r"[ ]+", eachline.strip())
                    name = text[4]
                    start, end = int(text[5]), int(text[6])
                    if name not in reads_repeat_interval:
                        reads_repeat_interval[name] = []
                    reads_repeat_interval[name].append([start, end])
            for key in reads_repeat_interval:
                interval_list = reads_repeat_interval[key]
                merged_interval_list = interval_combination(interval_list)
                repeat_length = 0
                for interval in merged_interval_list:
                    repeat_length += interval[1] - interval[0]
                reads_repeat_interval[key] = repeat_length
        else:
            print("Failed to read repeatmasker result, result in unknown format, please input it in fa.out format.")

    key_set = set(reads_repeat_interval.keys())
    total_reads, percent095_reads = 0, 0
    with gzip.open(fasta, "rt") as fi:
        t = False
        query_name, query_name_re, query_length = "", "", 0
        for eachline in fi:
            text = eachline.strip()
            if text[0] == ">":
                total_reads += 1
                if query_name != "":
                    percent = reads_repeat_interval[query_name_re]/query_length
                    reads_repeat_interval[query_name] = percent
                    if percent >= 0.95:
                        percent095_reads += 1
                    query_name, query_name_re, query_length = "", "", 0
                    t = False
                qn = text[1:]
                qnr = qn.split()[0]
                if qnr in key_set:
                    t = True
                    query_name, query_name_re = qn, qnr
                else:
                    reads_repeat_interval[qnr] = 0
            else:
                if t == True:
                    query_length += len(text)

        if query_name != "":
            percent = reads_repeat_interval[query_name_re] / query_length
            reads_repeat_interval[query_name] = percent
            total_reads += 1
            if percent >= 0.95:
                percent095_reads += 1

    return reads_repeat_interval, total_reads, percent095_reads



######SV detection
def is_bubble_or_betweenbubble(a_tuple):
    if isinstance(a_tuple, tuple) and all(isinstance(i, tuple) for i in a_tuple):
        return "betweenbubble"
    else:
        return "bubble"

def function_weight(coverage, x):
    '''
    Calcualting weight of edge from its coverage.
    :param coverage: from index_gfa
    :param x: sequencing depth, user input
    :return: weight(float)
    '''
    return math.exp(x - coverage)

def calculate_path_coverage(path, index_gfa):
    reads_length, path_length = 0, 0
    for seg in path:
        seg_coverage = index_gfa[seg][0]["seg"]
        seg_length = index_gfa[seg][2]
        reads_length += seg_coverage * seg_length
        path_length += seg_length
    return reads_length/path_length, path_length

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

def interval_combination_alt(interval_list):
    if not interval_list:
        return [], []
    else:
        interval_list.sort(key = lambda x: x[0])
        merged_interval_list = []
        interval_cluster_list = []
        for interval in interval_list:
            if not merged_interval_list:
                merged_interval_list.append(interval)
                interval_cluster_list.append([interval])
            else:
                if interval[0] > merged_interval_list[-1][1]:
                    merged_interval_list.append(interval)
                    interval_cluster_list.append([interval])
                else:
                    merged_interval_list[-1][1] = max(interval[1], merged_interval_list[-1][1])
                    interval_cluster_list[-1].append(interval)
        return merged_interval_list, interval_cluster_list

def interval_cluster_for_ins(interval_list):
    if not interval_list:
        return [], []
    else:
        interval_list.sort(key=lambda x: x[0])
        merged_interval_list = []
        interval_cluster_list = []
        for interval in interval_list:
            if not merged_interval_list:
                merged_interval_list.append(interval)
                interval_cluster_list.append([interval])
            else:
                if interval[0] - merged_interval_list[-1][1] > 50:
                    merged_interval_list.append(interval)
                    interval_cluster_list.append([interval])
                else:
                    merged_interval_list[-1][1] = max(interval[1], merged_interval_list[-1][1])
                    interval_cluster_list[-1].append(interval)
        return merged_interval_list, interval_cluster_list

def statistical_hypothesis_test(path_coverage, support_reads, misalignment_rate):
    '''

    :param path_coverage:
    :param support_reads:
    :param misalignment_rate:
    :return:
    '''
    p_value = 1 - binom.cdf(support_reads - 1, path_coverage, misalignment_rate)
    return p_value

def if_intersection_exist(interval1, interval2):
    max_0 = max(interval1[0], interval2[0])
    min_1 = min(interval1[1], interval2[1])
    if max_0 < min_1:
        return True
    else:
        return False

def from_interval_list_to_SV_interval_on_path_del(interval_list, graph, path, seg2base_dict, ref_seg2chr, chr2ref_seg):
    '''
    Information needed to be describe in vcf:
    chrom, pos, [id], ref, alt, [qual], [filter], info(on_ref/not_on_ref, if not in ref` give its path and postion on path, ins/del), GT, samples
    :param interval_list:
    :param graph:
    :param path:
    :param seg2base_dict:
    :param ref_seg2chr:
    :param chr2ref_seg:
    :return:
    '''
    start, end, support_reads_num = 0, 0, 0
    for interval in interval_list:
        start += interval[0]
        end += interval[1]
        support_reads_num += 1
    start = round(start/support_reads_num)
    end = round(end/support_reads_num)
    start_ref_seg, end_ref_seg = -1, -1
    i = 0
    for i in range(len(path)):
        seg = path[i]
        if seg in ref_seg2chr:
            start_ref_seg = i
        if if_intersection_exist(seg2base_dict[seg], (start, end)):
            break
    j = 0
    for j in range(len(path)):
        seg = path[len(path) - i - 1]
        if seg in ref_seg2chr:
            end_ref_seg = j
        if if_intersection_exist(seg2base_dict[seg], (start, end)):
            break
    path_seg = path[start_ref_seg : end_ref_seg + 1]

    seg1, strand1, pos1, length1, seg2, strand2, pos2, length2 = "", "", -1, -1, "", "", -1, -1      #For graph update!!!
    for seg in seg2base_dict:
        if start > seg2base_dict[seg][0] and start <= seg2base_dict[seg][1]:
            seg1, strand1 = seg, graph.nodes[seg]["strand"]
            pos1, length1 = start - seg2base_dict[seg][0], seg2base_dict[seg][1] - seg2base_dict[seg][0]
        if end >= seg2base_dict[seg][0] and end < seg2base_dict[seg][1]:
            seg2, strand2 = seg, graph.nodes[seg]["strand"]
            pos2, length2 = end - seg2base_dict[seg][0], seg2base_dict[seg][1] - seg2base_dict[seg][0]

    path_ = ""      #
    path_base = ""
    chrom = ref_seg2chr[path_seg[0]]      #

    path_length = 0
    for seg in path_seg:
        if graph.nodes[seg]["strand"] == "+":
            path_ += ">"
            path_ += seg
            path_base += graph.nodes[seg]["base"]
        else:
            path_ += "<"
            path_ += seg
            path_base += reverse_bases(graph.nodes[seg]["base"])
        path_length += len(graph.nodes[seg]["base"])
    length_before_start_ref_seg = 0
    k = 0
    for k in range(0, start_ref_seg):
        seg = path[k]
        length_before_start_ref_seg += len(graph.nodes[seg]["base"])
    length_before_end_ref_seg = length_before_start_ref_seg
    m = 0
    for m in range(start_ref_seg, end_ref_seg):
        seg = path[m]
        length_before_end_ref_seg += len(graph.nodes[seg]["base"])

    path_del_start_index, path_del_end_index = start - length_before_start_ref_seg, end - length_before_start_ref_seg      #

    on_or_not_on_ref = ""      #

    tuple_ = tuple(chr2ref_seg[chrom].keys())
    ref_path_seg = tuple_[tuple_.index(path_seg[0]) : tuple_.index(path_seg[-1]) + 1]
    ref_path_base = ""
    for seg in ref_path_seg:
        if graph.nodes[seg]["strand"] == "+":
            ref_path_base += graph.nodes[seg]["base"]
        else:
            ref_path_base += reverse_bases(graph.nodes[seg]["base"])
    if ref_path_seg == path_seg:
        on_or_not_on_ref = "OnRef"
    else:
        on_or_not_on_ref = "NotOnRef"

    pos, ref, alt = -1, "", ""      #

    if on_or_not_on_ref == "OnRef":
        pos = chr2ref_seg[chrom][path_seg[0]] + path_del_start_index
        ref = ref_path_base[path_del_start_index - 1 : path_del_end_index]
        alt = (ref_path_base[path_del_start_index - 1], "")
        return chrom, pos, "", ref, alt, 255, ".", ("OnRef", path_, path_del_start_index, path_del_end_index, "del"), "GT", "1/.", (seg1, strand1, pos1, length1, seg2, strand2, pos2, length2)
    elif on_or_not_on_ref == "NotOnRef":
        a = chr2ref_seg[chrom][tuple_[tuple_.index(path_seg[0]) + 1]] - chr2ref_seg[chrom][path_seg[0]]
        b = chr2ref_seg[chrom][tuple_[tuple_.index(path_seg[-1]) + 1]] - chr2ref_seg[chrom][path_seg[-1]]
        if path_del_start_index >= a:
            pos = chr2ref_seg[chrom][tuple_[tuple_.index(path_seg[0]) + 1]]
            if len(path_base) - path_del_end_index >= b:
                ref = ref_path_base[a - 1 : len(ref_path_base) - b]
                alt = (path_base[a - 1 : path_del_start_index], path_base[path_del_end_index : len(path_base) - b])
            else:
                ref = ref_path_base[a - 1 : len(ref_path_base) - (len(path_base) - path_del_end_index)]
                alt = (path_base[a - 1 : path_del_start_index], "")
        else:
            pos = chr2ref_seg[chrom][path_seg[0]] + path_del_start_index
            if len(path_base) - path_del_end_index >= b:
                ref = ref_path_base[path_del_start_index - 1 : len(ref_path_base) - b]
                alt = ("", path_base[path_del_end_index : len(path_base) - b])
            else:
                ref = ref_path_base[path_del_start_index - 1 : len(ref_path_base) - (len(path_base) - path_del_end_index)]
                alt = (ref[0], "")
        return chrom, pos, "", ref, alt, 255, ".", ("NotOnRef", path_, path_del_start_index, path_del_end_index, "del"), "GT", "1/.", (seg1, strand1, pos1, length1, seg2, strand2, pos2, length2)

def from_breakpoints_list_to_SV_interval_on_path_ins(breakpoints_list, graph, fasta, path, path_bases,
                                                     path_seg2base_dict, ref_seg2chr, chr2ref_seg,
                                                     breakpoint_location_dict, temporary_directory, sequencing_tech):
    dict_ = dict()
    total_count = 0
    for breakpoint in breakpoints_list:
        readsid = breakpoint[0]
        start_on_query, end_on_query, strand = breakpoint[5][0], breakpoint[5][1], breakpoint[5][2]
        path_reverse = breakpoint_location_dict[breakpoint][3]
        a = None
        if strand == "+":
            if path_reverse == False:
                a = False
            elif path_reverse == True:
                a = True
        elif strand == "-":
            if path_reverse == False:
                a = True
            elif path_reverse == True:
                a = False
        if readsid not in dict_:
            dict_[readsid] = list()
        dict_[readsid].append([start_on_query, end_on_query, a])
        total_count += 1

    reads_part = dict()
    for read in dict_:
        bases = fasta[read]
        for eachone in dict_[read]:
            soq, eoq, ifreverse = eachone[0], eachone[1], eachone[2]
            if soq - 1000 >= 0:
                soq = soq - 1000
            else:
                soq = 0
            if eoq + 1000 <= len(bases):
                eoq = eoq + 1000
            else:
                eoq = len(bases)
            if read not in reads_part:
                reads_part[read] = list()
            if ifreverse == False:
                reads_part[read].append(bases[soq: eoq])
            elif ifreverse == True:
                reads_part[read].append(reverse_bases(bases[soq: eoq]))
    '''
    with gzip.open(fasta, "rt") as fi:
        t = False
        reads = ""
        bases = ""
        for eachline in fi:
            text = eachline.strip()
            if text[0] == ">":
                name = text[1:]
                if t == True:
                    for eachone in dict_[reads]:
                        soq, eoq, ifreverse = eachone[0], eachone[1], eachone[2]
                        if soq - 1000 >= 0:
                            soq = soq - 1000
                        else:
                            soq = 0
                        if eoq + 1000 <= len(bases):
                            eoq = eoq + 1000
                        else:
                            eoq = len(bases)
                        if reads not in reads_part:
                            reads_part[reads] = list()
                        if ifreverse == False:
                            reads_part[reads].append(bases[soq : eoq])
                        elif ifreverse == True:
                            reads_part[reads].append(reverse_bases(bases[soq : eoq]))
                        c += 1
                    t, reads, bases = False, "", ""
                if c == total_count:
                    break
                if name in dict_:
                    t = True
                    reads = name
            else:
                if t == True:
                    bases += text
    '''
    unique_tmp_subdir = tempfile.mkdtemp(dir = temporary_directory)
    queryfa, cnsfa, reffa = "", os.path.join(unique_tmp_subdir, "output_prefix.ctg.fa"), ""
    with tempfile.NamedTemporaryFile(mode = "w", delete = False, dir = unique_tmp_subdir, suffix = ".fa") as query_fa:
        for reads in reads_part:
            co = 1
            for each in reads_part[reads]:
                query_fa.write(">{}_{}\n".format(reads, co))
                query_fa.write("{}\n".format(each))
                co += 1
        queryfa = query_fa.name
    with tempfile.NamedTemporaryFile(mode = "w", delete = False, dir = unique_tmp_subdir, suffix = ".fa") as ref_fa:
        ref_fa.write(">path\n{}\n".format(path_bases))
        reffa = ref_fa.name
    subprocess.run([
        "wtdbg2", "-x", sequencing_tech, "-g", "10k", "-e", "2", "-L", "500",
        "-i", queryfa, "-fo", "prefix"
    ], check = True, cwd = unique_tmp_subdir)
    subprocess.run([
        "wtpoa-cns", "-i", "prefix.ctg.lay.gz", "-fo", cnsfa
    ], check = True, cwd = unique_tmp_subdir)
    result = subprocess.run([
        "minimap2", "-ax", "asm5", reffa, cnsfa
    ], check = True, cwd = unique_tmp_subdir, stderr = subprocess.PIPE, stdout = subprocess.PIPE)
    aln_sam = result.stdout.decode("utf-8")
    start_on_ref_path, end_on_ref_path, start_on_query, end_on_query = -1, -1, -1, -1
    for line in aln_sam.splitlines():
        if line.startswith("@"):
            continue
        else:
            text_split = line.strip().split("\t")
            pos, cigar = int(text_split[3]), text_split[5]
            cigar_num_list = re.split(r"[MIDNSHPX=]+", cigar)
            del cigar_num_list[-1]
            cigar_op_list = re.split(r"[1234567890]+", cigar)
            del cigar_op_list[0]

            possible_sv_ins = []
            k = 0
            for k in range(len(cigar_op_list)):
                if cigar_op_list[k] == "I":
                    if int(cigar_num_list[k]) >= 10:
                        possible_sv_ins.append(k)
            i = 0
            while i <= len(possible_sv_ins) - 1:
                k_start = possible_sv_ins[i]
                length_ins = int(cigar_num_list[possible_sv_ins[i]])
                while (possible_sv_ins[i] + 2 == possible_sv_ins[i + 1]) and (
                        int(cigar_num_list[possible_sv_ins[i] + 1]) == 1):
                    i += 1
                    length_ins += int(cigar_num_list[possible_sv_ins[i]])
                k_end = possible_sv_ins[i]
                if length_ins < 30:
                    pass
                elif length_ins >= 30:
                    ref_length1 = count_query_and_ref_length_from_cigar(cigar_op_list, cigar_num_list, 0, k_start,
                                                                        "ref")
                    ref_length2 = count_query_and_ref_length_from_cigar(cigar_op_list, cigar_num_list, k_start,
                                                                        k_end + 1, "ref")
                    query_length1 = count_query_and_ref_length_from_cigar(cigar_op_list, cigar_num_list, 0, k_start,
                                                                          "query")
                    query_length2 = count_query_and_ref_length_from_cigar(cigar_op_list, cigar_num_list, k_start,
                                                                          k_end + 1, "query")

                    start_on_ref_path, end_on_ref_path = pos + ref_length1 - 1, pos + ref_length1 + ref_length2 -1
                    start_on_query, end_on_query = query_length1, query_length1 + query_length2
    if os.path.exists(unique_tmp_subdir):
        shutil.rmtree(unique_tmp_subdir)

    if start_on_ref_path == -1:
        return "", -1, "", "", "", -1, "", (), "", "", ()
    else:
        seg1, strand1, pos1, length1, seg2, strand2, pos2, length2 = "", "", -1, -1, "", "", -1, -1
        for seg in path_seg2base_dict:
            if start_on_ref_path > path_seg2base_dict[seg][0] and start_on_ref_path <= path_seg2base_dict[seg][1]:
                seg1, strand1 = seg, graph.nodes[seg]["strand"]
                pos1, length1 = start_on_ref_path - path_seg2base_dict[seg][0], path_seg2base_dict[seg][1] - path_seg2base_dict[seg][0]
            if end_on_ref_path >= path_seg2base_dict[seg][0] and end_on_ref_path < path_seg2base_dict[seg][1]:
                seg2, strand = seg, graph.nodes[seg]["strand"]
                pos2, length2 = end_on_ref_path - path_seg2base_dict[seg][0], path_seg2base_dict[seg][1] - path_seg2base_dict[seg][0]
        alter = ""      #
        with open(os.path.join(unique_tmp_subdir, cnsfa), "rt") as f_cns:
            cns = ""
            for eachline in f_cns:
                text = eachline.strip()
                if text[0] != ">":
                    cns += text
            alter = cns[start_on_query : end_on_query]
        start_ref_seg, end_ref_seg = -1, -1
        i = 0
        for i in range(len(path)):
            seg = path[i]
            if seg in ref_seg2chr:
                start_ref_seg = i
            if if_intersection_exist(path_seg2base_dict[seg], (start_on_ref_path, end_on_ref_path)):
                break
        j = 0
        for j in range(len(path)):
            seg = path[len(path) - i - 1]
            if seg in ref_seg2chr:
                end_ref_seg = j
            if if_intersection_exist(path_seg2base_dict[seg], (start_on_ref_path, end_on_ref_path)):
                break
        path_seg = path[start_ref_seg : end_ref_seg + 1]
        path_ = ""  #
        base_ = ""
        chrom = ref_seg2chr[path_seg[0]]  #

        path_length = 0
        for seg in path_seg:
            if graph.nodes[seg]["strand"] == "+":
                path_ += ">"
                path_ += seg
                base_ += graph.nodes[seg]["base"]
            else:
                path_ += "<"
                path_ += seg
                base_ += reverse_bases(graph.nodes[seg]["base"])
            path_length += len(graph.nodes[seg]["base"])
        length_before_start_ref_seg = 0
        k = 0
        for k in range(0, start_ref_seg):
            seg = path[k]
            length_before_start_ref_seg += len(graph.nodes[seg]["base"])
        length_before_end_ref_seg = length_before_start_ref_seg
        m = 0
        for m in range(start_ref_seg, end_ref_seg):
            seg = path[m]
            length_before_end_ref_seg += len(graph.nodes[seg]["base"])
        path_ins_start_index, path_ins_end_index = start_on_ref_path - length_before_start_ref_seg, end_on_ref_path - length_before_start_ref_seg  #
        on_or_not_on_ref = ""  #
        tuple_ = tuple(chr2ref_seg[chrom].keys())
        ref_path_seg = tuple_[tuple_.index(path_seg[0]): tuple_.index(path_seg[-1]) + 1]
        ref_path_base = ""
        for seg in ref_path_seg:
            if graph.nodes[seg]["strand"] == "+":
                ref_path_base += graph.nodes[seg]["base"]
            else:
                ref_path_base += reverse_bases(graph.nodes[seg]["base"])
        if ref_path_seg == path_seg:
            on_or_not_on_ref = "OnRef"
        else:
            on_or_not_on_ref = "NotOnRef"

        pos, ref = -1, ""  #
        if on_or_not_on_ref == "OnRef":
            pos = chr2ref_seg[chrom][path_seg[0]] + path_ins_start_index
            ref = ref_path_base[path_ins_start_index - 1: path_ins_end_index]
            alt = (ref_path_base[path_ins_start_index - 1], alter, "")
            return chrom, pos, "", ref, alt, 255, ".", ("OnRef", path_, path_ins_start_index, path_ins_end_index, "ins"), "GT", "1/.", (seg1, strand1, pos1, length1, seg2, strand2, pos2, length2)
        elif on_or_not_on_ref == "NotOnRef":
            a = chr2ref_seg[chrom][tuple_[tuple_.index(path_seg[0]) + 1]] - chr2ref_seg[chrom][path_seg[0]]
            b = chr2ref_seg[chrom][tuple_[tuple_.index(path_seg[-1]) + 1]] - chr2ref_seg[chrom][path_seg[-1]]
            if path_ins_start_index >= a:
                pos = chr2ref_seg[chrom][tuple_[tuple_.index(path_seg[0]) + 1]]
                if len(base_) - path_ins_end_index >= b:
                    ref = ref_path_base[a - 1: len(ref_path_base) - b]
                    alt = (base_[a - 1: path_ins_start_index], alter, base_[path_ins_end_index: len(base_) - b])
                else:
                    ref = ref_path_base[a - 1: len(ref_path_base) - (len(base_) - path_ins_end_index)]
                    alt = (base_[a - 1: path_ins_start_index], alter, "")
            else:
                pos = chr2ref_seg[chrom][path_seg[0]] + path_ins_start_index
                if len(base_) - path_ins_end_index >= b:
                    ref = ref_path_base[path_ins_start_index - 1: len(ref_path_base) - b]
                    alt = ("", alter, base_[path_ins_end_index: len(base_) - b])
                else:
                    ref = ref_path_base[path_ins_start_index - 1: len(ref_path_base) - (len(base_) - path_ins_end_index)]
                    alt = (ref[0], alter, "")
            return chrom, pos, "", ref, alt, 255, ".", ("NotOnRef", path_, path_ins_start_index, path_ins_end_index, "ins"), "GT", "1/.", (seg1, strand1, pos1, length1, seg2, strand2, pos2, length2)

def mapping_breakpoints_to_path(graph, path, breakpoint_list, reads_repeat_dict, ref_seg2chr, chr2ref_seg, path_coverage,
                                percent095_reads_percent, fasta, temporary_directory, sequencing_tech):
    '''
    Now path(s) in "bubble path"(genetic constitution) and breakpoints located in this "bubble path" are found.
    This function will map these breakpoints to the path(s), in order to find SVs.
    For insertion, fasta file will be read and local assembly will be processed to determine the inserted bases.

    If repeatmasker result is provided, also considering effect of repeat while detecting SV.
    :param graph:
    :param path:
    :param breakpoint_list:
    :param reads_repeat_dict:
    :param ref_seg2chr:
    :param chr2ref_seg:
    :param path_coverage: path_coverage != average_cover !!!
    :param percent095_reads_percent: from function read_repeatmasker()
    :return:
    '''
    bases = ""
    seg2base_dict = {}
    breakpoint_location = {}
    index = 0
    for seg in path:
        if graph.nodes[seg]["strand"] == "+":
            bases += graph.nodes[seg]["base"]
        elif graph.nodes[seg]["strand"] == "-":
            bases += reverse_bases(graph.nodes[seg]["base"])
        seg2base_dict[seg] = [index, len(bases)]      #[start on path_bases, end on paths bases] for each segment
        index = len(bases)
    seg2base_dict = dict(sorted(seg2base_dict.items(), key = lambda item: item[1][0]))
    mode = None
    if breakpoint_list[0][4] == "deletion":
        mode = "d"
    elif breakpoint_list[0][4] == "insertion":
        mode = "i"
    if mode == "d":
        #[readID, ref_path, start_on_ref_path, end_on_ref_path, "deletion", None]
        breakpoint_location = dict()
        for breakpoint in breakpoint_list:
            ref_path, start_on_ref_path, end_on_ref_path = breakpoint[1], breakpoint[2], breakpoint[3]
            ref_path_seg, ref_path_strand = re.split(r"[<>]+", ref_path)[1:], re.findall(r"[<>]+", ref_path)
            length_ref_path_seg = len(ref_path_seg)
            length_ref_path = 0
            reverse = False
            start, end = "", ""
            start_on_bases, end_on_bases = -1, -1
            m = 0
            for m in range(length_ref_path_seg):
                length_ref_path += len(graph.nodes[ref_path_seg[m]]["base"])
            i = 0
            for i in range(length_ref_path_seg):
                s, st = ref_path_seg[i], ref_path_strand[i]
                if s in path:
                    start = s
                    if st == graph.nodes[s]["strand"]:
                        pass
                    elif st != graph.nodes[s]["strand"]:
                        reverse = True
                    break
            j = 0
            for j in range(length_ref_path_seg):
                s, st = ref_path_seg[length_ref_path_seg - 1 - j], ref_path_strand[length_ref_path_seg - 1 - j]
                if s in path:
                    end = s
                    break
            if start != "" and end != "":
                if reverse == False:
                    if start == ref_path_seg[0]:
                        start_on_bases = seg2base_dict[start][0] + start_on_ref_path
                    else:
                        start_on_bases = seg2base_dict[start][0]
                    if end == ref_path_seg[-1]:
                        end_on_bases = seg2base_dict[end][1] - (length_ref_path - end_on_ref_path)
                    else:
                        end_on_bases = seg2base_dict[end][1]
                elif reverse == True:
                    if start == ref_path_seg[-1]:
                        end_on_bases = seg2base_dict[start][1] - start_on_ref_path
                    else:
                        end_on_bases = seg2base_dict[start][1]
                    if end == ref_path_seg[0]:
                        start_on_bases = seg2base_dict[end][0] + (length_ref_path - end_on_ref_path)
                    else:
                        start_on_bases = seg2base_dict[end][0]
                breakpoint_location[tuple(breakpoint)] = [start_on_bases, end_on_bases] + [1, reverse]

        if isinstance(reads_repeat_dict, dict):
            for breakpoint in breakpoint_location:
                read_id = breakpoint[0]
                if read_id in reads_repeat_dict:
                    if reads_repeat_dict[read_id] >= 0.95:
                        breakpoint_location[breakpoint][2] = 0.85

        locations = []
        for breakpoint in breakpoint_location:
            locations.append([breakpoint_location[breakpoint][0], breakpoint_location[breakpoint][1]])
        merged_interval_list, interval_cluster_list = interval_combination_alt(locations)

        deletionSV = []
        if merged_interval_list == []:
            return deletionSV
        else:
            i = 0
            for i in range(len(merged_interval_list)):
                low_repeat_reads_support, high_repeat_reads_support = 0, 0
                interval_cluster = interval_cluster_list[i]
                for breakpoint in breakpoint_location:
                    if [breakpoint_location[breakpoint][0], breakpoint_location[breakpoint][1]] in interval_cluster and breakpoint_location[breakpoint][2] == 1:
                        low_repeat_reads_support += 1
                    elif [breakpoint_location[breakpoint][0], breakpoint_location[breakpoint][1]] in interval_cluster and breakpoint_location[breakpoint][2] == 0.85:
                        high_repeat_reads_support += 1

                #chrom, pos, "", ref, alt, 255, ".", ["OnRef", "", "", "", "del"], "GT", "1/."
                #chrom, pos, "", ref, alt, 255, ".", ["NotOnRef", path_, path_del_start_index, path_del_end_index, "del"], "GT", "1/."
                if low_repeat_reads_support + high_repeat_reads_support > path_coverage:
                    chrom, pos, name, ref, alt, qual, filter, list1, format, sample, new_seg_link =\
                        from_interval_list_to_SV_interval_on_path_del(interval_cluster, graph, path, seg2base_dict, ref_seg2chr, chr2ref_seg)
                    alt = alt[0] + alt[1]
                    deletionSV.append((chrom, pos, name, ref, alt, qual, filter, list1, format, sample, new_seg_link))
                else:
                    part = (path_coverage - low_repeat_reads_support - high_repeat_reads_support) * (1 - percent095_reads_percent)
                    if int(part) == part:
                        part = int(part)
                    elif int(part) != part:
                        part = int(part) + 1
                    low_repeat_reads_total = part + low_repeat_reads_support
                    high_repeat_reads_total = path_coverage - low_repeat_reads_total
                    p_value = statistical_hypothesis_test(low_repeat_reads_total, low_repeat_reads_support, 0.065) * statistical_hypothesis_test(high_repeat_reads_total, high_repeat_reads_support, 0.15)
                    if p_value <= 0.05:
                        chrom, pos, name, ref, alt, qual, filter, list1, format, sample, new_seg_link =\
                            from_interval_list_to_SV_interval_on_path_del(interval_cluster, graph, path, seg2base_dict, ref_seg2chr, chr2ref_seg)
                        alt = alt[0] + alt[1]
                        deletionSV.append((chrom, pos, name, ref, alt, qual, filter, list1, format, sample, new_seg_link))
            return deletionSV

    elif mode == "i":
        print("start mapping")
        #[readID, ref_path, start_on_ref_path, end_on_ref_path, "insertion", (start_on_query, end_on_query, strand)]
        breakpoint_location = dict()
        for breakpoint in breakpoint_list:
            ref_path, start_on_ref_path, end_on_ref_path = breakpoint[1], breakpoint[2], breakpoint[3]
            ref_path_seg, ref_path_strand = re.split(r"[<>]+", ref_path)[1:], re.findall(r"[<>]+", ref_path)
            length_ref_path_seg = len(ref_path_seg)
            length_ref_path = 0
            reverse = False
            start, end = "", ""
            start_on_bases, end_on_bases = -1, -1
            m = 0
            for m in range(length_ref_path_seg):
                length_ref_path += len(graph.nodes[ref_path_seg[m]]["base"])
            i = 0
            for i in range(length_ref_path_seg):
                s, st = ref_path_seg[i], ref_path_strand[i]
                if s in path:
                    start = s
                    if st == graph.nodes[s]["strand"]:
                        pass
                    elif st != graph.nodes[s]["strand"]:
                        reverse = True
                    break
            j = 0
            for j in range(length_ref_path_seg):
                s, st = ref_path_seg[length_ref_path_seg - 1 - j], ref_path_strand[length_ref_path_seg - 1 - j]
                if s in path:
                    end = s
                    break
            if start != "" and end != "":
                if reverse == False:
                    if start == ref_path_seg[0]:
                        start_on_bases = seg2base_dict[start][0] + start_on_ref_path
                    else:
                        start_on_bases = seg2base_dict[start][0]
                    if end == ref_path_seg[-1]:
                        end_on_bases = seg2base_dict[end][1] - (length_ref_path - end_on_ref_path)
                    else:
                        end_on_bases = seg2base_dict[end][1]
                elif reverse == True:
                    if start == ref_path_seg[-1]:
                        end_on_bases = seg2base_dict[start][1] - start_on_ref_path
                    else:
                        end_on_bases = seg2base_dict[start][1]
                    if end == ref_path_seg[0]:
                        start_on_bases = seg2base_dict[end][0] + (length_ref_path - end_on_ref_path)
                    else:
                        start_on_bases = seg2base_dict[end][0]
                breakpoint_location[tuple(breakpoint)] = [start_on_bases, end_on_bases] + [1, reverse]

        if isinstance(reads_repeat_dict, dict):
            for breakpoint in breakpoint_location:
                read_id = breakpoint[0]
                if read_id in reads_repeat_dict:
                    if reads_repeat_dict[read_id] >= 0.95:
                        breakpoint_location[breakpoint][2] = 0.85

        locations = []
        for breakpoint in breakpoint_location:
            locations.append([breakpoint_location[breakpoint][0], breakpoint_location[breakpoint][1]])
        merged_interval_list, interval_cluster_list = interval_cluster_for_ins(locations)

        print("end mapping, start sv detecting")

        deletionSV = []
        if merged_interval_list == []:
            print("No breakpoints support")
            return deletionSV
        else:
            i = 0
            for i in range(len(merged_interval_list)):
                breakpoints_in_this_merged_interval = []      #
                low_repeat_reads_support, high_repeat_reads_support = 0, 0
                interval_cluster = interval_cluster_list[i]
                merged_interval = merged_interval_list[i]
                for breakpoint in breakpoint_location:
                    if [breakpoint_location[breakpoint][0], breakpoint_location[breakpoint][1]] in interval_cluster and \
                            breakpoint_location[breakpoint][2] == 1:
                        breakpoints_in_this_merged_interval.append(breakpoint)
                        low_repeat_reads_support += 1
                    elif [breakpoint_location[breakpoint][0], breakpoint_location[breakpoint][1]] in interval_cluster and \
                            breakpoint_location[breakpoint][2] == 0.85:
                        breakpoints_in_this_merged_interval.append(breakpoint)
                        high_repeat_reads_support += 1

                if low_repeat_reads_support + high_repeat_reads_support > path_coverage:
                    chrom, pos, name, ref, alt, qual, filter, list1, format, sample, new_seg_link = \
                        from_breakpoints_list_to_SV_interval_on_path_ins(breakpoints_in_this_merged_interval, graph, fasta,
                                                                         path, bases, seg2base_dict, ref_seg2chr, chr2ref_seg,
                                                                         breakpoint_location, temporary_directory, sequencing_tech)
                    if chrom != "":
                        alt = alt[0] + alt[1] + alt[2]
                        deletionSV.append((chrom, pos, name, ref, alt, qual, filter, list1, format, sample, new_seg_link))
                else:
                    part = (path_coverage - low_repeat_reads_support - high_repeat_reads_support) * (1 - percent095_reads_percent)

                    if int(part) == part:
                        part = int(part)
                    elif int(part) != part:
                        part = int(part) + 1
                    low_repeat_reads_total = part + low_repeat_reads_support
                    high_repeat_reads_total = path_coverage - low_repeat_reads_total
                    p_value = statistical_hypothesis_test(low_repeat_reads_total, low_repeat_reads_support, 0.065) * statistical_hypothesis_test(high_repeat_reads_total, high_repeat_reads_support, 0.15)
                    if p_value <= 0.05:
                        chrom, pos, name, ref, alt, qual, filter, list1, format, sample, new_seg_link = \
                            from_breakpoints_list_to_SV_interval_on_path_ins(breakpoints_in_this_merged_interval, graph,
                                                                             fasta, path, bases, seg2base_dict,
                                                                             ref_seg2chr, chr2ref_seg, breakpoint_location,
                                                                             temporary_directory, sequencing_tech)
                        if chrom != "":
                            alt = alt[0] + alt[1] + alt[2]
                            deletionSV.append((chrom, pos, name, ref, alt, qual, filter, list1, format, sample, new_seg_link))
                    else:
                        print("p value = {}".format(p_value))
            return deletionSV

def process_cluster(bubble_path, breakpoint_list, bubble2gfaID, gfa_file, index_gfa, offset_list, ref_seg2chr,
                    chr2ref_seg, average_cover, percent095_reads_percent, fasta, temporary_directory,
                    sequencing_tech, reads_repeat_dict = None):
    '''
    Picking sub-graph of each path, and detecting SV from previous information.
    For SV, calling its path in this subgraph, considering heterozygosis.
    If repeatmasker result is provided, also considering effect of repeat while detecting SV.
    一条path上可能存在多个SV
    :param bubble_path: key in dict "bubble_path_clustered"
    :param breakpoint_list: breakpoints belong to "bubble_path_clustered"
    :param bubble2gfaID
    :param gfa_file: pan-graph in gfa format, read by module mmap
    :param index_gfa: index of gfa
    :param offset_list: offset of eachline in gfa file
    :param fasta: fasta file so that this script can get the ALT of an insertion variant.
    :param reads_repeat_dict: if a repeatmasker out is provided,
    :return:
    '''
    #node in this "bubble path"
    print(bubble_path)
    time0 = time.time()
    segments, ref_segments = set(), set()
    start_seg, end_seg = "", ""
    i = 0

    time1 = time.time()
    for i in range(len(bubble_path)):
        element = bubble_path[i]
        if i == 0:
            print(element)
            timea = time.time()
            if isinstance(element, tuple) and all(isinstance(j, tuple) for j in element):
                list1, list2 = bubble2gfaID[element[0][0]][element[0][1]][1], \
                    bubble2gfaID[element[1][0]][element[1][1]][1]
                if list1[-1] == list2[0]:
                    start_seg = list2[0]
                    segments.add(start_seg)
                    if start_seg in ref_seg2chr:
                        ref_segments.add(start_seg)
                elif list1[0] == list2[-1]:
                    start_seg = list1[0]
                    segments.add(start_seg)
                    if start_seg in ref_seg2chr:
                        ref_segments.add(start_seg)
            else:
                list0 = bubble2gfaID[element[0]][element[1]][1]
                start_seg = list0[0]
                for seg in list0:
                    segments.add(seg)
                    if seg in ref_seg2chr:
                        ref_segments.add(seg)
            print("{}s".format(time.time() - timea))
        if i == len(bubble_path) - 1:
            print(element)
            timea = time.time()
            if isinstance(element, tuple) and all(isinstance(j, tuple) for j in element):
                list1, list2 = bubble2gfaID[element[0][0]][element[0][1]][1], \
                    bubble2gfaID[element[1][0]][element[1][1]][1]
                if list1[-1] == list2[0]:
                    end_seg = list1[-1]
                    segments.add(end_seg)
                    if end_seg in ref_seg2chr:
                        ref_segments.add(end_seg)
                elif list1[0] == list2[-1]:
                    end_seg = list2[-1]
                    segments.add(end_seg)
                    if end_seg in ref_seg2chr:
                        ref_segments.add(end_seg)
            else:
                list0 = bubble2gfaID[element[0]][element[1]][1]
                end_seg = list0[-1]
                for seg in list0:
                    segments.add(seg)
                    if seg in ref_seg2chr:
                        ref_segments.add(seg)
            print("{}s".format(time.time() - timea))
        if i > 0 and i < len(bubble_path) - 1:
            print(element)
            timea = time.time()
            if isinstance(element, tuple) and all(isinstance(j, tuple) for j in element):
                pass
            else:
                list0 = bubble2gfaID[element[0]][element[1]][1]
                for seg in list0:
                    segments.add(seg)
                    if seg in ref_seg2chr:
                        ref_segments.add(seg)
            print("{}s".format(time.time() - timea))

    print("start graph construction.")
    #graph construction
    segments_remained = copy.deepcopy(segments)      #set of segments whose edges haven't been added to graph
    segments_remained_2 = copy.deepcopy(segments)      #set of segments whose strand haven't been determined
    G = nx.DiGraph()
    for seg in segments:
        G.add_node(seg,  strand = "", base = "")
    with open(gfa_file, "rb") as fi:
        mm = mmap.mmap(fi.fileno(), 0, access=mmap.ACCESS_READ)

        now_segments = []
        for ss in ref_segments:
            G.nodes[ss]["strand"] = "+"
            segments_remained_2.discard(ss)
            segments_remained.discard(ss)
            seg_text_start = mm[offset_list[index_gfa[ss][1][0]] : offset_list[index_gfa[ss][1][0] + 1]].decode("utf-8")
            seg_text_start_split = seg_text_start.strip().split("\t")
            G.nodes[ss]["base"] = seg_text_start_split[2]
            link_index_start = index_gfa[ss][1][1:]
            for line in link_index_start:
                text = mm[offset_list[line]: offset_list[line + 1]].decode("utf-8")
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
                                G.add_edge(gfaID1, gfaID2, weight = function_weight(index_gfa[gfaID1][0][gfaID2], average_cover))
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
                                G.add_edge(gfaID2, gfaID1, weight=function_weight(index_gfa[gfaID1][0][gfaID2], average_cover))
                elif gfaID2 == ss:
                    if strand2 == "+":
                        if gfaID1 in segments:
                            if gfaID1 in segments_remained_2:
                                G.nodes[gfaID1]["strand"] = strand1
                                segments_remained_2.discard(gfaID1)
                                now_segments.append(gfaID1)
                            if not G.has_edge(gfaID1, gfaID2):
                                G.add_edge(gfaID1, gfaID2, weight = function_weight(index_gfa[gfaID1][0][gfaID2], average_cover))
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
                                G.add_edge(gfaID2, gfaID1, weight = function_weight(index_gfa[gfaID1][0][gfaID2], average_cover))
        while segments_remained != set():
            now_segments_ts = []
            for seg in now_segments:
                segments_remained.discard(seg)
                seg_text = mm[offset_list[index_gfa[seg][1][0]] : offset_list[index_gfa[seg][1][0] + 1]].decode("utf-8")
                seg_text_split = seg_text.strip().split("\t")
                G.nodes[seg]["base"] = seg_text_split[2]
                link_index = index_gfa[seg][1][1:]
                for line in link_index:
                    text = mm[offset_list[line]: offset_list[line + 1]].decode("utf-8")
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
                                    G.add_edge(gfaID1, gfaID2, weight = function_weight(index_gfa[gfaID1][0][gfaID2], average_cover))
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
                                    G.add_edge(gfaID2, gfaID1, weight = function_weight(index_gfa[gfaID1][0][gfaID2], average_cover))
                    elif gfaID2 == seg:
                        if strand2 == G.nodes[gfaID2]["strand"]:
                            if gfaID1 in segments:
                                if gfaID1 in segments_remained_2:
                                    G.nodes[gfaID1]["strand"] = strand1
                                    segments_remained_2.discard(gfaID1)
                                    now_segments_ts.append(gfaID1)
                                if not G.has_edge(gfaID1, gfaID2):
                                    G.add_edge(gfaID1, gfaID2, weight = function_weight(index_gfa[gfaID1][0][gfaID2], average_cover))
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
                                    G.add_edge(gfaID2, gfaID1, weight = function_weight(index_gfa[gfaID1][0][gfaID2], average_cover))
            now_segments = copy.deepcopy(now_segments_ts)

        mm.close()
    print("end graph construction. {}s".format(time.time() - time1))
    #path determination
    timep = time.time()

    SV_result = set()

    if start_seg not in G:
        print("\n{}\nstart seg does not in the graph\nstart seg:{}\nend seg:{}\nsegments:{}".format(bubble_path, start_seg, end_seg, segments))
    elif end_seg not in G:
        print("\n{}\nend seg does not in the graph\nstart seg:{}\nend seg:{}\nsegments:{}".format(bubble_path, start_seg, end_seg, segments))
    else:
        print("start determination of path")
        paths = nx.shortest_simple_paths(G, start_seg, end_seg, weight = "weight")
        first_path, first_path_coverage = next(paths, None), -1
        second_path, second_path_coverage = next(paths, None), -1
        print("path get")
        genetic_constitution = -1      #0 = Homo, 1 = Hetero, -1 = reads coverage is not enough
        if second_path == None:
            genetic_constitution = 0
        else:
            first_path_coverage, first_path_length = calculate_path_coverage(first_path, index_gfa)
            second_path_coverage, second_path_length = calculate_path_coverage(second_path, index_gfa)
            two_paths_intersection = set(first_path)
            two_paths_intersection.intersection(set(second_path))
            two_paths_intersection_length = 0
            for inter_seg in two_paths_intersection:
                two_paths_intersection_length += index_gfa[inter_seg][2]
            t1, t2 = None, None
            if first_path_coverage >= max(0.2*average_cover, 3):
                t1 = True
            else:
                t1 = False
            if second_path_coverage >= max(0.2*average_cover, 3):
                t2 = True
            else:
                t2 = False
            if t1 == False:
                pass
            elif t1 == True and t2 == False:
                genetic_constitution = 0
            elif t1 == True and t2 == True:
                if (two_paths_intersection_length/first_path_length >= 0.95) or (two_paths_intersection_length/second_path_length >= 0.95):
                    genetic_constitution = 0
                else:
                    genetic_constitution = 1
        print("path determination finished. {}s".format(time.time() - timep))
        if genetic_constitution == -1:
            pass
        elif genetic_constitution == 0:
            sv_list = mapping_breakpoints_to_path(G, first_path, breakpoint_list, reads_repeat_dict, ref_seg2chr,
                                                  chr2ref_seg, first_path_coverage, percent095_reads_percent, fasta,
                                                  temporary_directory, sequencing_tech)
            for sv in sv_list:
                if sv not in SV_result:
                    SV_result.add(sv)
        elif genetic_constitution == 1:
            sv_list1 = mapping_breakpoints_to_path(G, first_path, breakpoint_list, reads_repeat_dict, ref_seg2chr,
                                                  chr2ref_seg, first_path_coverage, percent095_reads_percent, fasta,
                                                  temporary_directory, sequencing_tech)
            for sv in sv_list1:
                if sv not in SV_result:
                    SV_result.add(sv)
            sv_list2 = mapping_breakpoints_to_path(G, second_path, breakpoint_list, reads_repeat_dict, ref_seg2chr,
                                                  chr2ref_seg, second_path_coverage, percent095_reads_percent, fasta,
                                                  temporary_directory, sequencing_tech)
            for sv in sv_list2:
                temp_list = copy.deepcopy(list(sv))
                if sv not in SV_result:
                    temp_list[9] = "./1"
                    SV_result.add(tuple(temp_list))
                elif sv in SV_result:
                    temp_list[9] = "1/1"
                    SV_result.remove(sv)
                    SV_result.add(tuple(temp_list))
    print("SV finding finished. {}s".format(time.time() - time0))
    return SV_result

#input bubble_path_del_clustered/bubble_path_ins_clustered
def detect_SV(bubble_path_clustered, cluster_result, thread,
              shared_bubble2gfaID, gfa_file, shared_index_gfa, shared_offset_list, shared_ref_seg2chr, shared_chr2ref_seg, average_cover,
              percent095_reads_percent, fasta, temporary_directory, sequencing_tech, reads_repeat_dict = None
              ):
    '''

    :param bubble_path_clustered: dict with bubble_path_clustered as key, bubble_path belonging to it as value
    :param cluster_result: dict with bubble_path as key, breakpoints as value
    :return:
    '''
    vcf_info = set()

    batch = []
    with ProcessPoolExecutor(max_workers=1) as executor:
        for key in bubble_path_clustered:
            tem_list = list()
            tem_list.append(key)
            breakpoint_l = list()
            for bubblep in bubble_path_clustered[key]:
                for breakp in cluster_result[bubblep]:
                    breakpoint_l.append(breakp)
            """
            count = 0
            i = 0
            for i in range(len(key)):
                element = key[i]
                if i == 0:
                    if isinstance(element, tuple) and all(isinstance(j, tuple) for j in element):
                        list1, list2 = shared_bubble2gfaID[element[0][0]][element[0][1]][1], \
                        shared_bubble2gfaID[element[1][0]][element[1][1]][1]
                        if list1[-1] == list2[0]:
                            count += len(list2)
                        elif list1[0] == list2[-1]:
                            count += len(list1)
                    else:
                        list0 = shared_bubble2gfaID[element[0]][element[1]][1]
                        count += len(list0)
                if i == len(key) - 1:
                    if isinstance(element, tuple) and all(isinstance(j, tuple) for j in element):
                        list1, list2 = shared_bubble2gfaID[element[0][0]][element[0][1]][1], \
                        shared_bubble2gfaID[element[1][0]][element[1][1]][1]
                        if list1[-1] == list2[0]:
                            count += len(list1)
                        elif list1[0] == list2[-1]:
                            count += len(list2)
                    else:
                        list0 = shared_bubble2gfaID[element[0]][element[1]][1]
                        count += len(list0)
                if i > 0 and i < len(key) - 1:
                    if isinstance(element, tuple) and all(isinstance(j, tuple) for j in element):
                        list1, list2 = shared_bubble2gfaID[element[0][0]][element[0][1]][1], \
                        shared_bubble2gfaID[element[1][0]][element[1][1]][1]
                        if list1[-1] == list2[0]:
                            count += len(list1)
                        elif list1[0] == list2[-1]:
                            count += len(list2)
            """
            tem_list.append(breakpoint_l)
            batch.append(tem_list)

        length_batch = len(batch)

        print("Total bubble path clustered: {}".format(length_batch))
        futures = [executor.submit(process_cluster, *args, shared_bubble2gfaID, gfa_file, shared_index_gfa,
                                   shared_offset_list, shared_ref_seg2chr, shared_chr2ref_seg, average_cover,
                                   percent095_reads_percent, fasta, temporary_directory, sequencing_tech,
                                   reads_repeat_dict) for args in batch]
        for future in as_completed(futures):
            try:
                result = future.result()      #process_cluster return a set
                vcf_info.update(result)
                length_batch -= 1
                print("remained: {}".format(length_batch))
            except Exception as e:
                print("task failed：{}".format(e))
    return vcf_info

######main
def main():
    parser = argparse.ArgumentParser(description="This program is writen for detect structural variants(SV)" \
                                                 "based on long reads alignment result and pan-genome graph." \
                                                 "")
    parser.add_argument("-i", "--input", type=str, help="alignment to pan-genome in gaf format")
    parser.add_argument("-fa", "--fasta_in", type=str, help="input reads in fasta format")
    parser.add_argument("-g", "--gfa_in", type=str, help="input pan-genome graph in gfa format")
    parser.add_argument("-o", "--output", type=str, help="output_SV")
    parser.add_argument("-t", "--thread", type=int, default=1, help="thread(default=1)")
    parser.add_argument("-tech", "--sequencing_technology", type=str,
                        choices=["rs", "sq", "ccs", "ont"], default="sq",
                        help="Sequencing technology(default=sq, choices=rs/sq/ccs/ont)\tPacbio RSII: rs\tPacbio Sequel: sq\tPacbio CCS: ccs\tOxford Nanopore: ont")
    parser.add_argument("-r", "--repeat", type=int, default=None,
                        help="Input repeatmasker result of reads to decrease the effect of repeats(optional, fa.out)")
    parser.add_argument("-temp", "--temporary_directory", type=str, default="./tmp", help="Temporary directory")

    args = parser.parse_args()

    # Check environment
    current_dir = os.getcwd()
    path_to_gfatools, path_to_wtdbg2, path_to_minimap2 = shutil.which("gfatools"), shutil.which("wtdbg2"), shutil.which("minimap2")
    if path_to_gfatools == None:
        print("Can not find gfatools.")
        sys.exit(1)
    if path_to_wtdbg2 == None:
        print("Can not find wtdbg2.")
        sys.exit(1)
    if path_to_minimap2 == None:
        print("Can not find minimap2.")
    if sys.version_info < (3, 7):
        print("Needs python >= 3.7.")
        sys.exit(1)
    os.environ["PATH"] += (
            os.pathsep + path_to_gfatools +
            os.pathsep + path_to_wtdbg2 +
            os.pathsep + path_to_minimap2
    )

    # Check path/directory
    def check_path(path):
        if os.path.exists(path):
            return True
        else:
            print(f"Path {path} does not exist.")
            return False

    if args.input:
        pass
    else:
        print("Need gaf input, see -i/--input.")
    check_path(args.input)

    if args.fasta_in:
        pass
    else:
        print("Need fasta input, see -fa/--fasta_in.")
    check_path(args.fasta_in)

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

    if args.repeat != None:
        check_path(args.repeat)

    if not os.path.exists(args.temporary_directory):
        os.makedirs(args.temporary_directory, exist_ok=True)
        print("Created temporary directory: {}".format(os.path.abspath(args.temporary_directory)))
    else:
        print("Set temporary directory: {}".format(os.path.abspath(args.temporary_directory)))

    gaf = args.input  # input gaf
    fasta = args.fasta_in  # input fa
    gfa = args.gfa_in  # input gfa
    out_dir = args.output  # output directory
    thread = args.thread  # thread
    sequencing_tech = args.sequencing_technology  # sequencing technology
    repeatmasker_out = args.repeat  # input repeatmasker result
    temporary_directory = os.path.abspath(args.temporary_directory)  # os.path.abspath(args.temporary_directory)

    ######read gfa file
    time1 = time.time()
    time2 = time.time()
    print(time1 - time2)

    print("Indexing gfa.")
    offset_list_gfa, index_dict_gfa, reference_seg2chr, chr2reference_seg = list(), dict(), dict(), dict()  #

    format1 = ""
    if gfa[len(gfa) - 4:] == ".gfa":
        format1 = "gfa"
    elif gfa[len(gfa) - 7:] == ".gfa.gz":
        format1 = "gfa.gz"
    else:
        print("Failed to read graph, graph in unknown format, please input graph in gfa format or gfa.gz format.")
    if format1 == "gfa":
        with open(gfa, "rb") as f_gfa:
            offset_list_gfa, index_dict_gfa, reference_seg2chr, chr2reference_seg = index_gfa(f_gfa)  #

    elif format1 == "gfa.gz":
        with gzip.open(gfa, "rb") as f_gfa:
            offset_list_gfa, index_dict_gfa, reference_seg2chr, chr2reference_seg = index_gfa(f_gfa)  #

    ######call bubbles
    time2 = time1
    time1 = time.time()
    print(time1 - time2)

    print("Calling bubbles.")
    bubble, gfaID2bubble = dict(), dict()      #

    with subprocess.Popen(["gfatools", "bubble", gfa], stderr=subprocess.PIPE, stdout=subprocess.PIPE,
                          universal_newlines=True) as proc:
        for line in proc.stdout:
            info = line.strip().split("\t")
            chr, start, end, segments = info[0], int(info[1]), int(info[2]), info[11].split(",")
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
        del line, info, chr, start, end, segments, segment

    for chrom in bubble:
        bubble[chrom] = dict(sorted(bubble[chrom].items()))
        tuple_tmp = tuple(bubble[chrom].keys())

        # fake head bubble
        index_h = tuple_tmp[0]
        seg_h = bubble[chrom][index_h][1][0]
        bubble[chrom][0] = [0, []]
        gfaID2bubble[seg_h].insert(0, (chrom, 0))

        # fake tail bubble
        index_t = tuple_tmp[-1]
        seg_t = bubble[chrom][index_t][1][-1]
        tail = bubble[chrom][index_t][0] + index_dict_gfa[seg_t][2]
        bubble[chrom][tail] = [tail, []]
        gfaID2bubble[seg_t].append((chrom, tail))

        bubble[chrom] = dict(sorted(bubble[chrom].items()))
    del chrom, tuple_tmp, index_h, index_t, seg_h, seg_t, tail

    with Manager() as manager:
        ######read gaf
        time2 = time1
        time1 = time.time()
        print(time1 - time2)

        print("Reading gaf.")

        breakpoint_list_del_gaf, breakpoint_list_ins_gaf, average_coverage = list(), list(), -1      #

        format2 = ""
        if gaf[len(gaf) - 4:] == ".gaf":
            format2 = "gaf"
        elif gaf[len(gaf) - 7:] == ".gaf.gz":
            format2 = "gaf.gz"
        else:
            print("Failed to read alignment result, alignment result in unknown format, please input it in gaf format or gaf.gz format.")
        if format2 == "gaf":
            with open(gaf, "rt") as f_gaf:
                breakpoint_list_del_gaf, breakpoint_list_ins_gaf, average_coverage = read_gaf(f_gaf, thread, index_dict_gfa, manager)

        elif format2 == "gaf.gz":
            with gzip.open(gaf, "rt", encoding="utf-8") as f_gaf:
                breakpoint_list_del_gaf, breakpoint_list_ins_gaf, average_coverage = read_gaf(f_gaf, thread, index_dict_gfa, manager)

        reads_with_ins_breakpoint = set()  # will be used later to pick a subset of fasta.
        for breakpoint in breakpoint_list_ins_gaf:
            reads_with_ins_breakpoint.add(breakpoint[0])

        ######Cluster breakpoints into bubbles
        ###Determine bubble path of each breakpoint.
        time2 = time1
        time1 = time.time()
        print(time1 - time2)

        print("Clustering breakpoints")

        shared_gfaID2bubble_dict = manager.dict(gfaID2bubble)      #需要修改函数breakpoint_cluster()
        del gfaID2bubble

        cluster_result_del = breakpoint_cluster(breakpoint_list_del_gaf, shared_gfaID2bubble_dict, thread)      #
        cluster_result_ins = breakpoint_cluster(breakpoint_list_ins_gaf, shared_gfaID2bubble_dict, thread)      #
        del breakpoint_list_del_gaf, breakpoint_list_ins_gaf
        del shared_gfaID2bubble_dict

        ######Decrease the effect of repeat(optional)
        time2 = time1
        time1 = time.time()
        print(time1 - time2)

        print("Reading repeatmasker result")

        reads_repeat_dict, total_reads_num, percent095_reads_num = None, -1, -1
        reads_repeat_dict_del, reads_repeat_dict_ins = None, None      #
        percent095_reads_percent = 0      #
        if repeatmasker_out == None:
            print("No repeatmasker result input, skip.")
        elif repeatmasker_out != None:
            reads_repeat_dict, total_reads_num, percent095_reads_num = \
                read_repeatmasker(repeatmasker_out, fasta)      #{readid:percent_repeat, ...}, {readid:bases, }

            percent095_reads_percent = percent095_reads_num / total_reads_num
            reads_repeat_dict_del, reads_repeat_dict_ins = dict(), dict()
            for key in cluster_result_del:
                for breakpoint in cluster_result_del[key]:
                    readID = breakpoint[0]
                    readID_re = readID.split()[0]
                    if readID_re in reads_repeat_dict:
                        reads_repeat_dict_del[readID] = reads_repeat_dict[readID_re]
            for key in cluster_result_ins:
                for breakpoint in cluster_result_ins[key]:
                    readID = breakpoint[0]
                    readID_re = readID.split()[0]
                    if readID_re in reads_repeat_dict:
                        reads_repeat_dict_ins[readID] = reads_repeat_dict[readID_re]
        del reads_repeat_dict, total_reads_num, percent095_reads_num

        ######SV detection
        time2 = time1
        time1 = time.time()
        print(time1 - time2)

        print("Detecting SV-clustering bubble path.")
        print("d")
        # bubble_path_clustered_deletion
        # For deletion, if several bubble paths have overlap, the breakpoints belonging to them may indicate the same
        # SV(deletion), so I combine these bubble path into a "bubble_path_del_clustered".
        # However, there still a possibility that the breakpoints belonging to them indicate several SV(deletion),
        # so this kind of situation should be taken into consideration.
        bubble_path_del_dict = {}      #{bubble_path: [ref_tag, start, end], ...}, ref_tag mains chromosome in most condition.
        for bubble_path in list(cluster_result_del.keys()):
            head_bubble, tail_bubble = bubble_path[0], bubble_path[-1]
            ref_tag, start, end = "", -1, -1
            if is_bubble_or_betweenbubble(head_bubble) == "betweenbubble":
                start = bubble[head_bubble[0][0]][head_bubble[0][1]][0]
                ref_tag = head_bubble[0][0]
            else:
                start = head_bubble[1]
                ref_tag = head_bubble[0]
            if is_bubble_or_betweenbubble(tail_bubble) == "betweenbubble":
                end = tail_bubble[1][1]
            else:
                end = bubble[tail_bubble[0]][tail_bubble[1]][0]
            bubble_path_del_dict[bubble_path] = (ref_tag, start, end)

        bubble_path_del_dict = dict(sorted(bubble_path_del_dict.items(), key=lambda item: (item[1][0], item[1][1], -item[1][2])))
        bubble_path_del_list = list(bubble_path_del_dict.keys())

        bubble_path_del_clustered = {}      #{bubble_path_clustered: [bubble_path1, ...], ...}      #
        i = 0
        now_bubble_path_clustered = []
        now_bubble_path = []
        for i in range(len(bubble_path_del_list)):
            key = bubble_path_del_list[i]
            if now_bubble_path_clustered == []:
                now_bubble_path_clustered = list(key)
                now_bubble_path = [key]
            elif now_bubble_path_clustered != []:
                set1 = set(now_bubble_path_clustered)
                set1.update(set(key))
                if len(set1) < len(now_bubble_path_clustered) + len(key):
                    for element in key:
                        if element not in now_bubble_path_clustered:
                            now_bubble_path_clustered.append(element)
                    now_bubble_path.append(key)
                elif len(set1) == len(now_bubble_path_clustered) + len(key):
                    bubble_path_del_clustered[tuple(now_bubble_path_clustered)] = now_bubble_path
                    now_bubble_path_clustered = list(key)
                    now_bubble_path = [key]
            if i == len(bubble_path_del_list) - 1:
                bubble_path_del_clustered[tuple(now_bubble_path_clustered)] = now_bubble_path
        del bubble_path_del_dict, bubble_path_del_list

        print("i")
        # bubble_path_clustered_insertion
        # For insertion, since the span on pan-genome of a breakpoint of a SV(insertion) is usually very small,
        # most of the time less than several bps(even 0 bp), the cluster strategy of deletion is not suitable for
        # insertion, however.
        # So I need a new cluster strategy.
        # First, I get the start and end of each bubble path, then cluster all the bubble path whose distance between
        # each other is less than 50bp.
        bubble_path_ins_dict = {}      #{bubble_path: [ref_tag, start, end], ...}
        for bubble_path in list(cluster_result_ins.keys()):
            head_bubble, tail_bubble = bubble_path[0], bubble_path[-1]
            ref_tag, start, end = "", -1, -1  # ref_tag means chromosome.
            if isinstance(head_bubble, tuple) and all(isinstance(i, tuple) for i in head_bubble):
                start = bubble[head_bubble[0][0]][head_bubble[0][1]][0]
                ref_tag = head_bubble[0][0]
            else:
                start = head_bubble[1]
                ref_tag = head_bubble[0]
            if isinstance(tail_bubble, tuple) and all(isinstance(i, tuple) for i in tail_bubble):
                end = tail_bubble[1][1]
            else:
                end = bubble[tail_bubble[0]][tail_bubble[1]][0]
            bubble_path_ins_dict[bubble_path] = (ref_tag, start, end)

        bubble_path_ins_dict = dict(sorted(bubble_path_ins_dict.items(), key=lambda item: (item[1][0], item[1][1], -item[1][2])))
        bubble_path_ins_list = list(bubble_path_ins_dict.keys())

        bubble_path_ins_clustered_ = {}      #{bubble_path_clustered: [bubble_path1, ...], ...}      #
        i = 0
        now_bubble_path_clustered = []
        now_ref_tag, now_start, now_end = "", -1, -1
        now_bubble_path = []
        for i in range(len(bubble_path_ins_list)):
            key = bubble_path_ins_list[i]
            if now_bubble_path_clustered == []:
                now_bubble_path_clustered = list(key)
                now_bubble_path = [key]
                now_ref_tag, now_start, now_end = bubble_path_ins_dict[key][0], bubble_path_ins_dict[key][1], \
                bubble_path_ins_dict[key][2]
            elif now_bubble_path_clustered != []:
                if bubble_path_ins_dict[key][0] == now_ref_tag:
                    if bubble_path_ins_dict[key][1] >= now_start and bubble_path_ins_dict[key][2] <= now_end:
                        now_bubble_path.append(key)
                    elif bubble_path_ins_dict[key][1] >= now_start and bubble_path_ins_dict[key][1] <= now_end and bubble_path_ins_dict[key][2] >= now_end:
                        for element in key:
                            if element not in now_bubble_path_clustered:
                                now_bubble_path_clustered.append(element)
                        now_end = bubble_path_ins_dict[key][2]
                        now_bubble_path.append(key)
                    elif bubble_path_ins_dict[key][1] >= now_end:
                        if bubble_path_ins_dict[key][1] - now_end <= 50:
                            now_end = bubble_path_ins_dict[key][2]
                            now_bubble_path.append(key)
                            list_1 = list(bubble[now_ref_tag].keys())
                            index1, index2 = -1, -1
                            end_now_bubble_path_clustered, start_key = now_bubble_path_clustered[-1], key[0]
                            if is_bubble_or_betweenbubble(end_now_bubble_path_clustered) == "bubble":
                                if is_bubble_or_betweenbubble(start_key) == "bubble":
                                    index1, index2 = list_1.index(end_now_bubble_path_clustered[1]), list_1.index(
                                        start_key[1])
                                    now_bubble_path_clustered.append(
                                        ((now_ref_tag, list_1[index1]), (now_ref_tag, list_1[index1 + 1])))
                                    j = 0
                                    for j in range(index1 + 1, index2):
                                        now_bubble_path_clustered.append((now_ref_tag, list_1[j]))
                                        now_bubble_path_clustered.append(
                                            ((now_ref_tag, list_1[j]), (now_ref_tag, list_1[j + 1])))
                                    for element in key:
                                        now_bubble_path_clustered.append(element)
                                elif is_bubble_or_betweenbubble(start_key) == "betweenbubble":
                                    index1, index2 = list_1.index(end_now_bubble_path_clustered[1]), list_1.index(
                                        start_key[0][1])
                                    j = 0
                                    for j in range(index1, index2):
                                        now_bubble_path_clustered.append(
                                            ((now_ref_tag, list_1[j]), (now_ref_tag, list_1[j + 1])))
                                        now_bubble_path_clustered.append((now_ref_tag, list_1[j + 1]))
                                    for element in key:
                                        now_bubble_path_clustered.append(element)
                            elif is_bubble_or_betweenbubble(end_now_bubble_path_clustered) == "betweenbubble":
                                if is_bubble_or_betweenbubble(start_key) == "bubble":
                                    index1, index2 = list_1.index(end_now_bubble_path_clustered[1][1]), list_1.index(
                                        start_key[1])
                                    j = 0
                                    for j in range(index1, index2):
                                        now_bubble_path_clustered.append((now_ref_tag, list_1[j]))
                                        now_bubble_path_clustered.append(
                                            ((now_ref_tag, list_1[j]), (now_ref_tag, list_1[j + 1])))
                                elif is_bubble_or_betweenbubble(start_key) == "betweenbubble":
                                    index1, index2 = list_1.index(end_now_bubble_path_clustered[1][1]), list_1.index(
                                        start_key[0][1])
                                    now_bubble_path_clustered.append((now_ref_tag, list_1[index1]))
                                    if index1 == index2:
                                        pass
                                    elif index1 < index2:
                                        j = 0
                                        for j in range(index1 + 1, index2 + 1):
                                            now_bubble_path_clustered.append(
                                                ((now_ref_tag, list_1[j - 1]), (now_ref_tag, list_1[j])))
                                            now_bubble_path_clustered.append((now_ref_tag, list_1[j]))
                                    for element in key:
                                        now_bubble_path_clustered.append(element)
                        else:
                            bubble_path_ins_clustered_[tuple(now_bubble_path_clustered)] = now_bubble_path
                            now_bubble_path_clustered = list(key)
                            now_bubble_path = [key]
                            now_ref_tag, now_start, now_end = bubble_path_ins_dict[key][0], bubble_path_ins_dict[key][
                                1], bubble_path_ins_dict[key][2]

                elif bubble_path_ins_dict[key][0] != now_ref_tag:
                    bubble_path_ins_clustered_[tuple(now_bubble_path_clustered)] = now_bubble_path
                    now_bubble_path_clustered = list(key)
                    now_bubble_path = [key]
                    now_ref_tag, now_start, now_end = bubble_path_ins_dict[key][0], bubble_path_ins_dict[key][1], \
                    bubble_path_ins_dict[key][2]
            if i == len(bubble_path_ins_list) - 1:
                bubble_path_ins_clustered_[tuple(now_bubble_path_clustered)] = now_bubble_path
        del bubble_path_ins_dict, bubble_path_ins_list

        print("e")
        # expand bubble_path for the subsequent local assembly and contig to part path alignment
        bubble_path_ins_clustered = {}
        bubble_path_clustered_list = list(bubble_path_ins_clustered_.keys())
        for bubble_path_clustered in bubble_path_clustered_list:
            '''bubble_path_clustered a tuple contain two kind of elements --> 
               (chr1, start1) represents a bubble;
               ((chr1, start1), (chr2, start2)) represents a ref segment between two bubbles.
            '''
            # first expand bubble_path_clustered to two side ending with segments between two bubbles((, ), (, ))
            expanded_bubble_path_clustered = list(bubble_path_clustered)
            head_bubble, tail_bubble = bubble_path_clustered[0], bubble_path_clustered[-1]
            ref_tag, start, end = "", -1, -1  # ref_tag means chromosome.
            if isinstance(head_bubble, tuple) and all(isinstance(i, tuple) for i in head_bubble):
                ref_tag = head_bubble[0][0]
            else:
                ref_tag = head_bubble[0]
            tuple_bubble_on_ref_tag = tuple(bubble[ref_tag].keys())  ###
            length_tuple_bubble_on_ref_tag = len(tuple_bubble_on_ref_tag)  ###
            if isinstance(head_bubble, tuple) and all(isinstance(i, tuple) for i in head_bubble):
                start = bubble[head_bubble[0][0]][head_bubble[0][1]][0]
            else:
                index0 = tuple_bubble_on_ref_tag.index(head_bubble[1])
                expanded_bubble_path_clustered.insert(0, (
                (ref_tag, tuple_bubble_on_ref_tag[index0 - 1]), (ref_tag, tuple_bubble_on_ref_tag[index0])))
                start = bubble[ref_tag][tuple_bubble_on_ref_tag[index0 - 1]][0]
            if isinstance(tail_bubble, tuple) and all(isinstance(i, tuple) for i in tail_bubble):
                end = tail_bubble[1][1]
            else:
                index8 = tuple_bubble_on_ref_tag.index(tail_bubble[1])
                expanded_bubble_path_clustered.append(
                    ((ref_tag, tuple_bubble_on_ref_tag[index8]), (ref_tag, tuple_bubble_on_ref_tag[index8 + 1])))
                end = tuple_bubble_on_ref_tag[index8 + 1]

            index1 = tuple_bubble_on_ref_tag.index(expanded_bubble_path_clustered[0][0][1])
            index2 = tuple_bubble_on_ref_tag.index(expanded_bubble_path_clustered[-1][1][1])
            while end - start < 5000:
                if index1 > 0:
                    expanded_bubble_path_clustered.insert(0, (ref_tag, tuple_bubble_on_ref_tag[index1]))
                    expanded_bubble_path_clustered.insert(0, (
                    (ref_tag, tuple_bubble_on_ref_tag[index1 - 1]), (ref_tag, tuple_bubble_on_ref_tag[index1])))
                    start = bubble[ref_tag][tuple_bubble_on_ref_tag[index1 - 1]][0]
                    index1 = index1 - 1
                if index2 < length_tuple_bubble_on_ref_tag - 1:
                    expanded_bubble_path_clustered.append((ref_tag, tuple_bubble_on_ref_tag[index2]))
                    expanded_bubble_path_clustered.append(
                        ((ref_tag, tuple_bubble_on_ref_tag[index2]), (ref_tag, tuple_bubble_on_ref_tag[index2 + 1])))
                    end = tuple_bubble_on_ref_tag[index2 + 1]
                    index2 = index2 + 1
                if index1 <= 0 and index2 >= length_tuple_bubble_on_ref_tag - 1:
                    break
            bubble_path_ins_clustered[tuple(expanded_bubble_path_clustered)] = copy.deepcopy(
                bubble_path_ins_clustered_[bubble_path_clustered])

        del bubble_path_ins_clustered_, bubble_path_clustered_list

        time2 = time1
        time1 = time.time()
        print(time1 - time2)

        print("Detecting SV-subfasta.")
        #write a dict containing readids and bases, only contain reads containing insertion breakpoint(s).
        dict_reads_with_ins_breakpoint_fasta = dict()

        total_count = len(reads_with_ins_breakpoint)

        with open("/lustre/home/acct-clswcc/clswcc-fsy/Shuyang/pangenomes/HPRC_MC_CHM13_pangraph/script/simulation_real_and_simu_SV/long_read_generation/reads/longreads.1.out/sim.srt.ont.test0.fa", "w") as fo:
            with gzip.open(fasta, "rt") as fi:
                t, reads, bases = False, "", ""
                c = 0
                for eachline in fi:
                    text = eachline.strip()
                    if text[0] == ">":
                        name = text[1:]
                        if t == True:
                            dict_reads_with_ins_breakpoint_fasta[reads] = bases
                            fo.write(">{}\n{}\n".format(reads, bases))
                            c += 1
                            t, reads, bases = False, "", ""
                        if c == total_count:
                            break
                        if name in reads_with_ins_breakpoint:
                            t = True
                            reads = name
                    else:
                        if t == True:
                            bases += text
                if t == True:
                    dict_reads_with_ins_breakpoint_fasta[reads] = bases
                    fo.write(">{}\n{}\n".format(reads, bases))

        time2 = time1
        time1 = time.time()
        print(time1 - time2)

        print("Detecting SV-detect sv.")
        # input bubble_path_del_clustered/bubble_path_ins_clustered to function detect_SV()
        shared_index_dict_gfa = manager.dict(index_dict_gfa)      #需要修改detect_SV()函数
        del index_dict_gfa
        shared_bubble2gfaID = manager.dict()
        for k, v in bubble.items():
            if isinstance(v, dict):
                shared_bubble2gfaID[k] = manager.dict(v)
            else:
                shared_bubble2gfaID[k] = v
        del bubble
        shared_offset_list = manager.list(offset_list_gfa)
        del offset_list_gfa
        shared_ref_seg2chr = manager.dict(reference_seg2chr)
        del reference_seg2chr
        shared_chr2ref_seg = manager.dict(chr2reference_seg)
        del chr2reference_seg
        shared_dict_reads_with_ins_breakpoint_fasta = manager.dict(dict_reads_with_ins_breakpoint_fasta)
        del dict_reads_with_ins_breakpoint_fasta

        ins_vcf_info = detect_SV(bubble_path_ins_clustered, cluster_result_ins, thread,
                                 shared_bubble2gfaID, gfa, shared_index_dict_gfa, shared_offset_list,
                                 shared_ref_seg2chr, shared_chr2ref_seg,
                                 average_coverage,
                                 percent095_reads_percent, shared_dict_reads_with_ins_breakpoint_fasta,
                                 temporary_directory, sequencing_tech,
                                 reads_repeat_dict_ins
                                 )
        del_vcf_info = detect_SV(bubble_path_del_clustered, cluster_result_del, thread,
                                 shared_bubble2gfaID, gfa, shared_index_dict_gfa, shared_offset_list, shared_ref_seg2chr, shared_chr2ref_seg,
                                 average_coverage,
                                 percent095_reads_percent, shared_dict_reads_with_ins_breakpoint_fasta, temporary_directory, sequencing_tech,
                                 reads_repeat_dict_del
                                 )


        del shared_offset_list, shared_bubble2gfaID

        ######output vcf
        time2 = time1
        time1 = time.time()
        print(time1 - time2)

        print("Outputing vcf/gfa.")
        seg_need_cut = dict()  # in this dict, all seg considered as forward, corresponding to their strand as segment in gfa
        new_link = set()
        new_seg = dict()  # insertions form new segments

        out_vcf_dir = os.path.join(out_dir, "NEW_SV.vcf")
        with open(out_vcf_dir, "w") as fo_vcf:
            fo_vcf.write("##fileformat=VCFv4.2\n")
            fo_vcf.write("##source=your_tool_or_script_name\n")
            fo_vcf.write("reference={}\n".format(gfa))
            fo_vcf.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
            fo_vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            fo_vcf.write(
                '##INFO=<ID=ORNOR,Number=1,Type=String,Description="Structural variant on reference path or not on reference path(in bubbles).">\n')
            fo_vcf.write('##INFO=<ID=PATH,Number=1,Type=String,Description="Path where structural variant occurs">\n')
            fo_vcf.write(
                '##INFO=<ID=START,Number=1,Type=Integer,Description="Start coordinate of structural variant on path">\n')
            fo_vcf.write(
                '##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate of structural variant on path">\n')
            fo_vcf.write(
                '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant, either ins or del">\n')
            for chromsome in shared_chr2ref_seg:
                final_seg = next(reversed(shared_chr2ref_seg[chromsome]))
                so, length_final_seg = shared_chr2ref_seg[chromsome][final_seg], shared_index_dict_gfa[final_seg][2]
                fo_vcf.write('##contig=<ID={},length={}>\n'.format(chromsome, so + length_final_seg))
            del shared_index_dict_gfa
            fo_vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
            # deletion
            for element in del_vcf_info:
                CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO_TUPLE, FORMAT, SAMPLE, NEW_SEG_LINK = element
                path_seg, path_strand = re.split(r"[<>]+", INFO_TUPLE[1]), re.findall(r"[<>]+", INFO_TUPLE[1])
                ID = "{}_{}_{}{}{}{}_{}_{}".format(CHROM, POS, path_strand[0], path_seg[1], path_strand[-1],
                                                   path_seg[-1], INFO_TUPLE[2], INFO_TUPLE[3])
                fo_vcf.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\tORNOR={};PATH={};START={};END={};SVTYPE={}\t{}\t{}\n".format(
                    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO_TUPLE[0], INFO_TUPLE[1], INFO_TUPLE[2], INFO_TUPLE[3],
                    INFO_TUPLE[4], FORMAT, SAMPLE
                ))

                seg1, strand1, pos1, length1, seg2, strand2, pos2, length2 = NEW_SEG_LINK
                if seg1 not in seg_need_cut:
                    seg_need_cut[seg1] = set()
                    seg_need_cut[seg1].add(0)
                    seg_need_cut[seg1].add(length1)
                if seg2 not in seg_need_cut:
                    seg_need_cut[seg2] = set()
                    seg_need_cut[seg2].add(0)
                    seg_need_cut[seg2].add(length2)
                if strand1 == "+":
                    seg_need_cut[seg1].add(pos1)
                elif strand1 == "-":
                    seg_need_cut[seg1].add(length1 - pos1)
                if strand2 == "+":
                    seg_need_cut[seg2].add(pos2)
                elif strand2 == "-":
                    seg_need_cut[seg2].add(length2 - pos2)
                new_link.add(NEW_SEG_LINK)
            # insertion
            for element in ins_vcf_info:
                CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO_TUPLE, FORMAT, SAMPLE, NEW_SEG_LINK = element
                path_seg, path_strand = re.split(r"[<>]+", INFO_TUPLE[1]), re.findall(r"[<>]+", INFO_TUPLE[1])
                ID = "{}_{}_{}{}{}{}_{}_{}".format(CHROM, POS, path_strand[0], path_seg[1], path_strand[-1],
                                                   path_seg[-1], INFO_TUPLE[2], INFO_TUPLE[3])
                fo_vcf.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\tORNOR={};PATH={};START={};END={};SVTYPE={}\t{}\t{}\n".format(
                    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO_TUPLE[0], INFO_TUPLE[1], INFO_TUPLE[2], INFO_TUPLE[3],
                    INFO_TUPLE[4], FORMAT, SAMPLE
                ))

                seg1, strand1, pos1, length1, seg2, strand2, pos2, length2 = NEW_SEG_LINK
                if seg1 not in seg_need_cut:
                    seg_need_cut[seg1] = set()
                    seg_need_cut[seg1].add(0)
                    seg_need_cut[seg1].add(length1)
                if seg2 not in seg_need_cut:
                    seg_need_cut[seg2] = set()
                    seg_need_cut[seg2].add(0)
                    seg_need_cut[seg2].add(length2)
                if strand1 == "+":
                    seg_need_cut[seg1].add(pos1)
                elif strand1 == "-":
                    seg_need_cut[seg1].add(length1 - pos1)
                if strand2 == "+":
                    seg_need_cut[seg2].add(pos2)
                elif strand2 == "-":
                    seg_need_cut[seg2].add(length2 - pos2)
                new_link.add(NEW_SEG_LINK)
                new_seg[NEW_SEG_LINK] = ALT

        for seg in seg_need_cut:
            seg_need_cut[seg] = tuple(sorted(seg_need_cut[seg]))


        out_gfa_expanded_dir = os.path.join(out_dir, "EXPANDED_GRAPH.gfa")
        with open(out_gfa_expanded_dir, "w") as fo_gfa:
            # new segments result from insertion
            for ele in new_seg:
                fo_gfa.write("S\tInsertionPart_{}_{}_{}_{}\t{}\n".format(
                    ele[0], ele[2], ele[4], ele[6], new_seg[ele]
                ))
            # new segments and links result from SV
            with open(gfa, "rt") as fi_gfa:
                for eachline in fi_gfa:
                    text = eachline.strip()
                    if text[0] == "S":
                        text_split = text.split("\t")
                        seg_id, seg_base = text_split[1], text_split[2]
                        if seg_id in seg_need_cut:
                            cut_pos = seg_need_cut[seg_id]
                            if seg_id in shared_ref_seg2chr:
                                SN = ""
                                SO = shared_chr2ref_seg[shared_ref_seg2chr[seg_id]][seg_id]
                                for part in text_split:
                                    if part.startswith("SN"):
                                        SN = part
                                        break
                                i = 0
                                for i in range(len(cut_pos) - 1):
                                    fo_gfa.write("S\t{}_{}_{}\t{}\t{}\tSO:i:{}\t{}\n".format(
                                        seg_id, cut_pos[i], cut_pos[i + 1], seg_base[cut_pos[i]: cut_pos[i + 1]],
                                        SN, SO + cut_pos[i], "SR:i:0"
                                    ))

                            else:
                                i = 0
                                for i in range(len(cut_pos) - 1):
                                    fo_gfa.write("S\t{}_{}_{}\t{}\n".format(
                                        seg_id, cut_pos[i], cut_pos[i + 1], seg_base[cut_pos[i]: cut_pos[i + 1]]
                                    ))

                        elif text_split[1] not in seg_need_cut:
                            fo_gfa.write(eachline)
                    elif text[0] == "L":
                        text_split = text.split("\t")
                        seg_s, strand_s, seg_e, strand_e = text_split[1], text_split[2], text_split[3], text_split[4]
                        seg_s_, seg_e_ = "", ""
                        if seg_s not in seg_need_cut and seg_e not in seg_need_cut:
                            fo_gfa.write(eachline)
                        elif seg_s in seg_need_cut and seg_e not in seg_need_cut:
                            if strand_s == "+":
                                seg_s_ = "{}_{}_{}".format(seg_s, seg_need_cut[seg_s][len(seg_need_cut[seg_s]) - 2],
                                                           seg_need_cut[seg_s][-1])
                            elif strand_s == "-":
                                seg_s_ = "{}_{}_{}".format(seg_s, seg_need_cut[seg_s][0], seg_need_cut[seg_s][1])
                            fo_gfa.write("L\t{}\t{}\t{}\t{}\tOM\n".format(seg_s_, strand_s, seg_e_, strand_e))
                        elif seg_s not in seg_need_cut and seg_e in seg_need_cut:
                            if strand_e == "+":
                                seg_e_ = "{}_{}_{}".format(seg_e, seg_need_cut[seg_e][0], seg_need_cut[seg_e][1])
                            elif strand_e == "-":
                                seg_e_ = "{}_{}_{}".format(seg_e, seg_need_cut[seg_e][len(seg_need_cut[seg_e]) - 2],
                                                           seg_need_cut[seg_e][-1])
                            fo_gfa.write("L\t{}\t{}\t{}\t{}\tOM\n".format(seg_s_, strand_s, seg_e_, strand_e))
                        elif seg_s in seg_need_cut and seg_e in seg_need_cut:
                            if strand_s == "+":
                                seg_s_ = "{}_{}_{}".format(seg_s, seg_need_cut[seg_s][len(seg_need_cut[seg_s]) - 2],
                                                           seg_need_cut[seg_s][-1])
                            elif strand_s == "-":
                                seg_s_ = "{}_{}_{}".format(seg_s, seg_need_cut[seg_s][0], seg_need_cut[seg_s][1])
                            if strand_e == "+":
                                seg_e_ = "{}_{}_{}".format(seg_e, seg_need_cut[seg_e][0], seg_need_cut[seg_e][1])
                            elif strand_e == "-":
                                seg_e_ = "{}_{}_{}".format(seg_e, seg_need_cut[seg_e][len(seg_need_cut[seg_e]) - 2],
                                                           seg_need_cut[seg_e][-1])
                            fo_gfa.write("L\t{}\t{}\t{}\t{}\tOM\n".format(seg_s_, strand_s, seg_e_, strand_e))
                    else:
                        break
            # new links in a segments which is cut.
            for seg_nc in seg_need_cut:
                cut_pos = seg_need_cut[seg_nc]
                if len(cut_pos) >= 3:
                    i = 0
                    for i in range(len(cut_pos) - 2):
                        fo_gfa.write("L\t{}_{}_{}\t+\t{}_{}_{}\t+\t0M\n".format(
                            seg_nc, cut_pos[i], cut_pos[i + 1], seg_nc, cut_pos[i + 1], cut_pos[i + 2]
                        ))
            # new links result from SV
            for elem in new_link:
                seg1, strand1, pos1, length1, seg2, strand2, pos2, length2 = elem
                start_, end_ = "", ""
                if strand1 == "+":
                    index1 = seg_need_cut[seg1].index(pos1)
                    start_ = "{}_{}_{}".format(seg1, index1 - 1, index1)
                elif strand1 == "-":
                    index1 = seg_need_cut[seg1].index(length1 - pos1)
                    start_ = "{}_{}_{}".format(seg1, index1, index1 + 1)
                if strand2 == "+":
                    index2 = seg_need_cut[seg2].index(pos2)
                    end_ = "{}_{}_{}".format(seg2, index2, index2 + 1)
                elif strand2 == "-":
                    index2 = seg_need_cut[seg2].index(length2 - pos2)
                    end_ = "{}_{}_{}".format(seg2, index2 - 1, index2)
                # del
                if elem not in new_seg:
                    fo_gfa.write("L\t{}\t{}\t{}\t{}\t0M\n".format(start_, strand1, end_, strand2))
                # ins
                elif elem in new_seg:
                    name_alt = "InsertionPart_{}_{}_{}_{}".format(elem[0], elem[2], elem[4], elem[6])
                    fo_gfa.write("L\t{}\t{}\t{}\t{}\t0M\n".format(start_, strand1, name_alt, "+"))
                    fo_gfa.write("L\t{}\t{}\t{}\t{}\t0M\n".format(name_alt, "+", end_, strand2))
        del shared_chr2ref_seg, shared_ref_seg2chr

        time2 = time1
        time1 = time.time()
        print(time1 - time2)

        print("Finish.")

if __name__ == "__main__":
    main()



