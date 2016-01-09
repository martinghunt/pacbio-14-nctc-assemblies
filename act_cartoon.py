#!/usr/bin/env python3

import os
import tempfile
import shutil
import argparse
import copy
import pyfastaq
import pymummer


def run_nucmer(ref, qry, outfile, min_id=98, simplify=True, min_length=250):
    n = pymummer.nucmer.Runner(
        ref,
        qry,
        outfile,
        min_id=min_id,
        min_length=250,
        breaklen=500,
        maxmatch=True,
        simplify=simplify,
        verbose=True,
    )
    n.run()


def load_nucmer_hits(infile, keep_all=False):
    '''Gets nucmer hits from file. Returns dict of query name -> list of nucmer hits'''
    hits = {}
    file_reader = pymummer.coords_file.reader(infile)

    for al in file_reader:
        if keep_all or al.hit_length_qry >= 5000 or (al.hit_length_qry >= 0.03 * al.qry_length):
            if al.qry_name not in hits:
                hits[al.qry_name] = []
            hits[al.qry_name].append(al)

    return hits


def assembly_to_list(filename):
    file_reader = pyfastaq.sequences.file_reader(filename)
    seqs = [copy.copy(x) for x in file_reader]
    for seq in seqs:
        seq.id = seq.id.split()[0]
    return seqs



def get_repeat_coords(infile, min_id=98):
    '''Returns dict of repeats sequences found by nucmer. key=sequence name. value=list of coords'''
    tmpdir = tempfile.mkdtemp(prefix='tmp.run_nucmer.', dir=os.getcwd())
    coords_file = os.path.join(tmpdir, 'nucmer.coords')
    run_nucmer(infile, infile, coords_file, min_id=98, simplify=False, min_length=100)
    nucmer_hits = load_nucmer_hits(coords_file, keep_all=True)
    shutil.rmtree(tmpdir)
    repeat_coords = {}
    for refseq in nucmer_hits:
        nucmer_hits[refseq] = [x for x in nucmer_hits[refseq] if x.ref_name == x.qry_name and x.ref_start != x.qry_start]
    return nucmer_hits


def make_ref_repeat_heatmap(ref_contigs, ref_file, min_id=98):
    repeat_hits = get_repeat_coords(ref_file, min_id=min_id)
    heatmap_data = get_ref_heatmap_data(ref_contigs, repeat_hits, binary=True)
    return heatmap_data


def match_qry_and_ref_contigs(hits):
    qry2ref = {}
    ref2qry = {}
    qry_to_reverse = set()

    for qry_name in hits:
        best_hit = None
        for hit in hits[qry_name]:
            if best_hit is None or hit.hit_length_qry > best_hit.hit_length_qry:
                best_hit = hit
        assert best_hit is not None
        qry2ref[qry_name] = best_hit
        if not best_hit.on_same_strand():
            qry_to_reverse.add(qry_name)

    for qry_name, hit in qry2ref.items():
        if hit.ref_name not in ref2qry:
            ref2qry[hit.ref_name] = set()
        ref2qry[hit.ref_name].add(qry_name)

    return qry2ref, ref2qry, qry_to_reverse


def get_qry_contig_order(ref_contigs, qry_contigs, ref2qry, qry2ref):
    ordered_qry_names = []
    unused_qry_names = set(qry_contigs.keys())

    for ref_contig in ref_contigs:
        ref_name = ref_contig.id
        if ref_name in ref2qry:
            qry_names = ref2qry[ref_name]
            qry_hits = [qry2ref[qry_name] for qry_name in qry_names]
            qry_hits = sorted(qry_hits, key=lambda hit: hit.ref_start)
            ordered_qry_names.extend([x.qry_name for x in qry_hits])
            for hit in qry_hits:
                try:
                    unused_qry_names.remove(hit.qry_name)
                except:
                    pass

    unused_qry_names = sorted(list(unused_qry_names), key=lambda x: len(qry_contigs[x]), reverse=True)
    return ordered_qry_names + unused_qry_names


def get_ref_contig_coords(ref_contigs, contig_space=5, x_offset=0):
    coords = {}
    running_total = x_offset
    total_length = sum([len(contig) for contig in ref_contigs]) + contig_space * (len(ref_contigs) - 1)
    max_x_value = 0
    for contig in ref_contigs:
        scaled_length = 100.0 * len(contig) / total_length
        coords[contig.id] = (running_total, running_total + scaled_length)
        max_x_value = max(max_x_value, running_total + scaled_length)
        running_total += scaled_length + contig_space
    return coords, max_x_value


def heatmap_to_rectangles(l):
    rectangles = [[1, 1, l[0]]]
    for i in range(len(l)):
        depth = l[i]
        if depth != rectangles[-1][-1]:
            rectangles.append([i, i+1, l[i]])
        else:
            rectangles[-1][1] = i
    return rectangles


def get_ref_heatmap_data(ref_contigs, nucmer_hits, binary=False):
    d = {}
    max_depth = 0
    for contig in ref_contigs:
        depths = [0] * len(contig)
        for hits in nucmer_hits.values():
            for hit in hits:
                if hit.ref_name == contig.id:
                    coords = hit.ref_coords()
                    for i in range(coords.start, coords.end + 1):
                        if binary:
                            depths[i] = 1
                        else:
                            depths[i] += 1

        max_depth = max(max_depth, max(depths))
        d[contig.id] = heatmap_to_rectangles(depths)

    return d


def heatmap_max_depth(heatmap):
    max_depth = 0
    for contig_name, rectangles in heatmap.items():
        for l in rectangles:
            max_depth = max(max_depth, l[-1])
    return max_depth


def normalise_heatmap(heatmap, max_depth):
    for contig_name, rectangles in heatmap.items():
        for l in rectangles:
            l[-1] /= max_depth


def normalise_heatmaps(heatmap1, heatmap2):
    max_depth = max(heatmap_max_depth(heatmap1), heatmap_max_depth(heatmap2))
    normalise_heatmap(heatmap1, max_depth)
    normalise_heatmap(heatmap2, max_depth)


def get_hit_and_contig_coords(hits, ref_contigs, qry_contigs, contig_space=2, lower_contigs_offset=0):
    qry2ref, ref2qry, qry_to_reverse = match_qry_and_ref_contigs(hits)
    qry_names_in_order = get_qry_contig_order(ref_contigs, qry_contigs, ref2qry, qry2ref)
    contig_coords = []
    hit_coords = []
    remaining_qry_names = set(qry_contigs.keys())
    total_ref_length = sum([len(x) for x in ref_contigs])
    total_qry_length = sum([len(x) for x in qry_contigs.values()])
    ref_contig_coords, x = get_ref_contig_coords(ref_contigs, x_offset=options.lower_contigs_offset)
    max_x_value = 0

    for qry_name in qry_names_in_order:
        remaining_qry_names.remove(qry_name)
        qry_length = len(qry_contigs[qry_name])
        scaled_qry_length = 100.0 * qry_length / total_ref_length
        if len(contig_coords) == 0:
            qry_contig_start = 0
        else:
            qry_contig_start = contig_coords[-1][-1] + contig_space
        qry_contig_start += lower_contigs_offset

        contig_coords.append((qry_name, qry_contig_start, qry_contig_start + scaled_qry_length))

        if qry_name in hits:
            for hit in hits[qry_name]:
                if qry_name in qry_to_reverse:
                    hit.reverse_query()
                qry_coords = hit.qry_coords()
                ref_coords = hit.ref_coords()
                ref_contig_start = ref_contig_coords[hit.ref_name][0]
                qry_start = qry_contig_start + 100.0 * qry_coords.start / total_ref_length
                qry_end = qry_start + 100.0 * len(qry_coords) / total_ref_length
                ref_start = ref_contig_start + 100.0 * ref_coords.start / total_ref_length
                ref_end = ref_start + 100.0 * len(ref_coords) / total_ref_length
                max_x_value = max(max_x_value, qry_end, ref_end)

                if hit.on_same_strand():
                    coords = (qry_start, qry_end, ref_end, ref_start, True)
                else:
                    coords = (qry_end, qry_start, ref_end, ref_start, False)

                hit_coords.append(coords)

                if qry_name in qry_to_reverse:
                    hit.reverse_query()
        else:
            max_x_value = max(max_x_value, qry_contig_start + scaled_qry_length)

    assert len(contig_coords) == len(qry_contigs)
    return  contig_coords, hit_coords, qry_to_reverse, max_x_value


# coords = list of tuples [(x1, y1), (x2, y2) ...]
def svg_polygon(coords, fill_colour, border_colour, border_width = 1, opacity=-1):
    return_string = '<polygon points="' + ' '.join([str(x[0])+','+str(x[1]) for x in coords]) + '" ' \
        + 'fill="' + fill_colour + '" '

    if opacity != -1:
        return_string += 'fill-opacity="' + str(opacity) + '" '

    return_string += 'stroke="' + border_colour + '" ' \
                     + 'stroke-width="' + str(border_width) + '" ' \
                     + '/>'
    return return_string


def write_lineplot(data, y_top, y_bottom, x_min, x_max, ref_contig_length, filehandle, colour):
    total_width = x_max - x_min
    total_height = y_bottom - y_top
    last_point = None
    polygon_points = []
    for start, end, height in data:
        plot_start = x_min + (total_width * start / ref_contig_length)
        plot_end = x_min + (total_width * end / ref_contig_length)
        y = y_bottom - (height * total_height)
        if len(polygon_points) == 0:
            polygon_points.append((plot_start, y_bottom))
        polygon_points.append((plot_start, y))
        polygon_points.append((plot_end, y))
        last_x = plot_end

    polygon_points.append((last_x, y_bottom))
    print(svg_polygon(polygon_points, colour, colour, border_width=0.1), file=filehandle)


def write_heatmap(data, y_bottom, y_top, x_min, x_max, ref_contig_length, filehandle, binary_blue=False):
    total_width = x_max - x_min
    for start, end, depth in data:
        plot_start = x_min + (total_width * start / ref_contig_length)
        plot_end = x_min + (total_width * end / ref_contig_length)
        colour_value = 255 - int(255 * depth)
        if binary_blue:
            if depth == 1:
                colours = [0, 0, 255]
            elif depth == 0:
                continue
            else:
                sys.exit('Unexpected error')
        else:
            colours = [colour_value] * 3

        colour = 'rgb(' + ','.join([str(x) for x in colours]) + ')'
        coords = [(plot_start, y_top), (plot_end, y_top), (plot_end, y_bottom), (plot_start, y_bottom)]
        print(svg_polygon(coords, colour, colour, border_width=0), file=filehandle)


def write_ref_contigs(ref_contigs, ref_coords, heatmap_top, heatmap_bottom, heatmap_middle, y_top, y_bottom, filehandle):
    total_height = y_bottom - y_top
    lineplot_height = 0.3 * total_height
    lineplot_space = 0.1 * total_height
    lineplot_top_top = y_top
    lineplot_top_bottom = lineplot_top_top + lineplot_height
    lineplot_middle_top = lineplot_top_bottom + lineplot_space
    lineplot_middle_bottom = lineplot_middle_top + lineplot_height
    lineplot_bottom_top = lineplot_middle_bottom + lineplot_space
    lineplot_bottom_bottom = lineplot_bottom_top + lineplot_height

    for ref_contig in ref_contigs:
        start, end = ref_coords[ref_contig.id]
        coords = [(start, y_top), (end, y_top), (end, y_bottom), (start, y_bottom)]
        #print(svg_polygon(coords, "white", "black", border_width=0.2), file=filehandle)
        write_lineplot(heatmap_top[ref_contig.id], lineplot_top_top, lineplot_top_bottom, start, end, len(ref_contig), filehandle, "black")
        write_lineplot(heatmap_middle[ref_contig.id], lineplot_middle_top, lineplot_middle_bottom, start, end, len(ref_contig), filehandle, "blue")
        write_lineplot(heatmap_bottom[ref_contig.id], lineplot_bottom_top, lineplot_bottom_bottom, start, end, len(ref_contig), filehandle, "black")


def load_ids_file(filename):
    if filename is None:
        return set()
    else:
        with open(filename) as f:
            ids = [x.rstrip() for x in f.readlines()]
        return ids


def write_assembly_contigs(contig_coords, y_bottom, y_top, filehandle, circular_contigs_file=None, merged_contigs_file=None, removed_contigs_file=None):
    circular_contigs = load_ids_file(circular_contigs_file)
    merged_contigs = load_ids_file(merged_contigs_file)
    removed_contigs = load_ids_file(removed_contigs_file)

    for (name, start, end) in contig_coords:
        coords = [(start, y_top), (end, y_top), (end, y_bottom), (start, y_bottom)]
        if name in circular_contigs:
            fill_colour = "lightgreen"
        elif name in merged_contigs and name in removed_contigs:
            fill_colour = "orange"
        elif name in removed_contigs:
            fill_colour = "red"
        elif name in merged_contigs:
            fill_colour = "yellow"
        else:
            fill_colour = "grey"
        print(svg_polygon(coords, fill_colour, fill_colour, border_width=0), file=filehandle)


def write_nucmer_hits(hit_coords, y_bottom, y_top, filehandle, y_flip=False):
    for qry_start, qry_end, ref_end, ref_start, same_strand in hit_coords:
        if y_flip:
            qry_start, ref_start = ref_start, qry_start
            qry_end, ref_end = ref_end, qry_end
        coords = [(qry_start, y_top), (qry_end, y_top), (ref_end, y_bottom), (ref_start, y_bottom)]
        if same_strand:
            border = "blue"
            fill = "lightblue"
        else:
            border = "orchid"
            fill = "pink"
        print(svg_polygon(coords, fill, border, border_width=0.1, opacity=0.5), file=filehandle)


parser = argparse.ArgumentParser(
    description = 'Makes cartoon figure comparing three genomes',
    usage = '%(prog)s <assembly_before.fasta> <assembly_after.fasta> <reference.fasta> <outprefix>')
parser.add_argument('--circular_contigs_file', help='Filename of names of circular contigs, one name per line')
parser.add_argument('--merged_contigs_file', help='Filename of names of merged contigs, one name per line')
parser.add_argument('--removed_contigs_file', help='Filename of names of contigs that were removed')
parser.add_argument('--min_id', type=int, help='Minimum identity when running nucmer [%(default)s]', default=98)
parser.add_argument('--lower_contigs_offset', type=float, help='Use to shift ref and final assembly contigs to the right', default=0)
parser.add_argument('assembly_before', help='Name of assembly before circularisation')
parser.add_argument('assembly_after', help='Name of assembly after circularisation')
parser.add_argument('ref_fasta', help='Reference genome')
parser.add_argument('outprefix', help='Prefix of output files')
options = parser.parse_args()

nucmer_before_v_ref = options.outprefix + '.nucmer_before_v_ref.coords'
nucmer_after_v_ref = options.outprefix + '.nucmer_after_v_ref.coords'

ref_contigs = assembly_to_list(options.ref_fasta)
contigs_before = {}
pyfastaq.tasks.file_to_dict(options.assembly_before, contigs_before)
contigs_after = {}
pyfastaq.tasks.file_to_dict(options.assembly_after, contigs_after)


ref_repeat_coords = get_repeat_coords(options.ref_fasta)
repeat_heatmap = make_ref_repeat_heatmap(ref_contigs, options.ref_fasta, min_id=options.min_id)
run_nucmer(options.ref_fasta, options.assembly_before, nucmer_before_v_ref, min_id=options.min_id)
run_nucmer(options.ref_fasta, options.assembly_after, nucmer_after_v_ref, min_id=options.min_id)
nucmer_hits_before_v_ref = load_nucmer_hits(nucmer_before_v_ref)
nucmer_hits_after_v_ref = load_nucmer_hits(nucmer_after_v_ref)
#os.unlink(nucmer_before_v_ref)
#os.unlink(nucmer_after_v_ref)

before_contig_coords, before_hit_coords, before_reversed_contigs, before_max_x = get_hit_and_contig_coords(nucmer_hits_before_v_ref, ref_contigs, contigs_before)
after_contig_coords, after_hit_coords, after_reversed_contigs, after_max_x = get_hit_and_contig_coords(nucmer_hits_after_v_ref, ref_contigs, contigs_after, lower_contigs_offset=options.lower_contigs_offset)

ref_contig_coords, ref_max_x = get_ref_contig_coords(ref_contigs, x_offset=options.lower_contigs_offset)
max_x = max(before_max_x, after_max_x, ref_max_x)
ref_heatmap_before = get_ref_heatmap_data(ref_contigs, nucmer_hits_before_v_ref)
ref_heatmap_after = get_ref_heatmap_data(ref_contigs, nucmer_hits_after_v_ref)
normalise_heatmaps(ref_heatmap_before, ref_heatmap_after)

svg_width=max_x
svg_height=max_x * 0.3
assembly_contig_height = svg_height * 0.04
ref_contig_height = 3 * assembly_contig_height
contig_to_hit_space = assembly_contig_height * 0.2

before_contig_top = 0
before_contig_bottom = assembly_contig_height
ref_contig_top = (0.5 * svg_height) - (0.5 * ref_contig_height)
ref_contig_bottom = (0.5 * svg_height) + (0.5 * ref_contig_height)
after_contig_top = svg_height - assembly_contig_height
after_contig_bottom = svg_height
before_v_ref_hits_top = before_contig_bottom + contig_to_hit_space
before_v_ref_hits_bottom = ref_contig_top - contig_to_hit_space
after_v_ref_hits_top = ref_contig_bottom + 4 * contig_to_hit_space
after_v_ref_hits_bottom = after_contig_top - contig_to_hit_space


svg_file = options.outprefix + '.svg'
f = pyfastaq.utils.open_file_write(svg_file)

print (r'''<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC " -//W3C//DTD SVG 1.0//EN" "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">
<svg width="''' + str(svg_width) + '" height="' + str(svg_height) + '">', file=f)


print('<!-- reference contigs -->', file=f)
write_ref_contigs(ref_contigs, ref_contig_coords, ref_heatmap_before, ref_heatmap_after, repeat_heatmap, ref_contig_top, ref_contig_bottom, f)
print('<!-- assembly contigs before -->', file=f)
write_assembly_contigs(before_contig_coords, before_contig_bottom, before_contig_top, f, merged_contigs_file=options.merged_contigs_file, removed_contigs_file=options.removed_contigs_file)
print('<!-- assembly contigs after -->', file=f)
write_assembly_contigs(after_contig_coords, after_contig_bottom, after_contig_top, f, circular_contigs_file=options.circular_contigs_file)
print('<!-- hits before vs ref -->', file=f)
write_nucmer_hits(before_hit_coords, before_v_ref_hits_bottom, before_v_ref_hits_top, f)
print('<!-- hits after vs ref -->', file=f)
write_nucmer_hits(after_hit_coords, after_v_ref_hits_bottom, after_v_ref_hits_top, f, y_flip=True)

print('</svg>', file=f)
pyfastaq.utils.close(f)

