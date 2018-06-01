
'''
This script creates json files which can be used to render Manhattan plots.
'''

from ..utils import chrom_order
from ..conf_utils import conf
from ..file_utils import VariantFileReader, write_json, common_filepaths
from .load_utils import MaxPriorityQueue, parallelize_per_pheno

import math


def run(argv):
    parallelize_per_pheno(
        get_input_filepaths = lambda pheno: common_filepaths['pheno'](pheno['phenocode']),
        get_output_filepaths = lambda pheno: common_filepaths['manhattan'](pheno['phenocode']),
        convert = make_json_file,
        cmd = 'manhattan',
    )


def make_json_file(pheno):
    BIN_LENGTH = int(3e6)
    NEGLOG10_PVAL_BIN_SIZE = 0.05 # Use 0.05, 0.1, 0.15, etc
    NEGLOG10_PVAL_BIN_DIGITS = 2 # Then round to this many digits

    with VariantFileReader(common_filepaths['pheno'](pheno['phenocode'])) as variants:
        variant_bins, unbinned_variants = bin_variants(
            variants,
            BIN_LENGTH,
            NEGLOG10_PVAL_BIN_SIZE,
            NEGLOG10_PVAL_BIN_DIGITS
        )
    label_peaks(unbinned_variants)
    rv = {
        'variant_bins': variant_bins,
        'unbinned_variants': unbinned_variants,
    }

    write_json(filepath=common_filepaths['manhattan'](pheno['phenocode']), data=rv)


def rounded_neglog10(pval, neglog10_pval_bin_size, neglog10_pval_bin_digits):
    return round(-math.log10(pval) // neglog10_pval_bin_size * neglog10_pval_bin_size, neglog10_pval_bin_digits)


def get_pvals_and_pval_extents(pvals, neglog10_pval_bin_size):
    # expects that NEGLOG10_PVAL_BIN_SIZE is the distance between adjacent bins.
    pvals = sorted(pvals)
    extents = [[pvals[0], pvals[0]]]
    for p in pvals:
        if extents[-1][1] + neglog10_pval_bin_size * 1.1 > p:
            extents[-1][1] = p
        else:
            extents.append([p,p])
    rv_pvals, rv_pval_extents = [], []
    for (start, end) in extents:
        if start == end:
            rv_pvals.append(start)
        else:
            rv_pval_extents.append([start,end])
    return (rv_pvals, rv_pval_extents)


# TODO?: unbin the top few variant in each peak, to avoid the top peak monopolizing all the unbinned variants
def bin_variants(variant_iterator, bin_length, neglog10_pval_bin_size, neglog10_pval_bin_digits):
    peak_best_variant = None
    peak_last_chrpos = None
    bins = {} # like {<chrom>: {<pos // bin_length>: [{chrom, startpos, neglog10pvals}]}}
    unbinned_variant_pq = MaxPriorityQueue()
    peak_pq = MaxPriorityQueue()

    def bin_variant(variant):
        chrom_idx = chrom_order[variant['chrom']]
        if chrom_idx not in bins: bins[chrom_idx] = {}
        pos_bin_id = variant['pos'] // bin_length
        if pos_bin_id not in bins[chrom_idx]:
            bins[chrom_idx][pos_bin_id] = {'chrom': variant['chrom'], 'startpos': pos_bin * bin_length, 'neglog10_pvals': set()}
        neglog10_pval = rounded_neglog10(variant['pval'], neglog10_pval_bin_size, neglog10_pval_bin_digits)
        bins[chrom_idx][pos_bin_id]["neglog10_pvals"].add(neglog10_pval)

    def maybe_bin_variant(variant):
        unbinned_variant_pq.add(variant, variant['pval'])
        if len(unbinned_variant_pq) > conf.manhattan_num_unbinned:
            old = unbinned_variant_pq.pop()
            bin_variant(old)

    def maybe_peak_variant(variant):
        peak_pq.add(variant)
        if len(peak_pq) > conf.manhattan_peak_max_count:
            old = peak_pq.pop()
            maybe_bin_variant(old)

    for variant in variant_iterator:
        if variant['pval'] < conf.manhattan_peak_pval_threshold: # part of a peak
            if peak_best_variant is None: # open a new peak
                peak_best_variant = variant
                peak_last_chrpos = (variant['chrom'], variant['pos'])
            elif peak_last_chrpos[0] == variant['chrom'] and peak_last_chrpos[1] + conf.manhattan_peak_sprawl_dist > variant['pos']: # extend current peak
                peak_last_chrpos = (variant['chrom'], variant['pos'])
                if variant['pval'] < peak_best_variant['pval']:
                    maybe_bin_variant(peak_best_variant)
                    peak_best_variant = variant
            else: # close old peak and open new peak
                maybe_peak_variant(peak_best_variant)
                peak_best_variant = variant
                peak_last_chrpos = (variant['chrom'], variant['pos'])
        else:
            maybe_bin_variant(variant)
    maybe_peak_variant(peak_best_variant)

    peaks = list(peak_pq.pop_all())
    for peak in peaks: peak['peak'] = True
    unbinned_variants = list(unbinned_variant_pq.pop_all())
    unbinned_variants = sorted(unbinned_variants + peaks, key=(lambda variant: variant['pval']))

    # unroll dict-of-dict-of-array `bins` into array `variant_bins`
    variant_bins = []
    for chrom_idx in sorted(bins.keys()):
        for pos_bin_id in sorted(bins[chrom_idx].keys()):
            b = bins[chrom_idx][pos_bin_id]
            assert len(b['neglog10_pvals']) > 0
            b['neglog10_pvals'], b['neglog10_pval_extents'] = get_pvals_and_pval_extents(b['neglog10_pvals'], neglog10_pval_bin_size)
            b['pos'] = int(b['startpos'] + bin_length/2)
            del b['startpos']
            variant_bins.append(b)

    return variant_bins, unbinned_variants


def bin_variants(variant_iterator, bin_length, neglog10_pval_bin_size, neglog10_pval_bin_digits):
    bins = {} # like {(chrom_key, pos//bin_length): [...]}
    unbinned_variant_pq = MaxPriorityQueue()
    chrom_n_bins = {}

    def bin_variant(variant):
        chrom_key = chrom_order[variant['chrom']]
        pos_bin = variant['pos'] // bin_length
        chrom_n_bins[chrom_key] = max(chrom_n_bins.get(chrom_key,0), pos_bin)
        if (chrom_key, pos_bin) in bins:
            bin = bins[(chrom_key, pos_bin)]

        else:
            bin = {"chrom": variant['chrom'],
                   "startpos": pos_bin * bin_length,
                   "neglog10_pvals": set()}
            bins[(chrom_key, pos_bin)] = bin
        bin["neglog10_pvals"].add(rounded_neglog10(variant['pval'], neglog10_pval_bin_size, neglog10_pval_bin_digits))

    # put most-significant variants into the priorityqueue and bin the rest
    for variant in variant_iterator:
        unbinned_variant_pq.add(variant, variant['pval'])
        if len(unbinned_variant_pq) > conf.manhattan_num_unbinned:
            old = unbinned_variant_pq.pop()
            bin_variant(old)

    unbinned_variants = list(unbinned_variant_pq.pop_all())

    # unroll bins into simple array (preserving chromosomal order)
    binned_variants = []
    for chrom_key in sorted(chrom_n_bins.keys()):
        for pos_key in range(int(1+chrom_n_bins[chrom_key])):
            b = bins.get((chrom_key, pos_key), None)
            if b and len(b['neglog10_pvals']) != 0:
                b['neglog10_pvals'], b['neglog10_pval_extents'] = get_pvals_and_pval_extents(b['neglog10_pvals'], neglog10_pval_bin_size)
                b['pos'] = int(b['startpos'] + bin_length/2)
                del b['startpos']
                binned_variants.append(b)

    return binned_variants, unbinned_variants

def label_peaks(variants):
    chroms = {}
    for v in variants:
        chroms.setdefault(v['chrom'], []).append(v)
    for vs in chroms.values():
        while vs:
            best_assoc = min(vs, key=lambda assoc: assoc['pval'])
            best_assoc['peak'] = True
            vs = [v for v in vs if abs(v['pos'] - best_assoc['pos']) > conf.within_pheno_mask_around_peak]
