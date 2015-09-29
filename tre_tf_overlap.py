#! /usr/bin/env python

# system
import os
import json
import argparse
import subprocess
import numpy as np

# graphing
import plotly.plotly as py
from plotly.graph_objs import *

def begin(args):
    # instantiate class
    run = TreTfOverlap(args)

    # execute
    run.get_gene_list()

class TreTfOverlap:

    def __init__(self, args):
        self.gene_list = args.gene_list

        # paths
        self.inpath = '/home/zachary/Zack/corelab/TRE_ChIP_analysis/'
        self.cwd = '/home/zachary/Zack/corelab/TRE_ChIP_analysis/data/'
        self.outpath = '/home/zachary/Zack/corelab/TRE_ChIP_analysis/data/out/'

    def get_gene_list(self):

        # Set which gene list to process.
        if self.gene_list == 'all':
            _up_ = 'TRE_order.bed'
            _down_ = None
            inpath = self.cwd

            # Bed file has to be rearranged somewhat.
            _up_ = self.smplify_bed(_up_)

        elif self.gene_list == 'reg':
            _up_ = 'TRE_upReg.txt'
            _down_ = 'TRE_downReg.txt'
            inpath = self.cwd

            # These files must be converted to bed and sorted before they can be processed.
            _up_, _down_ = self.txt_to_bed(_up_, _down_)

        elif self.gene_list == 'sig':
            _up_ = "all_RegTRE_upGenes_sorted.txt"
            _down_ = "all_RegTRE_downGenes_sorted.txt"
            inpath = self.inpath

        # execute #
        outpath = self.mkdir()
        self.execute_sig(_up_, _down_, inpath, outpath)

    def execute_sig(self, _up_, _down_, path, outpath, lbl=None):

        # ChIP-seq peaks for various transcription factors (TFs).

        # Set Up.
        inpath = self.inpath
        cwd = self.cwd
        tf_list = ["cMyc_mm9.bed", "CTCF_mm9.bed", "E2F1_mm9.bed", "Esrrb_mm9.bed", "Klf4_mm9.bed", "Nanog_mm9.bed",
                 "nMyc_mm9.bed", "Oct4_mm9.bed", "p300_mm9.bed", "Sox2_mm9.bed", "Stat3_mm9.bed", "Suz12_mm9.bed",
                 "Tcfcp2l1_mm9.bed", "zfx_mm9.bed"]
        up = '%s%s' % (path, _up_)
        down = '%s%s' % (path, _down_)
        tf_tre_overlap_up = open('%soverlap_up.txt' % outpath, 'w')
        tf_tre_overlap_down = open('%soverlap_down.txt' % outpath, 'w')
        plot_overlap_up = []
        plot_overlap_down = []

        # Total number of TREs.
        i = 0
        j = 0
        with open(up, 'r') as f:
            for line in f:
                i += 1

        if os.path.isfile(down):
            with open(down, 'r') as ff:
                for lline in ff:
                    j += 1

        total_up_tres = i
        total_down_tres = j

        if _down_:
            ls = [tf_tre_overlap_up, tf_tre_overlap_down]
        else:
            ls = [tf_tre_overlap_up]
            lbl = 'all'

        # Create table and graphs
        for _tre_ in ls:

            header = ''
            overlap = ''
            gene = ''
            promoter = ''
            intergenic = ''
            overlap = ''
            _gene_ = ''
            _promoter_ = ''
            _intergenic_ = ''
            plot_tf = []
            _plot_gene = []
            _plot_promoter = []
            _plot_intergenic = []

            for tf_file in tf_list:

                # set up
                tf = '%s%s' % (inpath, tf_file)
                _tf_ = tf_file.split('_')[0]
                plot_tf.append(_tf_)
                tf_out = '%s%s_extended.bed' % (cwd, _tf_)
                up_overlap = '%s%s_up_overlap.txt' % (cwd, _tf_)
                down_overlap = '%s%s_down_overlap.txt' % (cwd, _tf_)

                # create header for the overlap table.
                header += '%s\t' % _tf_

                # Extend the TF peaks 150bp on each side.
                # ----------------------------------------

                # load
                chr_ = np.loadtxt(tf, dtype=str, usecols=(0,))
                start = np.loadtxt(tf, dtype=int, usecols=(1,))
                stop = np.loadtxt(tf, dtype=int, usecols=(2,))

                # execute
                for k, v in enumerate(start):
                    start[k] = v - 150

                for k, v in enumerate(stop):
                    stop[k] = v + 150

                # write
                i = 0
                with open(tf_out, 'w') as f:
                    while i < len(chr_):
                        f.write('%s\t%s\t%s\n' % (chr_[i], start[i], stop[i]))
                        i += 1

                # Record the overlap between each transcription factor and the GRO-seq TREs.
                # --------------------------------------------------------------------------

                # set up
                if _tre_ == tf_tre_overlap_up:
                    _in_ = _up_
                    _out_write = open(up_overlap, 'a')
                    _out_bedtools = open(up_overlap, 'w')
                    _out_ = up_overlap
                    _total_ = total_up_tres
                    plot_overlap = plot_overlap_up
                    tf_tre_overlap = tf_tre_overlap_up
                    if not lbl:
                        lbl = 'up'
                elif _tre_ == tf_tre_overlap_down:
                    _in_ = _down_
                    _out_write = open(down_overlap, 'a')
                    _out_bedtools = open(down_overlap, 'w')
                    _out_ = down_overlap
                    _total_ = total_down_tres
                    plot_overlap = plot_overlap_down
                    tf_tre_overlap = tf_tre_overlap_down
                    lbl = 'down'

                else:
                    _in_ = _out_ = _total_ = _out_write = _out_bedtools = \
                        _plot_promoter_ = _plot_gene_ = _plot_intergenic_ = None
                    raise ValueError('Cannot open %s or %s' % (_up_, _down_))

                # overlap with bedtools intersect
                cmd = ["bedtools", "intersect", "-wa", "-a", "%s%s" % (path, _in_), "-b", tf_out]
                subprocess.call(cmd, stdout=_out_bedtools)

                # Determine the fraction of TREs that have overlapped with this TF and write to the overlap table.
                # ------------------------------------------------------------------------------------------------

                # read output
                i = 0
                tre_list = []
                region_list = []

                with open(_out_, 'r') as f:
                    for line in f:
                        try:
                            col = line.strip().split('\t')
                            tre = col[3]        # TRE ID
                            region = col[4]     # The genomic region in which that TRE resides
                            tre_list.append(tre)
                            region_list.append(region)
                        except IndexError:
                            print '%s\n%s' % (_tf_, line)
                            raise
                        i += 1

                # calculate overlap
                per_overlap = float(i) / float(_total_)
                overlap += '%s\t' % per_overlap

                # add to list for graph
                plot_overlap.append(per_overlap)

                # determine gene region fraction
                gene = 0
                promoter = 0
                intergenic = 0
                for r in region_list:
                    if r == 'gene':
                        gene += 1
                    elif r == 'promoter':
                        promoter += 1
                    elif r == 'intergenic':
                        intergenic += 1

                # calculate
                num_gene = float(gene) / float(i)
                num_promoter = float(promoter) / float(i)
                num_intergenic = float(intergenic) / float(i)

                # Add text
                _gene_ += '%s\t' % num_gene
                _promoter_ += '%s\t' % num_promoter
                _intergenic_ += '%s\t' % num_intergenic

                # Add to list for graph
                _plot_gene.append(num_gene)
                _plot_promoter.append(num_promoter)
                _plot_intergenic.append(num_intergenic)

                # Record to table.
                assert len(tre_list) == len(region_list), 'Some of the tres are missing their assigned regions.'
                i = 0
                _out_write.write('%s\n' % _tf_)  # TF header
                while i < len(tre_list):
                    try:
                        _out_write.write('%s\t%s\n' % (tre_list[i], region_list[i]))  # TREs
                        i += 1
                    except IndexError:
                        self.print_json(tre_list)
                        self.print_json(region_list)
                        raise

            # Write table and plot gene region graphs for up and down regulated genes.
            self.write_table(tf_tre_overlap, header, overlap, _gene_, _promoter_, _intergenic_)
            self.plot_genome_regions(plot_tf, _plot_gene, _plot_promoter, _plot_intergenic, lbl, outpath)

        # Plot results
        self.summary_plot(plot_tf, plot_overlap_up, plot_overlap_down, outpath)

    def print_json(self, msg):
        print json.dumps(msg, sort_keys=True, indent=4)

    def write_table(self, tre, header, overlap, _gene_, _promoter_, _intergenic_):
        tre.write('TF:\t%s\n'
                 'Percent_Overlap:\t%s\n'
                 'Percent_of_Overlap_in_Gene:\t%s\n'
                 'Percent_of_Overlap_in_Promoter:\t%s\n'
                 'Percent_of_Overlap_in_Intergenic:\t%s' % (header, overlap, _gene_, _promoter_, _intergenic_))

    # Graph the results.
    # -------------------
    def plot_genome_regions(self, plot_tf, _plot_gene_, _plot_promoter_, _plot_intergenic_, lbl, outpath):
        # Gene and Promoter regions
        plot_gene = Bar(
            x=plot_tf,
            y=_plot_gene_,
            name='Percent of TF-TRE overlap within gene regions'
        )

        plot_promoter = Bar(
            x=plot_tf,
            y=_plot_promoter_,
            name='Percent of TF-TRE overlap within promoter regions'
        )

        plot_intergenic = Bar(
            x=plot_tf,
            y=_plot_intergenic_,
            name='Percent of TF-TRE overlap within intergenic regions'
        )

        data = Data([plot_gene, plot_promoter, plot_intergenic])
        layout = Layout(barmode='stack', title='Genomic regions of TF-TRE overlaps for %s regulated genes' % lbl)
        fig = Figure(data=data, layout=layout)
        py.image.save_as(fig, filename='%sTRE_TF_overlap_by_genomic_region_%s.pdf' % (outpath, lbl))

    def summary_plot(self, plot_tf, plot_overlap_up, plot_overlap_down, outpath):
        # Overall overlap
        up = Bar(
            x=plot_tf,
            y=plot_overlap_up,
            name='Up regulated genes'
        )

        down = Bar(
            x=plot_tf,
            y=plot_overlap_down,
            name='Down regulated genes'
        )

        data = Data([up, down])
        layout = Layout(barmode='group', title='Percentage of TRE overlap with various TFs',
                        yaxis=YAxis(
                            title='Percent of TREs that overlap the TF peak'
                        ))
        fig = Figure(data=data, layout=layout)
        py.image.save_as(fig, filename='%sTRE_TF_overlap.pdf' % outpath)

    def txt_to_bed(self, up, down):

        # Create a bed file from the txt file.
        new = []
        inpath = self.inpath
        for reg_file in [up, down]:

            outfile = self.cwd + reg_file.replace('.txt', '.bed')

            with open(inpath + reg_file, 'r') as f:
                with open(outfile, 'w') as ff:
                    for line in f:
                        if line.startswith('TRE'):
                            col = line.strip().split('\t')

                            msg = [col[5], col[6], col[7], col[0], col[13], '+\n']
                            ff.write(str.join('\t', msg))

            # sort
            outfile_sorted = outfile.replace('.bed', '_sorted.bed')
            sorted_stdout = open(outfile_sorted, 'w')
            return_out = outfile_sorted.split('/')[-1]

            cmd = ["sort", "-k4", "-n", outfile]
            subprocess.call(cmd, stdout=sorted_stdout)

            new.append(return_out)

        return new[0], new[1]

    def smplify_bed(self, reg_file):

        # Create a bed file from the txt file.
        inpath = self.inpath
        outfile = reg_file.replace('.bed', '_rdy.bed')

        with open(inpath + reg_file, 'r') as f:
            with open(self.cwd + outfile, 'w') as ff:
                for line in f:
                    col = line.strip().split('\t')

                    msg = [col[0], col[1], col[2], col[3], col[6], '%s\n' % col[5]]
                    ff.write(str.join('\t', msg))

        return outfile

    def mkdir(self):

        if self.gene_list == 'all':
            dir_ = 'all/'
        elif self.gene_list == 'reg':
            dir_ = 'reg/'
        elif self.gene_list == 'sig':
            dir_ = 'sig/'

        dir_ = self.outpath + dir_

        if not os.path.isdir(dir_):
            cmd = ["mkdir", dir_]
            subprocess.call(cmd)

        return dir_

if __name__ == '__main__':

    # mode parser.
    parser = argparse.ArgumentParser()

    # get gene lists
    parser.add_argument('--gene-list', choices=['all', 'reg', 'sig'], default='sig',
                      help='all: signifies all 80,000 TREs\n'
                           'reg: signifies the up and down regulated TREs\n'
                           'sig: signfies the significant up and down regulated TREs')
    parser.set_defaults(func=begin)

    # parse args.
    args = parser.parse_args()
    args.func(args)


