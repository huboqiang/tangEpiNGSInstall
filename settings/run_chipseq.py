#-*- coding:utf-8 -*-
from __future__ import division
import re,sys,os

import cPickle as pickle
import numpy as np
from optparse   import OptionParser

import ChIP.frame.module01_mapping_from_raw as m01
import ChIP.frame.module02_call_peaks       as m02
import ChIP.frame.module03_geneDensity      as m03
import ChIP.utils.module_create_database    as m_db

def prepare_optparser():
    usage ="""usage: %s [options]
Reference fasta file and refGene file should be all put in dictionary
self.Database, which were defined in settings/projpath.py
If not put in self.Database, this program will downloading from UCSC.
Detail information could be reached in utils/module_create_database.py
Suport genome includes:
    http://hgdownload.soe.ucsc.edu/goldenPath
For species of inputs, if ref is hg19/hg38 or mm9/mm10, the default reference
for macs2 is hs or mm. Otherwise, please using the genome-size total length in
sam_peak.run_macs_rep, sam_peak.run_macs_mrg, sam_peak.run_macs_rep_broad or
sam_peak.run_macs_mrg_broad instead of M_species[ref].
Using -h or --help for more information
Example:
    python %s --ref hg19 --TSS_genebody_up 5000 --TSS_genebody_down 5000 --TSS_promoter_up 5000 --TSS_promoter_down 5000 --Body_extbin_len 50 --Body_bincnt 100 --TSS_bin_len 1 --top_peak_idr 100000 sample.20150719.xls
    """ % (sys.argv[0],sys.argv[0])

    description = " The ChIP seq analysis standard pipeline. "

    optparser = OptionParser(version="%s v0.2 20141130" % (sys.argv[0]),
        description=description,
        usage=usage,
        add_help_option=False
    )
    optparser.add_option(
        "-r", "--ref", default="mm9",
        help="\nReference genome. [default: %default]"
    )
    optparser.add_option("--TSS_genebody_up"  ,
        default=5000,
        help="\nTSS upstream   from a genebody [default: %default]"
    )
    optparser.add_option("--TSS_genebody_down",
        default=5000,
        help="\nTES downstream from a genebody [default: %default]"
    )
    optparser.add_option("--TSS_promoter_up"  ,
        default=5000,
        help="\nTSS upstream                   [default: %default]"
    )
    optparser.add_option("--TSS_promoter_down",
        default=5000,
        help="\nTSS downstream                 [default: %default]"
    )
    optparser.add_option("--Body_extbin_len"  ,
        default=50,
        help="\nGenebody extending length for a bin   [default: %default]"
    )
    optparser.add_option("--Body_bincnt"      ,
        default=100,
        help="\nGenebody divided bins          [default: %default]"
    )
    optparser.add_option("--TSS_bin_len"      ,
        default=1,
        help="\nTSS region bin-length          [default: %default]"
    )
    optparser.add_option("--top_peak_idr"     ,
        default=100000,
        help="\nUsing how many peaks for IDR analysis [default: %default]"
    )
    optparser.add_option("-h","--help", action="help",
        help="\nShow this help message and exit."
    )
    return optparser

M_species = {
    'hg19' : "hs",
    'hg38' : "hs",
    'mm9'  : "mm",
    'mm10' : "mm",
}

def main():
    prepare_optparser()
    (options,args) = prepare_optparser().parse_args()
    try:
        sam_file_chip = args[0]
        ref           = options.ref

        TSS_genebody_up     = int(options.TSS_genebody_up  )
        TSS_genebody_down   = int(options.TSS_genebody_down)
        TSS_promoter_up     = int(options.TSS_promoter_up  )
        TSS_promoter_down   = int(options.TSS_promoter_down)

        Body_extbin_len     = int(options.Body_extbin_len  )
        Body_bincnt         = int(options.Body_bincnt      )
        TSS_bin_len         = int(options.TSS_bin_len      )

        top_peak_idr        = int(options.top_peak_idr     )

    except IndexError:
        prepare_optparser().print_help()
        sys.exit(1)


    part0 = m_db.DataBaseInit(ref, sam_file_chip, is_debug=0)
    part0.check_files()

    samp_mapping = m01.Map_From_raw(sam_file_chip, ref, is_debug = 0)
    samp_mapping.run_QC(core_num=4)
    samp_mapping.run_bwa(core_num=1)
    samp_mapping.bam2repbed(ext_len=300, core_num=4)
    samp_mapping.mrgbed(core_num=4)
    samp_mapping.mrgbed_multi(core_num=4)

    samp_peak = m02.Macs2Peaks(sam_file_chip, ref, is_debug = 0)
    samp_peak.get_shift_size_rep(core_num=4)
    samp_peak.get_shift_size_mrg(core_num=4)
    samp_peak.run_macs_rep(pvalue=0.001, ref=M_species[ref], core_num=4)
    samp_peak.run_macs_mrg(pvalue=0.001, ref=M_species[ref], core_num=4)
    samp_peak.prepare_idr_input(top_peak = 100000, core_num=4)
    samp_peak.run_idr(core_num=4)
    samp_peak.get_idr_stat()
    samp_peak.get_idr_Peak(core_num=4)
    samp_peak.run_macs_rep_broad(pvalue=0.05, ref=M_species[ref], core_num=4)
    samp_peak.run_macs_mrg_broad(pvalue=0.05, ref=M_species[ref], core_num=4)
    samp_peak.sort_bdg(core_num=4)
#    samp_peak.make_igv_IDR()
#    samp_peak.make_igv_broad()

    peak_gene = m03.GeneBasedInfo(sam_file_chip,ref, is_debug = 0)
    #peak_gene.extend_gene_region(
    #    TSS_genebody_up,TSS_genebody_down,TSS_promoter_up,TSS_promoter_down
    #)
    #peak_gene.run_anno_peak(
    #    TSS_genebody_up,TSS_genebody_down,TSS_promoter_up,TSS_promoter_down,
    #    Body_extbin_len, Body_bincnt, TSS_bin_len
    #)
    #peak_gene.div_bed_to_bins_rep()
    #peak_gene.div_bed_to_bins_mrg()
    #peak_gene.merge_RPKM()
    peak_gene.density_baseLevel(core_num=4)

if __name__ == '__main__':
    main()