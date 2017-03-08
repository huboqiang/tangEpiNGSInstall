#-*- coding:utf-8 -*-
from __future__ import division
import re
import sys
import os
import subprocess
import cPickle as pickle
import time

import ChIP.utils.module_running_jobs as m_jobs
import ChIP.settings.projpath         as m_proj

class Scripts(m_proj.ProjInfo):
    def __init__(self):
        super(Scripts,self).__init__()


    def define_files(self,ref):
        dir_db = "%s/%s" % (self.Database, ref)
        self.ref           = ref
        self.genome_ref    = "%s/%s.fa"       % (dir_db, ref)
        self.refGeneTxt    = "%s/refGene.txt" % (dir_db)
        self.genome_gtf    = "%s/refGene.gtf" % (dir_db)
        self.rmsk_bed      = "%s/chrom.sort.bed"        % (dir_db)

    def db_01_DownloadRef(self):
        l_sh_info = []
        l_sh_info.append("ref=$1")
        l_sh_info.append("dir_database=%s/$ref" % (self.Database))
        l_sh_info.append("dir_path=%s"          % (self.path))
        l_sh_info.append("""
cd $dir_database

wget http://hgdownload.soe.ucsc.edu/goldenPath/${ref}/bigZips/chromFa.tar.gz

tar -zxvf $dir_database/chromFa.tar.gz

for i in {1..22} X Y M
do
    cat $dir_database/chr$i.fa
done  >$dir_database/${ref}.fa && rm $dir_database/chr*fa
        """)
        return l_sh_info

    def db_02_BuildRefIndex(self):
        getbin_cpp = "%s/div_bins/bed_read"% (self.bin)
        l_sh_info = []
        l_sh_info.append("ref=$1")
        l_sh_info.append("dir_database=%s/$ref" % (self.Database))
        l_sh_info.append("bwa_exe=%s"   % (self.sftw_bwa))
        l_sh_info.append("samtools_exe=%s" % (self.sftw_samtools))
        l_sh_info.append("div_bins_exe=%s" % (getbin_cpp))
        l_sh_info.append("""
$samtools_exe faidx $dir_database/${ref}.fa

$bwa_exe index $dir_database/${ref}.fa

$dix_bins_exe -b 100  $dir_database/${ref}.fa.fai $dir_database/columns.100.bed
$dix_bins_exe -b 1000 $dir_database/${ref}.fa.fai $dir_database/columns.1kb.bed

cut -f 1-2 $dir_database/${ref}.fa.fai >$dir_database/${ref}.fa.len
        """)
        return l_sh_info

    def db_03_RefGene(self):
        l_sh_info = []
        l_sh_info.append("ref=$1")
        l_sh_info.append("dir_database=%s/$ref" % (self.Database))
        l_sh_info.append("bedtools_exe=%s" % (self.sftw_bedtools))
        l_sh_info.append("ucsc_dir=%s"     % (self.sftw_ucsc_dir))
        l_sh_info.append("bin=%s"          % (self.bin))
        l_sh_info.append("dir_path=%s"     % (self.path))
        l_sh_info.append("""
cd $dir_database
wget http://hgdownload.soe.ucsc.edu/goldenPath/${ref}/database/refGene.txt.gz

### remove chromosome fragments(unassembled).
for i in {1..22} X Y M
do
    zcat $dir_database/refGene.txt.gz | grep -w chr$i
done >$dir_database/tmp
mv $dir_database/tmp    $dir_database/refGene.txt

# refGene.bed
cat $dir_database/refGene.txt                                              |\\
awk '{
    tag="noncoding";
    if($4~/^NM/){tag="protein_coding"};
    OFS="\\t";
    print $3,$5,$6,$2,$4,$10,$11,tag,$13
}' /dev/stdin                                                              |\\
python $bin/s03_genePred2bed.py /dev/stdin                                 |\\
$bedtools_exe sort -i /dev/stdin >$dir_database/refGene.bed               &&\\

# region.Intragenic.bed
# For novo lncRNA detection
$bin/find_ExonIntronIntergenic/find_ExonIntronIntergenic                    \\
    $dir_database/refGene.bed                                               \\
    $dir_database/${ref}.fa.fai >$dir_database/pos.bed                    &&\\

grep -v "Intergenic" $dir_database/pos.bed                                 |\\
    awk '{OFS="\t";print $1,$2,$3,"Intragenic"}' /dev/stdin                 \\
    >$dir_database/region.Intragenic.bed                                  &&\\

# refGene.gtf
# For mapping
zcat $dir_database/refGene.txt.gz                                          |\\
cut -f 2-                                                                  |\\
$ucsc_dir/genePredToGtf file stdin /dev/stdout                             |\\
grep -w exon                                                               |\\
$bedtools_exe sort -i /dev/stdin >$dir_database/refGene.gtf               &&\\
cat $dir_path/database/ERCC.gtf >>$dir_database/refGene.gtf
        """)
        return l_sh_info

    def db_04_rmsk(self):
        l_sh_info = []
        l_sh_info.append("ref=$1")
        l_sh_info.append("dir_database=%s/$ref" % (self.Database))
        l_sh_info.append("bedtools_exe=%s" % (self.sftw_bedtools))
        l_sh_info.append("ucsc_dir=%s"     % (self.sftw_ucsc_dir))
        l_sh_info.append("bin=%s"          % (self.bin))
        l_sh_info.append("dir_path=%s"     % (self.path))
        l_sh_info.append("""
cd $dir_database
wget http://hgdownload.soe.ucsc.edu/goldenPath/${ref}/database/rmsk.txt.gz

zcat $dir_database/rmsk.txt.gz                                             |\\
awk '{
    OFS="\\t";
    print $6,$7,$8,$2,".",".",".","("$9")",$10,$11,$12 "/" $13,$14,$15,$16,$17
}' /dev/stdin                                                              |\\
tail -n +2  /dev/stdin >$dir_database/chrom.bed

for i in {1..22} X Y M
do
    grep -w chr$i $dir_database/chrom.bed
done >$dir_database/tmp
mv $dir_database/tmp $dir_database/chrom.bed

$bedtools_exe sort -i   $dir_database/chrom.bed >$dir_database/chrom.sort.bed
        """)
        return l_sh_info

    def s01_QC(self):
        l_sh_info     = []
        l_sh_info.append("samp=$1")
        l_sh_info.append("data_type=$2")
        l_sh_info.append("pl_exe=%s"       % (self.sftw_pl) )
        l_sh_info.append("pl_QC=%s/QC.pl"  % (self.bin)     )
        l_sh_info.append("in_dir=%s"       % (self.dir_raw_data  ) )
        l_sh_info.append("out_dir=%s"      % (self.dir_clean_data) )
        l_sh_info.append("""""")
        l_sh_info.append("""
$pl_exe $pl_QC --indir $in_dir --outdir $out_dir --sample $samp             \\
    --end $data_type""")

        return l_sh_info

    def s02_bwa(self):
        l_sh_info = []
        l_sh_info.append("sam_name=$1")
        l_sh_info.append("brief_name=$2")
        l_sh_info.append("end=$3")
        l_sh_info.append("bwa_exe=%s"      % (self.sftw_bwa))
        l_sh_info.append("cln_dir=%s"      % (self.dir_clean_data))
        l_sh_info.append("bam_dir=%s"      % (self.dir_bam))
        l_sh_info.append("samtools_exe=%s" % (self.sftw_samtools))
        l_sh_info.append("genome=%s"       % (self.genome_ref))
        l_sh_info.append("""
bam_prefix=$bam_dir/$brief_name

if [ $end == "1" ]
  then $bwa_exe aln -i 15 -q 10 -t 4 $genome                                \\
          $cln_dir/$sam_name/1.cln.fq.gz >${bam_prefix}/1.sai            && \\
      $bwa_exe samse   $genome                                              \\
          ${bam_prefix}/1.sai $cln_dir/$sam_name/1.cln.fq.gz                \\
         >${bam_prefix}/${brief_name}.sam                                && \\
      $samtools_exe view -u -b -S -t $genome.fai                            \\
          ${bam_prefix}/${brief_name}.sam                                 | \\
      $samtools_exe sort -m 200000000 - ${bam_prefix}/${brief_name}      && \\
      rm ${bam_prefix}/1.sai && rm ${bam_prefix}/${brief_name}.sam
fi

if [ $end == "2" ]
  then $bwa_exe aln -i 15 -q 10 -t 4 $genome                                \\
          $cln_dir/$sam_name/1.cln.fq.gz >${bam_prefix}/1.sai            && \\
      $bwa_exe aln -i 15 -q 10 -t 4 $genome                                 \\
          $cln_dir/$sam_name/1.cln.fq.gz >${bam_prefix}/2.sai            && \\
      $bwa_exe sampe   $genome                                              \\
           ${bam_prefix}/1.sai            ${bam_prefix}/2.sai               \\
           $cln_dir/$sam_name/1.cln.fq.gz $cln_dir/$sam_name/2.cln.fq.gz    \\
          >${bam_prefix}/${brief_name}.sam                               && \\
      $samtools_exe view -u -b -S -t $genome.fai                            \\
           ${bam_prefix}/${brief_name}.sam                                | \\
      $samtools_exe sort -m 200000000 - ${bam_prefix}/${brief_name}      && \\
      rm ${bam_prefix}/*sai ${bam_prefix}/${brief_name}.sam
fi""")
        return l_sh_info

    def s03_bam2bedrep(self):
        l_sh_info = []
        l_sh_info.append("brief_name=$1")
        l_sh_info.append("end=$2")
        l_sh_info.append("ext_len=$3")
        l_sh_info.append("java_exe=%s"    % self.sftw_java    )
        l_sh_info.append("markDup_jar=%s" % self.sftw_MarkDup )
        l_sh_info.append("python_exe=%s"  % self.sftw_py      )
        l_sh_info.append("py_ExtRead=%s/ExtRead.py" % self.bin)
        l_sh_info.append("bam=%s"         % self.dir_bam      )
        l_sh_info.append("bed_rep=%s"     % self.dir_bed_rep  )
        l_sh_info.append("samtools_exe=%s"% self.sftw_samtools)
        l_sh_info.append("bedtools_exe=%s"% self.sftw_bedtools)
        l_sh_info.append("""
$java_exe -Xmx1g -jar $markDup_jar                                          \\
    INPUT=$bam/$brief_name/$brief_name.bam                                  \\
    OUTPUT=$bam/$brief_name/$brief_name.dedup.bam                           \\
    METRICS_FILE=$bam/$brief_name/$brief_name.picard_info.txt               \\
    REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT

$python_exe $py_ExtRead -l $ext_len -p $end                                 \\
    $bam/$brief_name/$brief_name.dedup.bam $bed_rep/$brief_name

$samtools_exe view -b -F 1548 -q 30 $bam/$brief_name/$brief_name.dedup.bam |\\
$bedtools_exe bamtobed -i -                                                |\\
    awk 'BEGIN{FS="\\t";OFS="\\t"}{$4="N"; print $0}'                       \\
    >$bed_rep/$brief_name.tagAlign

rm $bam/$brief_name/$brief_name.bam""")
        return l_sh_info

    def s04_bedmrg(self):
        l_sh_info = []
        l_sh_info.append("bed_mrg=$1")
        l_sh_info.append("sh_sortbed=%s/sort_bed.sh" % (self.bin))
        l_sh_info.append("bgzip_exe=%s" % self.sftw_bgzip  )
        l_sh_info.append("tabix_exe=%s" % self.sftw_tabix  )
        l_sh_info.append("samtools_exe=%s"% self.sftw_samtools)
        l_sh_info.append("bedtools_exe=%s"% self.sftw_bedtools)
        l_sh_info.append("""
shift

for i in $@
do
    cat $i.unique.bed
done >$bed_mrg.unique.bed

for i in $@
do
    cat $i.tagAlign
done >$bed_mrg.tagAlign

for i in $bed_mrg $@
do
    sh $sh_sortbed $i.unique.bed $i.unique.sort.bed $bgzip_exe $tabix_exe
done

for i in $bed_mrg $@
do
    sh $sh_sortbed $i.tagAlign $i.sort.tagAlign $bgzip_exe $tabix_exe
done""")

        return l_sh_info

    def s04_2_bedmrg_multi(self):
        l_sh_info = []
        l_sh_info.append("bed_mrg=$1")
        l_sh_info.append("sh_sortbed=%s/sort_bed.sh" % self.bin )
        l_sh_info.append("bgzip_exe=%s" % self.sftw_bgzip  )
        l_sh_info.append("tabix_exe=%s" % self.sftw_tabix  )
        l_sh_info.append("""
shift

for i in $@
do
    cat $i.multi.bed
done >$bed_mrg.multi.bed

for i in $bed_mrg $@
do
    sh $sh_sortbed $i.multi.bed $i.multi.sort.bed $bgzip_exe $tabix_exe
done""")
        return l_sh_info


    def s05_1_spp_rep_shiftSize(self):
        l_sh_info = []
        l_sh_info.append("samp=$1")
        l_sh_info.append("R_exe=%s"       % self.sftw_R     )
        l_sh_info.append("spp_exe=%s"     % self.sftw_spp   )
        l_sh_info.append("bed_rep_dir=%s" % self.dir_bed_rep)
        l_sh_info.append("spp_rep_dir=%s" % self.dir_spp_rep_shiftSize)
        l_sh_info.append("""
rm $spp_rep_dir/$samp/out/out.tab

$R_exe $spp_exe -c=$bed_rep_dir/$samp.sort.tagAlign.gz                      \\
    -odir=$spp_rep_dir/$samp/test -savp -rf                                 \\
    -out=$spp_rep_dir/$samp/out/out.tab""")

        return l_sh_info


    def s05_2_spp_mrg_shiftSize(self):
        l_sh_info = []
        l_sh_info.append("merge=$1")
        l_sh_info.append("R_exe=%s"       % self.sftw_R     )
        l_sh_info.append("spp_exe=%s"     % self.sftw_spp   )
        l_sh_info.append("bed_mrg_dir=%s" % self.dir_bed_mrg)
        l_sh_info.append("spp_mrg_dir=%s" % self.dir_spp_mrg_shiftSize)
        l_sh_info.append("""
rm $spp_mrg_dir/$merge/out/out.tab

$R_exe $spp_exe -c=$bed_mrg_dir/$merge.sort.tagAlign.gz                     \\
    -odir=$spp_mrg_dir/$merge/test -savp -rf                                \\
    -out=$spp_mrg_dir/$merge/out/out.tab""")

        return l_sh_info

    def s06_1_macs2PeakRep(self,ref):
        """ ref, effective genome size. It can be 1.0e+9 or 1000000000,
            or shortcuts:'hs' for human (2.7e9), 'mm' for mouse
            (1.87e9), 'ce' for C. elegans (9e7) and 'dm' for
            fruitfly (1.2e8), Default:hs """
        l_sh_info = []
        l_sh_info.append("brief_name=$1")
        l_sh_info.append("control=$2")
        l_sh_info.append("pvalue=$3")
        l_sh_info.append("shift_size=$4")
        l_sh_info.append("py_exe=%s"       % self.sftw_py     )
        l_sh_info.append("macs2_exe=%s"    % self.sftw_macs2  )
        l_sh_info.append("dir_bed_rep=%s"  % self.dir_bed_rep )
        l_sh_info.append("dir_bed_mrg=%s"  % self.dir_bed_mrg )
        l_sh_info.append("dir_Peak_rep=%s" % self.dir_Peak_rep)
        l_sh_info.append("get_psudoCount=%s" % self.sftw_get_psudoCount)
        l_sh_info.append("ref=%s" % ref)
        l_sh_info.append("""

treat_bed=${dir_bed_rep}/${brief_name}.sort.tagAlign.gz
ctrl_bed=${dir_bed_mrg}/${control}.sort.tagAlign.gz
name_prefix=${dir_Peak_rep}/${brief_name}/${brief_name}

$py_exe $macs2_exe callpeak -t $treat_bed                 -c $ctrl_bed      \\
    -f BED -n ${name_prefix}_VS_Input     -g $ref -p $pvalue --to-large     \\
    --nomodel --extsize $shift_size -B

sh $get_psudoCount $treat_bed

$py_exe $macs2_exe callpeak -t $treat_bed.pr1.tagAlign.gz -c $ctrl_bed      \\
    -f BED -n ${name_prefix}.pr1_VS_Input -g $ref -p $pvalue --to-large     \\
    --nomodel --extsize $shift_size -B

$py_exe $macs2_exe callpeak -t $treat_bed.pr2.tagAlign.gz -c $ctrl_bed      \\
    -f BED -n ${name_prefix}.pr2_VS_Input -g $ref -p $pvalue --to-large     \\
    --nomodel --extsize $shift_size -B""")
        return l_sh_info


    def s06_2_macs2PeakMrg(self,ref):
        """ ref, effective genome size. It can be 1.0e+9 or 1000000000,
            or shortcuts:'hs' for human (2.7e9), 'mm' for mouse
            (1.87e9), 'ce' for C. elegans (9e7) and 'dm' for
            fruitfly (1.2e8), Default:hs """

        l_sh_info = []
        l_sh_info.append("merge_name=$1")
        l_sh_info.append("control=$2")
        l_sh_info.append("pvalue=$3")
        l_sh_info.append("shift_size=$4")
        l_sh_info.append("py_exe=%s"         % self.sftw_py     )
        l_sh_info.append("macs2_exe=%s"      % self.sftw_macs2  )
        l_sh_info.append("dir_bed_mrg=%s"    % self.dir_bed_mrg )
        l_sh_info.append("dir_Peak_mrg=%s"   % self.dir_Peak_mrg)
        l_sh_info.append("get_psudoCount=%s" % self.sftw_get_psudoCount)
        l_sh_info.append("ref=%s"            % ref)

        l_sh_info.append("""
treat_bed=${dir_bed_mrg}/${merge_name}.sort.tagAlign.gz
ctrl_bed=${dir_bed_mrg}/${control}.sort.tagAlign.gz
name_prefix=${dir_Peak_mrg}/${merge_name}/${merge_name}

$py_exe $macs2_exe callpeak -t $treat_bed                 -c $ctrl_bed      \\
    -f BED -n ${name_prefix}_VS_Input     -g $ref -p $pvalue --to-large     \\
    --nomodel --extsize $shift_size -B

sh $get_psudoCount $treat_bed

$py_exe $macs2_exe callpeak -t $treat_bed.pr1.tagAlign.gz -c $ctrl_bed      \\
    -f BED -n ${name_prefix}.pr1_VS_Input -g $ref -p $pvalue --to-large     \\
    --nomodel --extsize $shift_size -B

$py_exe $macs2_exe callpeak -t $treat_bed.pr2.tagAlign.gz -c $ctrl_bed      \\
    -f BED -n ${name_prefix}.pr2_VS_Input -g $ref -p $pvalue --to-large     \\
    --nomodel --extsize $shift_size -B""")
        return l_sh_info


    def s07_1_IDR_prepare(self):
        l_sh_info = []
        l_sh_info.append("mrg_peak=$1")
        l_sh_info.append("top_peak=$2")
        l_sh_info.append("")
        l_sh_info.append("bgzip_exe=%s"    % self.sftw_bgzip  )
        l_sh_info.append("tabix_exe=%s"    % self.sftw_tabix  )
        l_sh_info.append("Peak_rep_dir=%s" % self.dir_Peak_rep)
        l_sh_info.append("Peak_mrg_dir=%s" % self.dir_Peak_mrg)
        l_sh_info.append("Peak_idr_dir=%s" % self.dir_Peak_idr)
        l_sh_info.append("""
prefix=$Peak_mrg_dir/${mrg_peak}/${mrg_peak}
shift
shift

for i in $@
do
   sort -k 8nr,8nr $Peak_rep_dir/${i}/${i}_VS_Input_peaks.narrowPeak      | \\
       head -n     $top_peak   |  $bgzip_exe -cf                            \\
       >$Peak_rep_dir/${i}/${i}_VS_Input_peaks.regionPeak.gz             && \\

   sort -k 8nr,8nr $Peak_rep_dir/${i}/${i}.pr1_VS_Input_peaks.narrowPeak  | \\
        head -n     $top_peak  |  $bgzip_exe -cf                            \\
       >$Peak_rep_dir/${i}/${i}.pr1_VS_Input_peaks.regionPeak.gz         && \\

   sort -k 8nr,8nr $Peak_rep_dir/${i}/${i}.pr2_VS_Input_peaks.narrowPeak  | \\
        head -n     $top_peak  |  $bgzip_exe -cf                            \\
       >$Peak_rep_dir/${i}/${i}.pr2_VS_Input_peaks.regionPeak.gz         && \\
   mv           $Peak_rep_dir/${i}/${i}_VS_Input_peaks.regionPeak.gz*       \\
                $Peak_idr_dir/${mrg_peak}                                && \\
   mv           $Peak_rep_dir/${i}/${i}.pr1_VS_Input_peaks.regionPeak.gz*   \\
                $Peak_idr_dir/${mrg_peak}                                && \\
   mv           $Peak_rep_dir/${i}/${i}.pr2_VS_Input_peaks.regionPeak.gz*   \\
                $Peak_idr_dir/${mrg_peak}
done

sort -k 8nr,8nr  ${prefix}_VS_Input_peaks.narrowPeak                      | \\
    head -n $top_peak  | $bgzip_exe -cf                                     \\
    >${prefix}_VS_Input_peaks.regionPeak.gz                              && \\
sort -k 8nr,8nr  ${prefix}.pr1_VS_Input_peaks.narrowPeak                  | \\
    head -n $top_peak  | $bgzip_exe -cf                                     \\
    >${prefix}.pr1_VS_Input_peaks.regionPeak.gz                          && \\
sort -k 8nr,8nr  ${prefix}.pr2_VS_Input_peaks.narrowPeak                  | \\
    head -n $top_peak  | $bgzip_exe -cf                                     \\
    >${prefix}.pr2_VS_Input_peaks.regionPeak.gz                          && \\
mv ${prefix}_VS_Input_peaks.regionPeak.gz                                   \\
    $Peak_idr_dir/${mrg_peak}                                            && \\
mv ${prefix}.pr1_VS_Input_peaks.regionPeak.gz                               \\
    $Peak_idr_dir/${mrg_peak}                                            && \\
mv ${prefix}.pr2_VS_Input_peaks.regionPeak.gz                               \\
    $Peak_idr_dir/${mrg_peak}""")
        return l_sh_info

    def s07_2_usingMacs2Peak(self):
        l_sh_info = []
        l_sh_info.append("in_sam1=$1")
        l_sh_info.append("in_sam2=$2")
        l_sh_info.append("merge_name=$3")
        l_sh_info.append("subtype=$4")
        l_sh_info.append("R_exe=%s"           % self.sftw_R       )
        l_sh_info.append("batchIDR_exe=%s"    % self.sftw_batchIDR)
        l_sh_info.append("genome_table=%s.len"% self.genome_ref       )
        l_sh_info.append("idr_dir=%s"         % self.dir_Peak_idr )
        l_sh_info.append("""
if [ ! -f $genome_table ]
then
    awk '{print $1 "\\t" $2}' %s.fai >$genome_table
fi

region1=$idr_dir/$merge_name/${in_sam1}_VS_Input_peaks.regionPeak.gz
region2=$idr_dir/$merge_name/${in_sam2}_VS_Input_peaks.regionPeak.gz
out=$idr_dir/$merge_name/$subtype/${in_sam1}_VS_${in_sam2}

$R_exe $batchIDR_exe $region1 $region2                                      \\
  -1  $out 0 F p.value $genome_table""")
        return l_sh_info

    def s07_3_IDR_passPeaks(self):
        l_sh_info = []

        l_sh_info.append("merge_name=$1")
        l_sh_info.append("peak=$2")
        l_sh_info.append("dir_Peak_idr=%s" % (self.dir_Peak_idr) )
        l_sh_info.append("bedtools_exe=%s" % self.sftw_bedtools)
        l_sh_info.append("bgzip_exe=%s"    % self.sftw_bgzip   )
        l_sh_info.append("tabix_exe=%s"    % self.sftw_tabix   )

        l_sh_info.append("""
prefix=${dir_Peak_idr}/${merge_name}/${merge_name}_VS_Input_peaks
infile=${prefix}.regionPeak.gz
outfile=${prefix}.conservative.regionPeak.gz

zcat $infile | sort -k7nr,7nr | head -n $peak |                             \\
    $bedtools_exe sort -i /dev/stdin |                                      \\
    $bgzip_exe -cf >$outfile                                             && \\
    $tabix_exe -f -p bed -s 1 -b 2 -e 3 -S 0 $outfile""")
        return l_sh_info

    def s08_macs2BroadPeakRep(self,ref="hs"):
        """ ref, effective genome size. It can be 1.0e+9 or 1000000000,
            or shortcuts:'hs' for human (2.7e9), 'mm' for mouse
            (1.87e9), 'ce' for C. elegans (9e7) and 'dm' for
            fruitfly (1.2e8), Default:hs """
        l_sh_info = []
        l_sh_info.append("brief_name=$1")
        l_sh_info.append("control=$2")
        l_sh_info.append("pvalue=$3")
        l_sh_info.append("py_exe=%s"      % self.sftw_py     )
        l_sh_info.append("macs2_exe=%s"   % self.sftw_macs2  )
        l_sh_info.append("dir_bed_rep=%s" % self.dir_bed_rep )
        l_sh_info.append("dir_bed_mrg=%s" % self.dir_bed_mrg )
        l_sh_info.append("get_psudoCount=%s" % self.sftw_get_psudoCount)
        l_sh_info.append("dir_broad_rep=%s"  % self.dir_BroadPeak_rep  )
        l_sh_info.append("ref=%s"         % ref)
        l_sh_info.append("""
treat_bed=${dir_bed_rep}/${brief_name}.sort.tagAlign.gz
ctrl_bed=${dir_bed_mrg}/${control}.sort.tagAlign.gz
out_dir=${dir_broad_rep}/${brief_name}

$py_exe $macs2_exe callpeak  --broad  -t $treat_bed       -c $ctrl_bed      \\
  -f BED -n ${brief_name}_VS_Input --outdir $out_dir -g $ref -p $pvalue -B""")

        return l_sh_info


    def s09_macs2BroadPeakMrg(self,ref="hs"):
        """ ref, effective genome size. It can be 1.0e+9 or 1000000000,
            or shortcuts:'hs' for human (2.7e9), 'mm' for mouse
            (1.87e9), 'ce' for C. elegans (9e7) and 'dm' for
            fruitfly (1.2e8), Default:hs """
        l_sh_info = []
        l_sh_info.append("merge_name=$1")
        l_sh_info.append("control=$2")
        l_sh_info.append("pvalue=$3")
        l_sh_info.append("py_exe=%s"      % self.sftw_py     )
        l_sh_info.append("macs2_exe=%s"   % self.sftw_macs2  )
        l_sh_info.append("dir_bed_mrg=%s" % self.dir_bed_mrg )
        l_sh_info.append("get_psudoCount=%s"% self.sftw_get_psudoCount)
        l_sh_info.append("dir_broad_mrg=%s"  % self.dir_BroadPeak_mrg  )
        l_sh_info.append("ref=%s"         % ref)
        l_sh_info.append("""
treat_bed=${dir_bed_mrg}/${merge_name}.sort.tagAlign.gz
ctrl_bed=${dir_bed_mrg}/${control}.sort.tagAlign.gz
out_dir=${dir_broad_mrg}/${merge_name}

$py_exe $macs2_exe callpeak  --broad  -t $treat_bed       -c $ctrl_bed      \\
  -f BED -n ${merge_name}_VS_Input --outdir $out_dir -g $ref -p $pvalue -B""")

        return l_sh_info


    def s10_sortbdg(self):
        l_sh_info = []
        l_sh_info.append("sam=$1")
        l_sh_info.append("dir=$2")
        l_sh_info.append("sort_bdg_exe=%s" % self.sftw_sort_bdg)
        l_sh_info.append("""
sh $sort_bdg_exe $dir/$sam/${sam}_VS_Input_treat_pileup.bdg                 \\
    $dir/$sam/${sam}_VS_Input_treat_pileup.sort.bdg                      && \\
sh $sort_bdg_exe $dir/$sam/${sam}_VS_Input_control_lambda.bdg               \\
    $dir/$sam/${sam}_VS_Input_control_lambda.sort.bdg""")
        return l_sh_info

    def s11_makeIGV(self):
        l_sh_info = []
        l_sh_info.append("sam_in=$1")
        l_sh_info.append("sam_out=$2")
        l_sh_info.append("dir_in=$3")
        l_sh_info.append("dir_out=%s"      % self.dir_Peak_TDF  )
        l_sh_info.append("igvtools_exe=%s" % self.sftw_igvtools )
        l_sh_info.append("genome_fa=%s"    % self.genome_ref        )
        l_sh_info.append("""
ln -s $dir_in/$sam_in/${sam_in}_VS_Input_treat_pileup.bdg                   \\
    $dir_in/$sam_in/${sam_in}_VS_Input_treat_pileup.bedGraph
ln -s $dir_in/$sam_in/${sam_in}_VS_Input_control_lambda.bdg                 \\
    $dir_in/$sam_in/${sam_in}_VS_Input_control_lambda.bedGraph

if [ ! -d $dir_out/${sam_out} ]
    then mkdir $dir_out/${sam_out}
fi

$igvtools_exe toTDF                                                         \\
    $dir_in/$sam_in/${sam_in}_VS_Input_treat_pileup.bedGraph                \\
    $dir_out/${sam_out}/${sam_in}_treat_pileup.tdf   $genome_fa &&          \\

$igvtools_exe toTDF                                                         \\
    $dir_in/$sam_in/${sam_in}_VS_Input_control_lambda.bedGraph              \\
    $dir_out/${sam_out}/${sam_in}_control_lambda.tdf $genome_fa

rm  $dir_in/$sam_in/${sam_in}_VS_Input_treat_pileup.bedGraph                \\
    $dir_in/$sam_in/${sam_in}_VS_Input_control_lambda.bedGraph""")

        return l_sh_info

    def s11_makeIGV_broad(self):
        l_sh_info = []
        l_sh_info.append("sam_in=$1")
        l_sh_info.append("sam_out=$2")
        l_sh_info.append("dir_in=$3")
        l_sh_info.append("dir_out=%s"      % self.dir_BroadPeak_TDF  )
        l_sh_info.append("igvtools_exe=%s" % self.sftw_igvtools )
        l_sh_info.append("genome_fa=%s"    % self.genome_ref        )
        l_sh_info.append("""
ln -s $dir_in/$sam_in/${sam_in}_treat_pileup.bdg                            \\
      $dir_in/$sam_in/${sam_in}_treat_pileup.bedGraph
ln -s $dir_in/$sam_in/${sam_in}_control_lambda.bdg                          \\
      $dir_in/$sam_in/${sam_in}_control_lambda.bedGraph

if [ ! -d $dir_out/${sam_out} ]
    then mkdir $dir_out/${sam_out}
fi

$igvtools_exe toTDF $dir_in/$sam_in/${sam_in}_treat_pileup.bedGraph         \\
    $dir_out/${sam_out}/${sam_in}_treat_pileup.broadPeaks.tdf  $genome_fa &&\\

$igvtools_exe toTDF $dir_in/$sam_in/${sam_in}_control_lambda.bedGraph       \\
    $dir_out/${sam_out}/${sam_in}_control_lambda.broadPeaks.tdf $genome_fa

ln -s $dir_in/$sam_in/${sam_in}_peaks.broadPeak                             \\
    $dir_in/$sam_in/${sam_in}_peaks.broadPeak.bed

$igvtools_exe count $dir_in/$sam_in/${sam_in}_peaks.broadPeak.bed           \\
    $dir_out/${sam_out}/${sam_in}_peaks.broadPeak.tdf $genome_fa         && \\

rm  $dir_in/$sam_in/${sam_in}_treat_pileup.bedGraph                         \\
    $dir_in/$sam_in/${sam_in}_control_lambda.bedGraph                       \\
    $dir_in/$sam_in/${sam_in}_peaks.broadPeak.bed""")
        return l_sh_info


    def s12_PeakGeneRegion(self,
            TSS_genebody_up,TSS_genebody_down,TSS_promoter_up,
            TSS_promoter_down,ext_binlen=50,body_bincnt=100,
            tss_binlen=1):

        prefix  = ".".join(self.refGeneTxt.split(".")[:-1])
        genebody_region = "%s.up%d_down%d.%sBsorted.longestTid.bed"           %\
                     (prefix, TSS_genebody_up, TSS_genebody_down, "genebody.")

        getbin_genebody_cpp = "%s/get_bin_ignore_peak/get_bin" % (self.bin)

        l_sh_info = []
        l_sh_info.append("merge_name=$1")
        l_sh_info.append("dir_Peak_mrg=%s"       % self.dir_Peak_mrg      )
        l_sh_info.append("dir_Peak_mrg_TSS=%s"   % self.dir_Peak_mrg_TSS  )
        l_sh_info.append("dir_Peak_mrg_Gene=%s"  % self.dir_Peak_mrg_Gene )
        l_sh_info.append("bedtools_exe=%s"       % self.sftw_bedtools     )
        l_sh_info.append("genebody_region=%s"    % genebody_region        )
        l_sh_info.append("getbin_genebody_cpp=%s"% getbin_genebody_cpp    )
        l_sh_info.append("TSS_genebody_up=%d"    % TSS_genebody_up        )
        l_sh_info.append("TSS_genebody_down=%d"  % TSS_genebody_down      )
        l_sh_info.append("TSS_promoter_up=%d"    % TSS_promoter_up        )
        l_sh_info.append("TSS_promoter_down=%d"  % TSS_promoter_down      )
        l_sh_info.append("ext_binlen=%d"         % ext_binlen             )
        l_sh_info.append("body_bincnt=%d"        % body_bincnt            )
        l_sh_info.append("tss_binlen=%d"         % tss_binlen             )

        l_sh_info.append("""
prefix_Peak=${dir_Peak_mrg}/${merge_name}/${merge_name}
prefix_P=${dir_Peak_mrg_TSS}/${merge_name}/${merge_name}
prefix_G=${dir_Peak_mrg_Gene}/${merge_name}/${merge_name}

count_bdg=${prefix_Peak}_VS_Input_treat_pileup.sort.bdg.gz
stat1=$prefix_P.genebody.up${TSS_genebody_up}_down${TSS_genebody_down}.xls
stat2=$prefix_G.genebody.up${TSS_promoter_up}_down${TSS_promoter_down}.xls

zcat $count_bdg | tail -n +2                                               |\\
    $bedtools_exe intersect -sorted -wo -a /dev/stdin -b $genebody_region  |\\
    $getbin_genebody_cpp -U $TSS_genebody_up -D $TSS_genebody_down          \\
    -T $TSS_promoter_up -b $ext_binlen -B $body_bincnt /dev/stdin           \\
    $stat1 $stat2""")

        return l_sh_info



    def s13_RPM_density_rep(self):
        getbin_readcnt_cpp = "%s/chip01.bed_reads_unique/bed_read"% (self.bin)
        genome_fai = "%s.fai" % (self.genome_ref)

        l_sh_info = []
        l_sh_info.append("brief_name=$1")
        l_sh_info.append("mapped_reads=$2")
        l_sh_info.append("getbin_readcnt_cpp=%s" % getbin_readcnt_cpp    )
        l_sh_info.append("genome_fai=%s"         % genome_fai            )
        l_sh_info.append("dir_bed_rep=%s"        % self.dir_bed_rep      )
        l_sh_info.append("dir_RPM_bins_rep=%s"   % self.dir_RPM_bins_rep )
        l_sh_info.append("""
unique_bed=${dir_bed_rep}/${brief_name}.sort.tagAlign.gz
out_cnt=${dir_RPM_bins_rep}/${brief_name}/${brief_name}.cnt.uniq
out_RPKM=${dir_RPM_bins_rep}/${brief_name}/${brief_name}.RPKM.uniq

zcat $unique_bed                                                           |\\
    $getbin_readcnt_cpp -b 100  $genome_fai /dev/stdin $out_cnt.100

awk '{ print $1*1000000*1000/(100 *'$mapped_reads') }' $out_cnt.100         \\
    >$out_RPKM.100

zcat $unique_bed                                                           |\\
    $getbin_readcnt_cpp -b 1000 $genome_fai /dev/stdin $out_cnt.1kb

awk '{ print $1*1000000*1000/(1000*'$mapped_reads') }' $out_cnt.1kb         \\
    >$out_RPKM.1kb""")

        return l_sh_info



    def s14_RPM_density_mrg(self):
        getbin_readcnt_cpp = "%s/chip01.bed_reads_unique/bed_read"% (self.bin)
        genome_fai = "%s.fai" % (self.genome_ref)

        l_sh_info = []
        l_sh_info.append("merge_name=$1")
        l_sh_info.append("mapped_reads=$2")
        l_sh_info.append("getbin_readcnt_cpp=%s" % getbin_readcnt_cpp    )
        l_sh_info.append("genome_fai=%s"         % genome_fai            )
        l_sh_info.append("dir_bed_mrg=%s"        % self.dir_bed_mrg      )
        l_sh_info.append("dir_RPM_bins_mrg=%s"   % self.dir_RPM_bins_mrg )
        l_sh_info.append("""
unique_bed=${dir_bed_mrg}/${merge_name}.sort.tagAlign.gz
out_cnt=${dir_RPM_bins_mrg}/${merge_name}/${merge_name}.cnt.uniq
out_RPKM=${dir_RPM_bins_mrg}/${merge_name}/${merge_name}.RPKM.uniq

zcat $unique_bed                                                           |\\
    $getbin_readcnt_cpp -b 100  $genome_fai /dev/stdin $out_cnt.100

awk '{ print $1*1000000*1000/(100 *'$mapped_reads') }' $out_cnt.100         \\
    >$out_RPKM.100

zcat $unique_bed                                                           |\\
    $getbin_readcnt_cpp -b 1000 $genome_fai /dev/stdin $out_cnt.1kb

awk '{ print $1*1000000*1000/(1000*'$mapped_reads') }' $out_cnt.1kb         \\
    >$out_RPKM.1kb""")

        return l_sh_info


    def s15_merge_RPKM(self):
        l_sh_info = []
        l_sh_info.append("header=$1")
        l_sh_info.append("window=$2")
        l_sh_info.append("type=$3")
        l_sh_info.append("dir_RPM_mrg=%s" % self.dir_RPM_mrg )
        l_sh_info.append("bgzip_exe=%s"   % self.sftw_bgzip  )
        l_sh_info.append("tabix_exe=%s"   % self.sftw_tabix  )
        l_sh_info.append("column=%s/%s/columns.${window}.bed" % \
                                            (self.Database, self.ref))
        l_sh_info.append("""
shift
shift
shift

prefix=${dir_RPM_mrg}/Merge_RPKM_bin${window}.uniq.${type}

echo -e "$header"  >${prefix}_${window}.tmp                              && \\
paste $column $@  >>${prefix}_${window}.tmp                              && \\
$bgzip_exe -c -f ${prefix}_${window}.tmp >${prefix}_${window}.bed.gz     && \\
$tabix_exe -f -p bed -s 1 -b 2 -e 3 -S 1  ${prefix}_${window}.bed.gz
rm ${prefix}_${window}.tmp""")

        return l_sh_info


    def s16_densityBaselv(self):
        genome_fai = "%s.fai" % (self.genome_ref)

        l_sh_info = []
        l_sh_info.append("sam=$1")
        l_sh_info.append("ref=$2")
        l_sh_info.append("dir=$3")
        l_sh_info.append("bin=%s" % (self.bin))
        l_sh_info.append("ucsc_dir=%s" % (self.sftw_ucsc_dir))
        l_sh_info.append("python_exe=%s" % self.sftw_py)
        l_sh_info.append("genome_fai=%s" % genome_fai)
        l_sh_info.append("bedtools_exe=%s" % (self.sftw_bedtools))
        l_sh_info.append("column=%s/%s/columns.1kb.bed" % (self.Database, \
            self.ref))
        l_sh_info.append("""
$python_exe $bin/norm_density.py                                            \\
    -r $ref $dir/$sam/${sam}_VS_Input_treat_pileup.sort.bdg.gz           && \\
$python_exe $bin/norm_density.py                                            \\
    -r $ref $dir/$sam/${sam}_VS_Input_control_lambda.sort.bdg.gz         && \\

$bedtools_exe intersect -wo -sorted                                         \\
    -a $dir/$sam/${sam}_VS_Input_treat_pileup.sort.norm.bedGraph            \\
    -b $dir/$sam/${sam}_VS_Input_control_lambda.sort.norm.bedGraph         |\\
    $python_exe $bin/minus_control_bdg.py /dev/stdin                        \\
    >$dir/$sam/${sam}_treat_minus_control.sort.norm.bedGraph

$ucsc_dir/bedGraphToBigWig                                                  \\
    $dir/$sam/${sam}_treat_minus_control.sort.norm.bedGraph                 \\
    $genome_fai $dir/$sam/${sam}_treat_minus_control.sort.norm.bw        && \\

rm $dir/$sam/${sam}_VS_Input_control_lambda.sort.norm.bedGraph              \\
   $dir/$sam/${sam}_VS_Input_treat_pileup.sort.norm.bedGraph

awk '{print $0 "\\t" $1 ":" $2 "-" $3 }' $column                           |\\
$ucsc_dir/bigWigAverageOverBed                                              \\
    $dir/$sam/${sam}_treat_minus_control.sort.norm.bw /dev/stdin            \\
    /dev/stdout | awk '{print $6}' >$dir/$sam/${sam}.1kb.norm_avg.xls

""")

        return l_sh_info

    def s17_merge_RPKM(self):
        l_sh_info = []
        l_sh_info.append("header=$1")
        l_sh_info.append("window=$2")
        l_sh_info.append("type=$3")
        l_sh_info.append("dir_RPM_mrg=%s" % self.dir_RPM_mrg )
        l_sh_info.append("bgzip_exe=%s"   % self.sftw_bgzip  )
        l_sh_info.append("tabix_exe=%s"   % self.sftw_tabix  )
        l_sh_info.append("column=%s/%s/columns.${window}.bed" % \
                                            (self.Database, self.ref))
        l_sh_info.append("""
shift
shift
shift

prefix=${dir_RPM_mrg}/Merge_density.${type}

echo -e "$header"  >${prefix}_${window}.tmp                              && \\
paste $column $@  >>${prefix}_${window}.tmp                              && \\
$bgzip_exe -c -f ${prefix}_${window}.tmp >${prefix}.${window}.bed.gz     && \\
$tabix_exe -f -p bed -s 1 -b 2 -e 3 -S 1  ${prefix}.${window}.bed.gz
rm ${prefix}_${window}.tmp""")

        return l_sh_info