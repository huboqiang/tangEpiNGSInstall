### Check for required parameters

if [ -z ${ref} ]
then
    echo "ERROR: Reference must be specified."
    exit 1
fi


if [ ${type} == "ChIP" ]
then
    for i in `tail -n +2 /fastq/sample.tab.xls | awk '{print $1}'`
    do
      mkdir -p /home/analyzer/project/ChIP_test/00.0.raw_data/$i &&\
      ln -s /fastq/$i/* /home/analyzer/project/ChIP_test/00.0.raw_data/$i
    done
    cd /home/analyzer/project/ChIP_test
    
    cp /settings/run_chipseq.py     /home/analyzer/module/ChIP/
    cp /settings/scripts_chipseq.py /home/analyzer/module/ChIP/settings/scripts.py
    
    cp /fastq/sample.tab.xls ./
    python /home/analyzer/module/ChIP/run_chipseq.py --ref ${ref} sample.tab.xls
    
    mkdir -p result result/peaks result/bigwig result/tables
    cp 03.2.Peak_mrg/*/*_treat_minus_control.sort.norm.bw result/bigwig
    cp 03.3.Peak_idr/*/*.conservative.regionPeak.gz*      result/peaks
    cp StatInfo/* result/tables
fi

if [ ${type} == "RNA" ]
then
    for i in `tail -n +2 /fastq/sample.tab.xls | awk '{print $1}'`
    do
      mkdir -p /home/analyzer/project/RNA_test/00.0.raw_data/$i &&\
      ln -s /fastq/$i/* /home/analyzer/project/RNA_test/00.0.raw_data/$i
    done
    cd /home/analyzer/project/RNA_test
    ls /home/analyzer/bin
    
    cp /settings/run_mRNA.py     /home/analyzer/module/RNA_v2/
    cp /settings/scripts_mRNA.py /home/analyzer/module/RNA_v2/settings/scripts.py
    
    cp /fastq/sample.tab.xls ./
    python /home/analyzer/module/RNA_v2/run_mRNA.py --ref ${ref} --given_GTF /home/analyzer/database_RNA/mm10/refGene.gtf sample.tab.xls
    mkdir -p result result/count result/fpkm result/repeat result/table
    cp 02.HTSeq_result/*                             result/count
    cp 05.2.cufflinksMerge_known/merge.FPKM.gene.xls result/fpkm
    cp 07.merge_repFPKM/merge.Repeat.*               result/repeat
    cp StatInfo/*                                    result/table
fi