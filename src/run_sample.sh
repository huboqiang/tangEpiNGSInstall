for i in `tail -n +2 /fastq/sample.tab.xls | awk '{print $1}'`
do
  mkdir -p /home/analyzer/project/ChIP_test/00.0.raw_data/$i &&\
  ln -s /fastq/$i/* /home/analyzer/project/ChIP_test/00.0.raw_data/$i
done
cd /home/analyzer/project/ChIP_test

cp /fastq/sample.tab.xls ./
python /home/analyzer/module/ChIP/run_chipseq.py --ref mm10 sample.tab.xls

mkdir -p result result/peaks result/bigwig result/tables
cp 03.2.Peak_mrg/*/*_treat_minus_control.sort.norm.bw result/bigwig
cp 03.3.Peak_idr/*/*.conservative.regionPeak.gz*      result/peaks
cp StatInfo/* result/tables
