# Run HeraT
echo "\time -v ../../bin/HeraT count -x ../../reference/HeraT_v2/"$2"/"$2" -t 56 -o "$1"/Hera-T -p "$1" -1 "$(ls $1/raw_read/*_R1_*)" -2 "$(ls $1/raw_read/*_R2_*)
\time -v ../../bin/HeraT count -x ../../reference/HeraT_v2/$2/$2 -t 56 -o $1/Hera-T -p $1 -1 $(ls $1/raw_read/*_R1_*) -2 $(ls $1/raw_read/*_R2_*) &> $1/Hera-T.performance

# Run alevin
mkdir $1/alevin
echo "../../bin/salmon alevin -l ISR --chromium -i ../../reference/alevin_v2/"$2"/ -p 56 -o "$1"/alevin --tgMap ../../reference/alevin_v2/"$2"/trans2gene.tsv --dumpCsvCounts -1 "$(ls $1/raw_read/*_R1_*)" -2 "$(ls $1/raw_read/*_R2_*)
\time -v ../../bin/salmon/bin/salmon alevin -l ISR --chromium -i ../../reference/alevin_v2/$2/ -p 56 -o $1/alevin --tgMap ../../reference/alevin_v2/$2/trans2gene.tsv --dumpCsvCounts -1 $(ls $1/raw_read/*_R1_*) -2 $(ls $1/raw_read/*_R2_*) &> $1/alevin.perfomance
mkdir $1/alevin/result
awk '{print $1" gene"}' $1/alevin/alevin/quants_mat_cols.txt > $1/alevin/result/genes.tsv
cp $1/alevin/alevin/quants_mat_rows.txt $1/alevin/result/barcodes.tsv
python ../../scripts/parse_alevin_mat.py $1/alevin/alevin/quants_mat.csv $1/alevin/result/matrix.mtx

# Run Cellranger
cd $1/
\time -v ../../../bin/cellranger-2.1.1/cellranger count --id cellranger_rerun --fastqs=raw_read/ --transcriptome=../../../reference/cellranger_v2/$2/ --nosecondary &> cellranger.performance
