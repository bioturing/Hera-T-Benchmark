mkdir $1
cd $1
echo "Download fastqs"
wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/$1/$1_fastqs.tar
tar -xvf $1_fastqs.tar
rm $1_fastqs.tar
mv ./* raw_read

echo "Download cellranger result"
wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/$1/$1_filtered_feature_bc_matrix.tar.gz
tar -xzvf $1_filtered_feature_bc_matrix.tar.gz
mv filtered_feature_bc_matrix/$2 cellranger
rm -r filtered_feature_bc_matrix
rm $1_filtered_feature_bc_matrix.tar.gz
cd ../
