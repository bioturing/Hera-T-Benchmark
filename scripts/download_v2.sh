mkdir $1
cd $1
echo "Download fastqs"
wget http://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/2.1.0/$1/$1_fastqs.tar
tar -xvf $1_fastqs.tar
rm $1_fastqs.tar
mv ./* raw_read

echo "Download cellranger result"
wget http://cf.10xgenomics.com/samples/cell-exp/2.1.0/$1/$1_filtered_gene_bc_matrices.tar.gz
tar -xzvf $1_filtered_gene_bc_matrices.tar.gz
mv filtered_gene_bc_matrices/$2 cellranger
rm -r filtered_gene_bc_matrices
rm $1_filtered_gene_bc_matrices.tar.gz
cd ../
