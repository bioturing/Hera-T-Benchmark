cd data/v2
awk '{system("sh ../../scripts/download_v3.sh "$1" "$2)}' list
awk '{system("sh ../../scripts/count_v3.sh "$1" "$2)}' list
cd ../../
mkdir result/v2
awk '{system("scripts/benchmark data/v3/"$1"/cellranger/ data/v3/"$1"/Hera-T/")}' data/v3/list &> result/v3/cellranger_vs_Hera-T
awk '{system("scripts/benchmark data/v3/"$1"/cellranger/ data/v3/"$1"/alevin/")}' data/v3/list &> result/v3/cellranger_vs_alevin
mkdir result/distribution
awk '{system("scripts/benchmark_dis data/v3/"$1"/cellranger/ data/v3/"$1"/Hera-T/ > result/v3/distribution/"$1)}' data/v3/list 
