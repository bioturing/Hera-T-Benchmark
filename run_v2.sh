cd data/v2
awk '{system("sh ../../scripts/download_v2.sh "$1" "$2)}' list
awk '{system("sh ../../scripts/count_v2.sh "$1" "$2)}' list
cd ../../
mkdir result/v2
awk '{system("scripts/benchmark data/v2/"$1"/cellranger/ data/v2/"$1"/Hera-T/")}' data/v2/list &> result/v2/cellranger_vs_Hera-T
awk '{system("scripts/benchmark data/v2/"$1"/cellranger/ data/v2/"$1"/alevin/")}' data/v2/list &> result/v2/cellranger_vs_alevin
mkdir result/distribution
awk '{system("scripts/benchmark_dis data/v2/"$1"/cellranger/ data/v2/"$1"/Hera-T/ > result/v2/distribution/"$1)}' data/v2/list 
