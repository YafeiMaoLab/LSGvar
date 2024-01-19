source config.json
nowdic=$(pwd)
minimap2 -t 12 -cx asm20 --secondary=no --eqx -Y -K 8G -s 1000 $ref_path $hap1_path  -o hm_prihap1.paf
python "${tool_path}scripts/one2multi_filter.py" -m ${mappingtsv} -f  hm_prihap1.paf -1 6 -2 1 |awk '{print $6,$7,$1,$2}' >p_c_chrlen1.txt

python "${tool_path}scripts/one2multi_filter.py" -m ${mappingtsv} -f  hm_prihap1.paf -1 6 -2 1 \
	| rustybam trim-paf \
	| rustybam break-paf --max-size 5000 \
	| rustybam filter --paired-len 100000 \
	| awk '$10 >= 20000' \
	> hm_prihap1.flt.paf
if [ ! -d "${nowdic}/saffireh1" ]; then
	mkdir "${nowdic}/saffireh1"
fi
if [ "$3" = "ctn" ]; then
    Rscript "${tool_path}scripts/chaos_filt.r" "/home/jmhan/SDR/run_chimpanzee/p_c_chrlen1.txt" "${nowdic}/saffireh1/" "/home/jmhan/SDR/run_chimpanzee/hm_prihap1.flt.paf" "${nowdic}/afterchaos_hap1.flt.paf" $1 $2 $3
elif [ "$3" = "cts" ]; then
    Rscript "${tool_path}scripts/chaos_filt.r" "/home/jmhan/SDR/run_chimpanzee/p_c_chrlen1.txt" "${nowdic}/saffireh1/" "/home/jmhan/SDR/run_chimpanzee/hm_prihap1.flt.paf" "${nowdic}/afterchaos_hap1.flt.paf" $1 $2 $3 ${centro} ${telome}
else
    echo "please input ctn/cts parameter"
fi

awk 'BEGIN{OFS="\t"}{print $6, $8, $9, ($8+$9)/2, $1, $3, $4, ($3+$4)/2, $5}' "${nowdic}/afterchaos_hap1.flt.paf"  |sort -k1,1 -k2,2n > "${nowdic}/h1syntenic_blocks.tsv"
if [ ! -d "${nowdic}/result" ]; then
	mkdir "${nowdic}/result"
fi
Rscript "${tool_path}scripts/SDR.r" "${tool_path}scripts/SDRfun.r" "${nowdic}/h1syntenic_blocks.tsv" "${nowdic}/p_c_chrlen1.txt" "${nowdic}/result/" 
cd result/
cat ./*  |sort -r |uniq >end.txt