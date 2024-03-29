# LSGvar



### Dependencies




- **minimap2**
- **rustybam**

- **dbscan**

- **dplyr**

- **data.table**

- **IRanges**

## WorkFlow




### data:

**query**: Four primates(10 haplotypes): Orangutans(orangutan1:Sumatran orangutan,orangutan2:Bornean orangutan)、Chimpanzee、gorilla、bonobo

*datasourse:[marbl/Primates: Complete assemblies of non-human primate genomes (github.com)](https://github.com/marbl/Primates?tab=readme-ov-file)*

**reference**:T2T-CHM13v2.0

*datasourse: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz*

(note：Due to extensive reversals in some of the primitive sequences of primates, processing has been carried out.)



### Run scripts

Download the SDR script through the following steps:

```shell
git clone https://github.com/xiaomiyongyuan/SDR.git
```

Please create a new empty folder to store the run results and navigate into that folder. Copy the `SDR/run.sh` script into this folder.

Configure `config.json`: To start the run, three essential input files are required :`reference genome`, `aligned genome`, and the `artificial chromosome pairing file` (refer to examples for guidance). Additionally, the `centromere and telomere file` is optional (alignments in this region will be filtered out during structural variation computation)."



```shell
tool_path="/home/SDR/"
ref_path="/home/chm13v2.0.fa"
hap1_path="/home/query_h1.fa"
mappingtsv="/home/SDR/examples/chromosome_mapping.tsv"
centro="/home/SDR/examples/T2Tdatabase/hm_centroend.tsv"
telome="/home/SDR/examples/T2Tdatabase/hm_teloend.tsv"
```

Run the shell script(The current initial version of the code has not been updated to the Snakemake workflow yet.)

```shell
bash run.sh 200000 300000 cts
```

Explanation for the following three parameters: The first two parameters are filtering criteria. The first one is the clustering parameter for filtering, where a smaller value results in a stricter filter(Recommended parameters：200000), capable of removing more segments. The second parameter(Recommended parameters：300000) is the desired deletion length for alignments, where a larger value enforces a stricter filter,the last one you can choose cts or ctn.

'`cts`' indicates inputting telomere and centromere fragments from the reference genome for filtering.
'`ctn`' indicates no input of telomere and centromere files. (If the reference genome is T2T-CHM13, refer to the examples/T2Tdatabase for telomere and centromere files).

### SV-annotation：

The final result can be found in `/result/end.txt`.

|                   | SV-annotation                                                |
| ----------------- | ------------------------------------------------------------ |
| SV_(length<10k)   | DEL(deletion)、INS(insertion)、DUP(duplication)、TRANS(translocation)、INV(inversion)、NM(no-matched) |
| SDR_(length>=10k) | DEL(deletion)、INS(insertion)、DUP(duplication)、TRANS(translocation)、INV(inversion)、NM(no-matched) |
| COMPLEX           | complex regions                                              |



( The initial version may still have a few small issues, for reference.)
