# Shweta_et_al_2023
Regulation of replicative histone RNA metabolism by the histone chaperone ASF1  

<br/><br/>

Scripts must be executed from the root of the repository folder.

## RNA-seq processing and read counts

`Scripts/process.RNA.sh` perform mapping of the reads on the genome reference and the counting of reads/fragments per gene annotation. Usually ran on a computer cluster.  
Usage: bash process.RNA.sh <run>
> <run>    accession number of the FASTQ file(s) in a $PWD/data directory
>
>          i.e.  $PWD/data/<run>.fastq.gz (single-end)
>                $PWD/data/<run>.R[12].fastq.gz (paired-end)

Requirements:
- hisat2
- samtools
- featureCounts
- GRCh38.104 homo sapiens transcriptome index, build with hisat2 using genome sequence and transcriptome annotations from EnsEMBL.


## Differential Gene Expression analysis

> ### Evaluate differentially expressed genes (DEGs)
> 
> `Scripts/dge_analysis_0.2.py` perform differential expression analysis and find the DEGs.


> ### Volcano plot of genes in asynchronous cells
> 
> `Scripts/Volcano_plot.r` build the volcano plot of Fig2a.


> ### Heatmap of DEGs
> 
> `Scripts/Heatmaps.py` build the heatmaps for Fig2b, Fig2c, Fig3b, Fig3c. Gene order is set by hierarchical clustering performed on gene expression in synchronized cell lines.


> ### clustering of DEGs in synchronized cells
> 
> `Scripts/hclust_analysis_0.2_less_clusters.py` extracts the six clusters and build the boxplot panel shown in Fig2d and Fig3a. It is as an extension of [Heatmap of differentially expressed genes](#Heatmap-of-differentially-expressed-genes) section.


> ### classification of genes based on gene expression
> 
> `Scripts/second_decile_expression_genes_syncCells.r` define deciles of genes based on gene expression in synchronized cells.

<br/><br/>


## Poly-adenylation of replicative histone pre-mRNA

> ### count reads in replicative histone genes, or covering the stem loop (SL), or the histone downstream element (HDE)
> 
> `Scripts/reads_in_annot_job.sh` performs counting of reads mapped on a set of annotations. Annotations sets are `Data/annotations/histone_stem_loop.gtf` and `Data/annotations/histone_downstream_element.gtf`.  


> ### Evaluate the level of transcripts which are not cleaved in 3' of the SL
> 
> `Scripts/analysis_of_clived-unclived_replicative_histone_mRNA_Retriever.sh` count the number of transcript fragments cover a 49nt region that end at the HDE of each replicative histone gene.  
> 
> `Scripts/analysis_of_clived-unclived_replicative_histone_mRNA_Retriever.r` normalize the counts and build plots for Fig4b left. DESeq2 data were generated with `Scripts/DEG_analysis_ASF1only.r`


> ### Evaluate level of poly-adenylated transcripts
> 
> `Scripts/replHist_polyA_analysis.py` counts the number of read pairs carrying a poly-adelynation signature for each annotations. The input BAM file must have kept the unmapped reads, and the reads sorted by name.  
> 
> `Scripts/replHist_polyA_analysis.sh` orchestrate this counting on every sample and for several annotation sets.
> 
> `Scripts/replHist_polyA_analysis.r` yield the level of poly-adenylated transcripts for the set of replicative histones and build the panel in Fig4b right. This work uses the output of `replHist_polyA_analysis.py`.  


<!--  -->
