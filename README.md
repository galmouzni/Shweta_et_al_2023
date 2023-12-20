# Shweta_et_al_2023
Regulation of replicative histone RNA metabolism by the histone chaperone ASF1  

This Repository is under license GNU-PLv3.

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

`env_bin.sh` define path to useful executables

<br/><br/>


## Differential Gene Expression analysis

> ### Expression of histone, replicative and variant, genes
> 
> `Scripts/expression_analysis.r` build the heatmap of expression relative to mean level in siCtrl samples.


> ### Evaluate differentially expressed genes (DEGs)
> 
> `Scripts/dge_analysis_0.2.py` perform differential expression analysis and find the DEGs.


> ### Volcano plot of genes in asynchronous cells
> 
> `Scripts/Volcano_plot.r` build the volcano plot of Fig2a and SuppFig2b.


> ### Heatmap of DEGs
> 
> `Scripts/Heatmaps.py` build the heatmaps for Fig2c, Fig2d, Fig3b, Fig3c, SuppFig3a and SuppFig3d. Gene order is set by hierarchical clustering performed on gene expression in synchronized cell lines.


> ### clustering of DEGs in synchronized cells
> 
> `Scripts/hclust_analysis_0.2_less_clusters.py` extracts the six clusters shown in Fig2c and Fig3b, and build the boxplot panels shown in Fig2e, Fig3a and SuppFig3e. It is as an extension of [Heatmap of differentially expressed genes](#Heatmap-of-differentially-expressed-genes) section.


> ### classification of genes based on gene expression
> 
> `Scripts/second_decile_expression_genes_syncCells.r` define deciles of genes based on gene expression in synchronized cells.

<br/><br/>


## Poly-adenylation of replicative histone pre-mRNA

> ### count reads in replicative histone genes, or covering the stem loop (SL), or the histone downstream element (HDE)
> 
> `Scripts/reads_in_annot_job.sh` performs counting of reads mapped on a set of annotations. Annotations sets are `Data/annotations/histone_stem_loop.gtf` and `Data/annotations/histone_downstream_element.gtf`.  


> ### Signal in ORF and in region downstream of the SL
> 
> `Scripts/SpikeInControls_ORFdsSLsignal_remappedPolyAreads_pubHeLaExpression.r` normalize by spike-in controls and sequencing depth the counts, and build plots for Fig3d and SuppFig5c.
> 
> `Scripts/metagene_analysis.r` normalize and plots metagene foldchange for Fig4a.


> ### Evaluate the level of transcripts which are not cleaved in 3' of the SL
> 
> `Scripts/analysis_of_clived-unclived_replicative_histone_mRNA_Retriever.sh` count the number of transcript fragments that cover a 49nt region which end at the HDE of each replicative histone gene.  
> 
> `Scripts/analysis_of_clived-unclived_replicative_histone_mRNA_Retriever.r` normalize the counts and build plots for SuppFig5b. DESeq2 data were generated with `Scripts/expression_analysis.r`


> ### Evaluate level of poly-adenylated transcripts
> 
> `Scripts/replHist_polyA_analysis.py` counts the number of read pairs carrying a poly-adelynation signature for each annotations. The input BAM file must have kept the unmapped reads, and the reads sorted by name.  
> 
> `Scripts/replHist_polyA_analysis.sh` orchestrate this counting on every sample and for several annotation sets.
> 
> `Scripts/replHist_polyA_analysis.r` yield the level of poly-adenylated transcripts for the set of replicative histones and build the panel SuppFig5b. This work uses the output of `replHist_polyA_analysis.py`.  

<br/><br/>


## Statistical analyses of RT-qPCR results

`Scripts/RT-qPCR_statistical_tests.r` and `Scripts/RT-qPCR_statistical_tests_oriented.r` compute the statistical test for every condition from every RT-qPCR experiment compared to an estimated reference distribution, or between conditions. The two scripts respectively perform "two.sided" comparisons and considering "lower than" or "greater than" comparison according to the needs in the study; in Fig1e and Fig4b.


<br/><br/>

## In Data

> `raw.csv.gz` are the raw read counts per gene and per sample
> `design.csv` is the design of the experiment
> `GRCh38.104.sorted.geneOnly.gtf` are the gene annotations downloaded from EnsEMBL
> `histones.gff` is derived from `GRCh38.104.sorted.geneOnly.gtf`, containing a unique reference annotation for every histone gene (replicative and other variants).
> 
> `annotations` folder contains the annotations used in this study:
> - `genes.GRCh38.95.csv` are the gene annotations downloaded from EnsEMBL first used for processing the RNA-seq data.
> - `histone_downstream_element.gtf` are the HDE annotations established by Sébastien Lemaire in this study for all the replicative histone.
> - `histone_stem_loop.gtf` are the SL annotations established by Sébastien Lemaire in this study for all the replicative histone.
> - `replicative_histone_genes.gff` derived from `genes.GRCh38.95.csv` for annotations of the replicative histone genes only.
> - `replicative_histone_CDS.gtf` is derived from `genes.GRCh38.95.csv` for CDS annotations of the replicative histone genes.
> - `variant_histone_and_histone_chaperone_genes.gff` is derived from `genes.GRCh38.95.csv` for CDS annotations of the non-replicative histone genes and of the genes encoding histone chaperones.
> 
> `ERCC92` folder for spike-in controls related files:
> - `ERCC_Controls_Analysis.txt` and `ERCC92.csv` are metadata about spike-in controls
> - `ERCC92.fa` are the sequences of the spike-in controls.
> 
> `GtexTcgaGeneExpr`:
> `histone_chaperone_complete_list.csv` is a table of all the histone genes (replicative and other variants) and of the histone chaperone genes.


<br/><br/>


## In Results

> `sequencing_depth_norm_factor.tsv`, `spike_norm_factor.tsv`, `spike_and_seqDepth_norm_factor.tsv` are the normalization factors respectively based on sequencing depth alone, spike-in controls alone, and the combination of the two.


<!--  -->
