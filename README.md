# Medaka-pituitary-scRNA-seq

## Background

The pituitary is a master endocrine gland in vertebrates, which controls a variety of physiological functions including growth, metabolism, homeostasis, reproduction, and response to stress. These functions are modulated by the secretion of several protein hormones, e.g. growth hormone (Gh), luteinizing hormone (Lh) and follicle-stimulating hormone (Fsh). In teleost fish (a group comprising nearly all ray-finned fishes), each hormone is presumably produced by a specific cell type. Unfortunately, key details of hormone production and its regulation are still poorly understood.  

## Data Analysis steps

#### Intended pipeline:

* CellRanger/STAR alignment (CR version 3.0.2, STAR 2.5.1b)
* DropEst quantification (version 0.8.5)
* Scone normalization
* Seurat exploration
* DoubletFinder doublet identification
* Final clustering
* Differentially expressed (DE) genes
* Ballon plot


#### Data:
* sample 1 = females
* sample 2 = males

#### Genome:
Oryzias_latipes.ASM223467v1.94.plusGFP_withchanges_13112018.gtf 
Oryzias_latipes.ASM223467v1.Ensembl94_plusGFP.fasta

#### Alignment:

Sample 1
```
/opt/cellranger-3.0.2/cellranger mkgtf genome13112018/Oryzias_latipes.ASM223467v1.94.plusGFP_withchanges_13112018 .gtf genome13112018/medaka_cellranger.gtf -- attribute=gene_biotype:protein_coding
/opt/cellranger-3.0.2/cellranger mkref --genome medaka_cellranger --fasta /media/sf_F_DRIVE/scmedaka/genome13112018/Oryzias_latipes.ASM223467v1.Ense mbl94_plusGFP.fasta --genes /media/sf_F_DRIVE/scmedaka/genome13112018/medaka_cellranger.gtf --nthreads 6
/opt/cellranger-3.0.2/cellranger count --id sample1_11072019 --fastqs /media/sf_F_DRIVE/scmedaka/fastq-files/sample1/ --transcriptome medaka_cellranger --chemistry threeprime --expect-cells 2500 --indices=SI- GA-C7
```

![average_plot_CTCF](https://github.com/sikh09/Medaka-pituitary-scRNA-seq/blob/master/Cell_ranger_results.png)

#### Quantification
```
dropest -f -g ../../genome13112018/medaka_cellranger/genes/ medaka_genes_fix_25112019.gtf -c 10x.xml –M -V -L eEBA -o sample1_25112019 ../outs/possorted_genome_bam.bam
```
output

* Unfiltered: 20419 genes, 74338 cells. Filtered: 18786 genes, 4945 cells.
* Exons: 17284 genes, 4945 cells.
* Introns: 15036 genes, 4945 cells.
* Intron/exon spanning: 13766 genes, 4945 cells.

#### computational analysis in R:
01_initial_analysis.R is used for initial quality control (QC) of scRNA seq data.
* Used DropEst to find the accurate estimation of moleular counts per cell. 

#### QC:
02_Remove_doublets.R is used to remove the douplets from scRNA seq data and updated quality control step.
* I ended up with 2592 cells, 18377 genes in sample 1
* I ended up with 3804 cells, 18660 genes in sample 2

#### Downstream Analysis:
3_Downstream_analysis.R is used for normalization, initial clustering and sub-clustering, final clustering,DE gene, heatmap and ballon plot

