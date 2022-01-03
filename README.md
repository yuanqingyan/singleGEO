## SingleGEO: Query and Integrative Analysis of High-throughput Single Cell Sequencing Data from GEO

Single cell next generation sequencing technologies have advanced the research to examine the genomic alterations from individual cells. Single cell Assay for Transposase Accessible Chromatin with high-throughput sequencing (scATAC-seq) enables the research of chromatin-accessibility signatures across the entire genome. Single cell RNAseq (scRNAseq) studies the single-cell transcriptomics to better understand the differential cell compositions as wells as the molecular profile changes. Currently a huge amount of single cell genetic data have been deposited in GEO and the number is keeping on exponentially increasing. These single cell high-throughput sequencing data come from different tissues, health status, disease stage and drug treatments. Various machine learning algorithms, such as canonical correlation analysis and reciprocal principal component analysis, have been adapted in single cell data analysis. The implementation of such algorithms correct the potential batch effects and make the data integration between different datasets feasible. SingleGEO is developed to efficiently query, download high-throughput single cell sequencing data from GEO and perform and integrative analysis. The program will initially download GEO meta database into local drive and implement a fast and light-weight searching engine to query the database. By feeding the program with the keywords relevant to the user's research interests, the query results with specific meta information will be returned for further investigation. Once the targeted GSE IDs have been selected, the single cell sequencing data deposited in GEO supplement can be downloaded into local derive. In addition, singleGEO provides pipelines to easily generate Seurat objects from raw/normalized data and perform an integrative analysis with multiple datasets using Seurat framework. Novel findings can be obtained through the integrative analysis from different datasets. 

