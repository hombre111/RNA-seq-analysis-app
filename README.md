# RNA-seq-analysis-app

The app is available as docker image at docker hub: hombre111/rna-seq_analysis_app

The app for analysis of bulk RNA-seq data made in shiny environment. Requires to upload a file with count matrix. The first column should contain Ensmbl gene IDs. The column names should be in following format: GROUP1, GROUP2 etc, where group corresponds to group name and number corresponds to individual sample. For example CONTROL1, CONTROL2, LPS1, LPS2 etc.

The app allows for quick exploratory data analysis, providing access to PCA results, correlation analysis of samples or easy plotting of gene expression with boxplots. The app also allows for differential gene expression analysis with edgeR and visualisation of results with ggplot.

So far, only annotations for rat genome- rnor6 are supported.
