# MAFTools

General Information and Setup

Code Involves the input of Text files titled TST_META.txt and TST_Mutation.txt. These Txt files were downloaded from Microsoft Excel, after data formatting and cleainng. It is important to note, that many columns within the text files act as dummy columns, and were not used within our final analysis. These columns include Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2. Such data inputs were combined and created into MAF object, titled MAF_TST170. TST_META.txt acted as clinical data, providing information on cancer variation, and unique patient sample ID. TST_Mutation.txt acted as primary data, which was analyzed. All corresponding files have been commited to the master branch of the MAFtools repository. 


Library Information

To conduct thorough analysis of the MAF_TST170 Object, multiple libraries were implemented to allow for thorough analysis. A list of libraries are included below, along with breif descriptions: 

1. library(maftools) - a package for analyzing and visualizing somatic mutation data in cancer genomes. It provides functions to import, summarize, visualize and annotate mutations from multiple sources.

2. library(tidyverse) - a collection of packages for data manipulation and visualization in R, including ggplot2, dplyr, tidyr, and others. It emphasizes a consistent and coherent set of data manipulation functions with a "grammar of data manipulation."

3. library(data.table) - a package that extends the basic functionality of R's data.frame with additional features for large datasets, such as fast aggregation and join operations, efficient sorting, and memory optimization.

4. library(ggplot2) - a package for creating high-quality graphics and visualizations in R. It is based on the grammar of graphics and provides a wide range of options for customizing plots.

5. library("mclust") - a package for clustering multivariate data using model-based methods. It provides functions for fitting mixture models with different types of covariance structures and for selecting the optimal number of clusters based on various criteria.


General Maftools Analysis

In terms of General Maftools Analysis, implementation of the following functions occurred to allow for effective data analyses and provide insight into the underlying causative mechansims of tumorigenesis and oncological research. 

1. Lollipop plots:
Lollipop plots are a type of graph that are commonly used to visualize genomic data, such as somatic mutations or copy number alterations. They consist of a horizontal line representing the genomic position of a gene or genomic region, with vertical lines (lollipops) indicating the presence and type of mutations or alterations at that position. The height of the lollipop represents the frequency or magnitude of the mutation or alteration.

Lollipop plots are useful for identifying recurrent mutations or alterations in a gene or genomic region and for visualizing their spatial distribution along the genome. They are often used in cancer genomics to identify driver mutations and to assess the clonal evolution of tumors.

Sample Code: 

lollipopPlot(
  MAF_TST170,
  gene = "TP53",
  AACol = 'Protein_Change',
  labelPos = NULL,
  labPosSize = 0.9,
  showMutationRate = TRUE,
  showDomainLabel = TRUE,
  cBioPortal = TRUE,
  refSeqID = NULL,
  proteinID = NULL,
  roundedRect = TRUE,
  repel = TRUE,
  collapsePosLabel = TRUE,
  showLegend = TRUE,
  legendTxtSize = 0.8,
  labPosAngle = 0,
  domainLabelSize = 0.8,
  axisTextSize = c(1, 1),
  printCount = FALSE,
  colors = NULL,
  domainAlpha = 1,
  domainBorderCol = "black",
  bgBorderCol = "black",
  labelOnlyUniqueDoamins = TRUE,
  defaultYaxis = FALSE,
  titleSize = c(1.2, 1),
  pointSize = 1.5
)

2. Oncoplots:
Oncoplots are a type of graph that are commonly used in cancer genomics to visualize the distribution and frequency of mutations or alterations across multiple genes or genomic regions. They consist of a grid of horizontal lines representing individual genes or genomic regions, with vertical bars indicating the presence and type of mutations or alterations at each position.

Oncoplots are useful for identifying patterns of co-occurring or mutually exclusive mutations or alterations across multiple genes or pathways, and for assessing their potential impact on cancer progression and treatment. They can also be used to stratify patients based on their genomic profiles and to identify potential targets for therapy.

Sample Code: 
oncoplot (MAF_TST170, top = 20,
         minMut = NULL,
         genes = NULL,
         altered = FALSE,
         drawRowBar = TRUE,
         drawColBar = TRUE,
         leftBarData = NULL,
         leftBarLims = NULL,
         rightBarData = NULL,
         rightBarLims = NULL,
         topBarData = NULL,
         logColBar = FALSE,
         includeColBarCN = TRUE,
         clinicalFeatures = NULL,
         annotationColor = NULL,
         annotationDat = NULL,
         pathways = NULL,
         selectedPathways = NULL,
         draw_titv = FALSE,
         showTumorSampleBarcodes = FALSE,
         barcode_mar = 4,
         barcodeSrt = 90,
         gene_mar = 5,
         anno_height = 1,
         legend_height = 4,
         sortByAnnotation = FALSE,
         groupAnnotationBySize = TRUE,
         annotationOrder = NULL,
         sortByMutation = FALSE,
         keepGeneOrder = FALSE,
         GeneOrderSort = TRUE,
         sampleOrder = NULL,
         additionalFeature = NULL,
         additionalFeaturePch = 20,
         additionalFeatureCol = "gray70",
         additionalFeatureCex = 0.9,
         genesToIgnore = NULL,
         removeNonMutated = TRUE,
         fill = TRUE,
         cohortSize = NULL,
         colors = NULL,
         bgCol = "#CCCCCC",
         borderCol = "white",
         annoBorderCol = NA,
         numericAnnoCol = NULL,
         drawBox = FALSE,
         fontSize = 0.8,
         SampleNamefontSize = 1,
         titleFontSize = 1.5,
         legendFontSize = 1.2,
         annotationFontSize = 1.2,
         sepwd_genes = 0.5,
         sepwd_samples = 0.25,
         writeMatrix = FALSE,
         colbar_pathway = FALSE,
         showTitle = TRUE,
         titleText = NULL,
)


3. Somatic interactions:
Somatic interactions refer to the complex interactions that occur between somatic mutations and alterations in cancer cells, including both genetic and epigenetic changes. These interactions can lead to the activation or inhibition of specific pathways and networks, and can have important implications for cancer progression and response to therapy.

Somatic interaction analysis involves the identification and characterization of these interactions, often using computational methods such as network analysis or pathway enrichment analysis. This can help to identify key driver mutations and pathways in cancer, as well as potential targets for therapy.

Sample Code: 
somaticInteractions(MAF_TST170, top = 20, pvalue = c(0.05, 0.1), fontSize = 0.55)

4. Oncogenic pathways:
Oncogenic pathways are cellular pathways that are dysregulated in cancer, leading to uncontrolled cell growth and proliferation. These pathways can be activated by a variety of genetic and epigenetic changes, including somatic mutations and alterations.

Oncogenic pathway analysis involves the identification and characterization of these dysregulated pathways, often using computational methods such as pathway enrichment analysis or network analysis. This can help to identify key driver pathways in cancer and to develop targeted therapies that specifically inhibit these pathways. It can also provide insights into the mechanisms underlying cancer progression and resistance to therapy.

Sample Code: OncogenicPathways(MAF_TST170, pathways = NULL, fontSize = 1, panelWidths = c(2, 4, 4))


Drug-Gene Interaction Database Analysis

Conducted Drug-Gene Interaction Database Analysis, to analyze the interaction between genes and specific therapeutic targets. The MAFtools package uses the Drug-Gene Interaction Database for clinical information. 

Sample Code:
dnmt3a.dgi = drugInteractions(genes = "TP53", drugs = TRUE)
dnmt3a.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]


Note: For additional Information, additional comments have been included throughout the code, to allow fo effective comprehension and organization of syntax. 


