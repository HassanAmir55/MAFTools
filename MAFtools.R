library(maftools)
library(tidyverse)
library(data.table)
library(ggplot2)
library("mclust")

#READ IN TEXT FILES
TST_Mutation <- read_tsv("TST_Mutation.txt")
TST_Meta <- read_tsv("TST_Meta.txt")

#Convert to MAF File
MAF_TST170 <- read.maf(maf = TST_Mutation,clinicalData = TST_Meta, isTCGA = FALSE,)

#Developed TP53,BARD1,NOTCH1,CHEK2 and ATM 
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

lollipopPlot(
  MAF_TST170,
  gene = "BARD1",
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

lollipopPlot(
  MAF_TST170,
  gene = "NOTCH1",
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

lollipopPlot(
  MAF_TST170,
  gene = "CHEK2",
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

lollipopPlot(
  MAF_TST170,
  gene = "ATM",
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

#Plotted Lollipop Plots for other signifigant genes
lollipopPlot(
  MAF_TST170,
  gene = "AR",
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

lollipopPlot(
  MAF_TST170,
  gene = "ALK",
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

lollipopPlot(
  MAF_TST170,
  gene = "APC",
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


lollipopPlot(
  MAF_TST170,
  gene = "ARID1A",
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

lollipopPlot(
  MAF_TST170,
  gene = "MSH3",
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


lollipopPlot(
  MAF_TST170,
  gene = "PIK3CA",
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

lollipopPlot(
  MAF_TST170,
  gene = "PTCH1",
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

lollipopPlot(
  MAF_TST170,
  gene = "TET2",
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

lollipopPlot(
  MAF_TST170,
  gene = "BRCA1",
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

lollipopPlot(
  MAF_TST170,
  gene = "BRCA2",
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

lollipopPlot(
  MAF_TST170,
  gene = "BRIP1",
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

lollipopPlot(
  MAF_TST170,
  gene = "EGFR",
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

lollipopPlot(
  MAF_TST170,
  gene = "EP300",
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

lollipopPlot(
  MAF_TST170,
  gene = "KMT2A",
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

lollipopPlot(
  MAF_TST170,
  gene = "RB1",
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

lollipopPlot(
  MAF_TST170,
  gene = "TSC2",
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


#Development of oncoplot
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


#Developed Oncogene Pathways Data Visualization
OncogenicPathways(MAF_TST170, pathways = NULL, fontSize = 1, panelWidths = c(2, 4, 4))

#Plotted Oncogenic Pathways
PlotOncogenicPathways(MAF_TST170, pathways = "Cell_Cycle")
PlotOncogenicPathways(MAF_TST170, pathways = "RTK-RAS")
PlotOncogenicPathways(MAF_TST170, pathways = "PI3K")
PlotOncogenicPathways(MAF_TST170, pathways = "NOTCH")
PlotOncogenicPathways(MAF_TST170, pathways = "TP53")
PlotOncogenicPathways(MAF_TST170, pathways = "MYC")
PlotOncogenicPathways(MAF_TST170, pathways = "WNT")


#Developed Somatic Interactions Data Visualization
somaticInteractions(MAF_TST170, top = 20, pvalue = c(0.05, 0.1), fontSize = 0.55)

#Plotted Variant Allele Frequency (VAF)
plotVaf(MAF_TST170, vafCol = "Vaf_Percent")

#Developed Drug-Gene Interactions Visualization
dgi = drugInteractions(MAF_TST170, top = 5, fontSize = 0.75)
                
dnmt3a.dgi = drugInteractions(genes = "TP53", drugs = TRUE)
dnmt3a.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]

dnmt3a.dgi = drugInteractions(genes = "NOTCH1", drugs = TRUE)
dnmt3a.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]

dnmt3a.dgi = drugInteractions(genes = "BRCA1", drugs = TRUE)
dnmt3a.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]

dnmt3a.dgi = drugInteractions(genes = "BRCA2", drugs = TRUE)
dnmt3a.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]

dnmt3a.dgi = drugInteractions(genes = "EGFR", drugs = TRUE)
dnmt3a.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]

dnmt3a.dgi = drugInteractions(genes = "ALK", drugs = TRUE)
dnmt3a.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]

dnmt3a.dgi = drugInteractions(genes = "ROS1", drugs = TRUE)
dnmt3a.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]

dnmt3a.dgi = drugInteractions(genes = "BRAF", drugs = TRUE)
dnmt3a.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]