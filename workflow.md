This Rbbt workflow is in charge of finding the steady state of cell-lines. The steady state is represented as a pattern of activations of protein activities.

##Data sources

The data used for determining the steady states of cell-lines comes from a few
different resources that describe them in terms of molecular characteristics.
These are the resources along with the data that can potentially consider:

* CCLE: Cancer cell-line encyclopedia (CNV, expression)
* GDSC: Genomics of Drug Sensitivity in Cancer (CNV, expression)
* Achilles: gene essentially in cell lines (gene essentiality)
* MCLP: MD Anderson Cell-Line Project (phospho-proteomics)
* COSMIC: Catalogue of Somiatic Mutations In Cancer (somatic mutations)
* COREAD: Colo-rectal Adenocarcinoma from Roumeliotis (MS-phospho) 

These resources each deal with a subset of all possible cell-lines. There are
large overlaps, but not all cell-lines are in all resources. For instance only
around 700 cell-lines are in common between MCLP and CCLE. Since each datasets
refers to the cell-lines in slightly different ways, we use a term
normalization engine to match the name of the cell lines that is queried to the
name in each database.

##Approach

The different resources provide with different molecular characterizations of
the cell-lines. From these we need to determine the activity of the different
proteins. To do this we make use of the Paradigm software tool

###Paradigm

Paradigm takes in a topology of relationships between molecular entities and
across the different levels of the central dogma. The central dogma implies
that CNV affect expression, expression affects protein abundance, protein
abundance affects protein activity. Beyond the central dogma the topology
establishes relationships across these entities, such as that MDMA abundance
inhibits TP53 activity, or that FOXA1 activity promotes expression of its
target genes.

Paradigm takes this topology and builds a discrete Bayesian graph which is
imputed with a series of observations that we need to provide it with. These
observations refer to the different levels of the central dogma. That is, they
can be CNV, gene expression, protein abundance or protein activity
observations, which we derive from the different data sources.

####Paradigm topology

We construct the topology that we feed to paradigm by first considering the
core topology of interest, which in our work represents the relevant
cell-signaling around the cell cycle control. In order to better exploit the
gene expression data we introduce the necessary components to link the
signaling proteins to transcription factors and these to their target genes.
This allows us to infer the activity of transcription factors out of the
expression levels of the target genes. However, we do not actually use the
Paradigm topology to link the transcription factors to the target genes, as
this will increase tremendously the size of the topology, instead we will
separately calculate the activity of the transcription factors using a
different tools and introduce that as observations to Paradigm, which will then
be able to use those activity to infer the activities of the signal
transduction proteins that regulate these.

####Calculating Transcription Factor activity

We take the network of transcription factors and target genes, which we call
the *regulome*, and use it in the TF-activity inference tool of our choice
(either Viper or ROMA). This regulome is derived form the ExTRI tool (formerly
FNL), which is a resource that combines DBs of transcription regulation with
text-mined regulation events. From ExTRI we take all TF-TG pairs that are
listed in at least two databases with high-quality (see the section on ExTRI
for more information). The entire mRNA matrix for all cell-lines is
run through the inference tool to be run collectively, rendering activity
vectors for TFs for each cell-lines.

####Paradigm Observation data

As we mentioned above, Paradigm needs to ground his Bayesian graph with some
observation data to make its inferences. It can take observations at all the
different levels of the central dogma. In CLSS we can use CNVs for the genome
level, mRNA expression for the expression level, and TF activity and MCLP RPPA
protein activity for the activity level. One does not have to include all
levels in the Paradigm tool. We provide two versions, one that only uses CNV
and mRNA, and one that include TF-activity as well. One could easily explore
including RPPA for protein activity and also for protein abundance (it measures
both)

Note that for mRNA used as observation data we perform a normalization that
transforms the original values into estimates of high and low expression. To this
end we use the R method 'densityMclust' from the 'mclust' package.

In the current implementation we use CCLE data for CNV and mRNA expression
(used directly as observation data and as input to Viper to form TF activities)

