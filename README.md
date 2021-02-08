# WebEFS: ensemble feature selection methods for analysis of molecular data

## Welcome to WebEFS-tools
** WebEFS is tool that allows the user to: **
* filter the most informative biomarkers from molecular data generated from high-throughput microarray experiments that could be a new diagnostic/prognostic markers or therapeutic targets;
* establish the selected parameters for predictive models, such as the best number of most informative variables/biomarkers;
* examine the impact of correlation between informative features on the predictive power of predictive model;
* evaluate stability of biomarker sets and performance of predictive models;
* find information about gene collection (gene ontology, pathways, tissue specificity, miRNA targets, regulatory motif, protein complexes, disease phenotypes) in several biological databases.

It can be applied to two-class problems. WebEFS based on the several fil-ter feature selection algorithms, such as the U-test, the Monte Carlo FeatureSelection (MCFS), the MultiDimensional Feature Selection (MDFS) and theMinimum Redundancy Maximum Relevance (MRMR) for discovering the mostimportant biomarkers and used the machine learning algorithms to evaluatequality of the set of variables. Predictive models are built using the Random Forest algorithm.
The information about each of the biomarkers was obtained from the biological databases, namely Gene Ontology (molecular function, cellular component, biological process), Kyoto Encyclopedia of Genes and Genomes (pathways), Reactome (pathways), WikiPathways (pathways), Transfac (regulatory motif), miRNA targets (miRTarBase), Human Protein Atlas (tissue specificity), CORUM protein complexes, and Human Phenotype Ontology (human disease phenotypes).
