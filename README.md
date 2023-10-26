# EnsembleFS: ensemble feature selection methods for analysis of molecular data

## Welcome to EnsembleFS
**EnsembleFS is a tool that allows the user to:**
* filter the most informative biomarkers from molecular data generated from high-throughput experiments that could be a new diagnostic/prognostic markers or therapeutic targets;
* establish the selected parameters for predictive models, such as the best number of most informative variables/biomarkers;
* remove redundant and correlated features from the obtained feature subsets;
* evaluate the quality of the feature set by using the random forest classification algorithm and machine learning validation techniques;
* offer various stability measures for assessing the stability of a given method on a given data set in comparison with other methods.
* establish the selected parameters for predictive models, such as the top N number of most informative variables;
* save and visualize the data and model results in the form of plots and tables;
* find information about gene collection (gene ontology, pathways, tissue specificity, miRNA targets, regulatory motif, protein complexes, and disease phenotypes) in several biological databases.

It can be applied to two-class problems. EnsembleFS is based on several filter feature selection algorithms, such as the U-test, the Monte Carlo FeatureSelection (MCFS), the MultiDimensional Feature Selection (MDFS), and the Minimum Redundancy Maximum Relevance (MRMR) for discovering the most important biomarkers and is used the machine learning algorithms to evaluate the quality of the set of variables. Predictive models are built using the Random Forest algorithm.
The information about each of the biomarkers was obtained from the biological databases, namely the Gene Ontology (molecular function, cellular component, biological process), the Kyoto Encyclopedia of Genes and Genomes (pathways), the Reactome (pathways), the WikiPathways (pathways), the Transfac (regulatory motif), the miRNA targets (miRTarBase), the Human Protein Atlas (tissue specificity), the CORUM protein complexes, and the Human Phenotype Ontology (human disease phenotypes).

**Tutorial**
[Watch the video](https://www.youtube.com/embed/ENf3LEMb56E)
