# EnsembleFS web-based tool for a filter ensemble feature selection of molecular omics data
## Project developed by:
Aneta Polewko-Klim and Paweł Grablis

**Citing EnsembleFS in publications and presentations**

Polewko-Klim, A., Grablis, P., Rudnicki, W. (2024). EnsembleFS: an R Toolkit and a Web-Based Tool for a Filter Ensemble Feature Selection of Molecular Omics Data. In: Franco, L., de Mulatier, C., Paszynski, M., Krzhizhanovskaya, V.V., Dongarra, J.J., Sloot, P.M.A. (eds) Computational Science – ICCS 2024. ICCS 2024. Lecture Notes in Computer Science, vol 14835. Springer, Cham. https://doi.org/10.1007/978-3-031-63772-8_7


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

## Web demo server
[https://uco.uwb.edu.pl/apps/EnsembleFS/](https://uco.uwb.edu.pl/apps/EnsembleFS/)

[https://ensemblefs.shinyapps.io/ensemblefs/](https://ensemblefs.shinyapps.io/ensemblefs/)


## Tutorial
[Watch the video](https://www.youtube.com/embed/ENf3LEMb56E)

## Software project documentation

A software user manual: 
[UserDocumentationEnsembleFS.pdf](https://github.com/biocsuwb/EnsembleFS/blob/main/User%20Documentation%20EnsembleFS.pdf)

## Main functionality modules of EnsembleFS web app: Feature Selection tab and Gene information tab

![Fig.1](https://github.com/biocsuwb/Images/blob/main/S2.jpg?raw=true)
Fig.1 Main functionality of the Feature Selection module of EnsembleFS web app.

![Fig.2](https://github.com/biocsuwb/Images/blob/main/S3.jpg?raw=true)
Fig.1 Main functionality of the Gene Information module of EnsembleFS web app.

