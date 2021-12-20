#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)
library(ensembleFS)
library(grid)
library(gridExtra)
library(venn)
library(gprofiler2)
library(plotly)
library(ggplot2)
library(slickR)

ui <- fluidPage(theme = shinytheme("sandstone"),
                titlePanel(title=div(img(src="banner.gif",  align = "center", height="60%", width="100%"))),
                navbarPage("", collapsible = T, id="demo",
                           tabPanel(h4('HOME'), icon = icon("home", lib = "glyphicon", "fa-2x"),
                                    mainPanel(width = 11,
                                              h1(strong('Welcome to  EnsembleFS ')),
                                              hr(),
                                              div(
                                                h4("EnsembleFS is tool that allows the user to: "),
                                                h4('- filter the most informative features (biomarkers) by using ensemble feature selection approach from molecular data generated from high-throughput molecular 
                                   biology experiments;'),
                                                h4('- establish the selected parameters for predictive models, such as the number of top N informative features;'),
                                                h4('- remove redundant features by building a the Spearman correlation matrix that identifies highly correlated features;'),
                                                h4('- evaluate the stability of feature subsets and performance of predictive models;'),
                                                h4('- find information about selected biomarkers (gene ontology, pathways, tissue specificity, miRNA targets, regulatory motif, protein complexes, disease phenotypes) 
                                   in several biological databases.'),
                                                h4('It can be applied to two-class problems. EnsembleFS is based on several filter feature selection algorithms, such as the test U Manna-Whitneya, 
                                   the Monte Carlo Feature Selection (MCFS), the MultiDimensional Feature Selection (MDFS), and the Minimum Redundancy Maximum Relevance (MRMR) for discovering the most important biomarkers and used the machine learning algorithms to evaluate the quality of feature sets. Predictive models are built using the Random Forest algorithm.'),
                                                h4('The information about each of the biomarkers is extracted from diverse biological databases, namely the Gene Ontology, the Kyoto Encyclopedia of Genes and Genomes, 
                                   the Reactome, the WikiPathways, the Transfac, the miRNA targets, the Human Protein Atlas, the CORUM, and the Human Phenotype Ontology.'),
                                                h4(strong('For details on feature selection and classification algorithms, please refer to the Help -> Terminology')),
                                                
                                                h1(strong('Overview')),
                                                hr(),
                                                img(src="Overview_func.png", align = "center",height='900px',width='2200px'),
                                                h4('Figure 1. Flow chart for EsembleFS: A) the feature selection process and model evaluation process, B) a scheme of ensemble-based feature selection method,
                                   C) a scheme for biological information collection and integration about biomarkers. For details on used notation, please refer to the Help -> Terminology'),
                                                
                                                h3(strong('The feature selection process and the model validation')),
                                                hr(),
                                                h4('The FS process is based on the combination of the heterogeneous ensemble approach for feature selection and machine learning. As shown in Figure 1A, 
                                   the basic modeling procedure involves the following steps:'),
                                                h4('(1) split randomly the samples into training and test set;'),
                                                h4('(2) select and rank informative features using the feature filter on training set;'),
                                                h4('(3) remove correlated features with the training set;'),
                                                h4('(4) build model on the training set;'),
                                                h4('(5) estimate the stability of a FS algorithm and quality of the predictive model on test set.'),
                                                h4('To evaluate quality of the feature set, EnsambleFS applies the supervised machine mearning procedure, namely, the predictive models are built with top N features
                                using the random forest algorithm and the stratified 3-fold cross-validation procedure or the 0.3 random sampling is used to evaluate classification models. 
                                The mean values of ACC, AUC, MCC, and ASM metric are calculated for each of FS methods. To minimize the collinearity of features, the correlated variables are removed 
                                using the correlation coefficient method from the train dataset. The most informative features are identified as those that appear most consistently in top N features 
                                selected by FS method in n resampling operations. Herein, the best feature set is selected by the majority voting method (ie. morethan half) for each basic feature 
                                selector. As shown in Figure 1B, the total feature set for further biological analysis included all the best features or overlapping best features selected by basic filters.'),
                                                hr(),
                                                h4('For more details on the used procedure for building predictive models, please refer to [1]'),
                                                h5('[1] A. Polewko-Klim, W.R. Rudnicki. Analysis of Ensemble Feature Selection for Correlated High-Dimensional RNA-Seq Cancer Data. In: Krzhizhanovskaya V. et al. (eds) Comp.Scien. ICCS 2020. Lecture Notes in Computer Science 12139 (2020), 525-538, Springer, Cham.'),
                                                uiOutput("url.art.EnsembleFS"),
                                                
                                                h3(strong('The biological gene information collection')),
                                                hr(),
                                                h4('EnsembleFS allows the users to access to fundamental biological information for finally selected biomakers from several databases. 
                                As shown in Figure 1C, user research analysis may include intersection or union of the most informative biomarker sets with up to five FS methods.
                                For each biomarker, the biological information is obtained from the following databases: the Gene Ontology (the molecular function, the cellular component, 
                                and the biological process), the Kyoto Encyclopedia of Genes and Genomes (the pathways), the Reactome (the pathways), 
                                and WikiPathways (the pathways), Transfac (the regulatory motif), miRTarBase (the miRNA targets), 
                                the Human Protein Atlas (the tissue specificity), the CORUM (the protein complexes), 
                                and the Human Phenotype Ontology (the human disease phenotypes).'),
                                                
                                                
                                                h1(strong('Workflow')),
                                                hr(),
                                                img(src="Workflow.png", align = "center",height='500px',width='2200px'),
                                                h4('Figure 2. Main functionality modules of EnsembleFS: A) Feature Selection tab, B) Gene information tab. The cuboids represent the interaction between EnsembleFS and the user, and the ellipses represent EnsembleFS processes.'),
                                                
                                                h1(strong('Tutorial')),
                                                hr(),
                                                h4(strong('For examples on how to use EnsembleFS, please refer to the Help -> Tutorial')),
                                                
                                                h3(strong('Cite as')),
                                                hr(),
                                                h4('When using this web server, please cite the following references:'),
                                                
                                                h3(strong('EnsembleFS R toolkit')),
                                                hr(),
                                                uiOutput("home.url_app"),
                                                uiOutput("home.url_package"),
                                                
                                                h3(strong('Contact')),
                                                hr(), 
                                                uiOutput("home.email"),
                                                hr(), 
                                                h4('Developed by Pavel Hrablis & Aneta Polewko-Klim'),
                                                h4('Institute of Computer Science, University of Bialystok, Bialystok, Poland')
                                              ),
                                    )),
               tabPanel(h4('FEATURE SELECTION'),icon = icon("hand-pointer", "fa-2x"),
                        sidebarPanel(
                            textOutput('info.load.main.data'),
                            fileInput('file1', 'Load file',
                                      multiple = FALSE,
                                      accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                            
                            
                            checkboxInput('data.default', label = h4('Load demo data'), value = FALSE),
                            
                            numericInput("num",
                                         label = h4("Column number of decision variable"),
                                         value = 1, min = 1),
                            
                            checkboxGroupInput("methods",
                                               label = h4("Feature selection methods (FS)"), 
                                               choices = list("U-test" = 'fs.utest',
                                                              "MRMR" = 'fs.mrmr',
                                                              "MCFS" = 'fs.mcfs',
                                                              "MDFS-1D" = 'fs.mdfs.1D',
                                                              "MDFS-2D" = 'fs.mdfs.2D'),
                                               selected = 'fs.utest'),
                            
                            conditionalPanel(
                                condition = "input.methods.includes('fs.utest') ||  input.methods.includes('fs.mdfs.1D') || input.methods.includes('fs.mdfs.2D')",
                                selectInput("adjust", label = h4("Multitest correction"),
                                             choices = list("Benjamini & Hochberg" = "BH",
                                                            "Benjamini & Yekutieli" = "BY",
                                                            "Bonferroni" = "bonferroni",
                                                            "FDR" = "fdr",
                                                            "Hochberg" = "hochberg",
                                                            "Holm" = 'holm',
                                                            "None" = "none"), 
                                             selected = 'fdr')
                            ),
                            
                            
                            conditionalPanel(
                              condition = "input.methods.includes('fs.mrmr')",
                              numericInput("nvar", label = h4("Number of relevant variables"), value = 100),
                              helpText("Note: hyperparameter of MRMR,
                                         setup no more than the number of all variables")
                            ),
                            
                            
                            conditionalPanel(
                              condition = "input.methods.includes('fs.mcfs')",
                              selectInput("cutoff", label = h4("Cutoff method"),
                                          choices = list("kmeans" = "kmeans",
                                                         "permutations" = "permutations",
                                                         "criticalAngle" = "criticalAngle"), 
                                          selected = 'kmeans'),
                              helpText("Note: hyperparameter of MCFS,
                                          The methods of finding cutoff value between important and unimportant attributes")
                            ),
                            
                            
                            sliderInput("level.cor", label = h4("Correlation coefficient"), min = 0, max = 1, value = 0.75),
                            
                            radioButtons("cv", label = h4("Validation methods"),
                                         choices = list("random sample (test set 30%)" = 'rsampling', "3-fold cross-validation" = 'kfoldcv'), 
                                         selected = 'rsampling', inline = TRUE),
                            
                            sliderInput("niter", label = h4("Number iteration"), min = 1, max = 30, value = 10),
                            
                            uiOutput("run"),

                        ),
                        mainPanel(
                            
                            tabsetPanel(
                                type = "tabs",
                                tabPanel("Data",  dataTableOutput('contents')),
                                tabPanel("Ranking List", dataTableOutput("ranking"), uiOutput('downloadRanking')),
                                tabPanel("FS Stability",  dataTableOutput("stability"), uiOutput('downloadAsm')),
                                tabPanel("Model Accuracy",  dataTableOutput("model"), uiOutput('downloadModel')),
                                tabPanel("Plots",
                                         plotlyOutput("plot.stab"),
                                         plotlyOutput("plot.model.acc"),
                                         plotlyOutput("plot.model.auc"),
                                         plotlyOutput("plot.model.mcc"),
                                         br(),
                                         uiOutput('downloadPlot')),
                                tabPanel("Download Zip", br(), uiOutput('downloadZip'))
                            ))),
               
               tabPanel(h4('GENE INFORMATION'), icon = icon("dna", "fa-2x"),
                        sidebarPanel(
                            textOutput('info.load.data'),        
                            hr(), 
                            h4("Number of relevant biomarkers"),
                            hr(),
                            sliderInput("geneNumber",
                                        label = "Top N variables with FS filter",
                                        min = 1, max = 100 , value = 100),
                                        #choices = list(5,10,15,20,30,40,50,75,100), selected = 100),
                            helpText("Note: biomarkers are analyzed that occur in at least half of the feature subsets"),
                            hr(),
                            uiOutput('checkbox.condition.method'),
                            hr(),
                            uiOutput('get.information'),
                            hr(),
                            radioButtons("typeBase", label = "Data bases",
                                         choices = list("Molecular function (GO:MF)" = 1, 
                                                        "Cellular component (GO:CC)" = 2, 
                                                        "Biological process (GO:BP)" = 3,
                                                        "KEGG" 	= 4,
                                                        "Reactome (REAC)" = 5,	
                                                        "WikiPathways (WP)" = 6,
                                                        "Transfac (TF)" = 7,
                                                        "miRTarBase (MIRNA)" = 8,
                                                        "Human Protein Atlas (HPA)" = 9,
                                                        "CORUM protein complexes" = 10,
                                                        "Human Phenotype Ontology (HP)" = 11,
                                                        "All" = 12),
                                         selected = 1),

                            hr(),
                            uiOutput('save.information'),
                            hr()
                            
                        ),
                        
                        mainPanel(
                            fluidRow(splitLayout(cellWidths = c("50%", "50%"),
                                 plotOutput("graph.venn.methods"), 
                                 plotOutput("graph.venn.result" ))  
                            ),
                            dataTableOutput('information')
                        )),
               tabPanel(h4(' HELP '), icon = icon("question-circle", "fa-2x"), 
                        tabsetPanel(
                          type = "tabs",
                          tabPanel('Terminology', 
                                   h1(strong('Feature selection algorithms')),
                                   hr(),
                                   h3(strong('Mann-Whitney U-test (U-test)')),
                                   h4('U-test is a nonparametric statistical test [1] that assigns a probability to the hypothesis that two samples corresponding'),
                                   h4('to two decision classes (e.g. normal and tumor tissue) are drawn from populations with the same average value.'),
                                   hr(),
                                   h5('[1] H. B. Mann and D. R. Whitney. Controlling the false discovery rate: A practical and powerful approach to multiple testing. Ann. Math. Statist., 18(1):50–60, 1947.'),
                                   h3(strong('MultiDimensional Feature Selection (MDFS)')),
                                   h4('MDFS method measures the decrease of the information entropy of the decision variable due to knowledge of k-dimensional tuples of variables and measures the influence of each variable in the tuple.'),
                                   h4('It performs an exhaustive search over all possible k-tuples and assigns to each variable a maximal information gain due to a given variable that was achieved in any of the k-tuple that included this variable.'),
                                   h4('MDFS-1D is one-dimensional version of this algorithm. MDFS-2D is two-dimensional version of this algorithm.[2]'),
                                   hr(),
                                   h5('[2] K. Mnich and W.R. Rudnicki. All-relevant feature selection using multidimensional filters with exhaustive search. Information Sciences, 524:277-297, 2017'),
                                   uiOutput("url.MDFS"),
                                   h3(strong('Monte Carlo Feature Selection (MCFS)')),
                                   h4('MCFS method relies on a Monte Carlo approach to select informative features. The MCFS-ID algorithm is capable of incorporating inter-dependencies between features.'),
                                   h4('The EnsambleFS uses the one of three different cut-off methods for discerning informative features and non-informative features, namely critical angle, permutations, and k-means.[3]'),
                                   h4('For example k-means method (the default cut-off method) clusters the relative importance values into two clusters and sets the cutoff where two clusters are separated.'),
                                   hr(),
                                   h5('[3] M. Draminski, M. Kierczak, J. Koronacki, and J. Komorowski. Monte Carlo Feature Selection and Interdependency Discovery in Supervised Classification, Advances in Machine Learning II. Studies in Computational Intelligence, vol 263. Springer, Berlin, Heidelberg'),
                                   uiOutput("url.MCFS"),
                                   h3(strong('Minimum Redundancy Maximum Relevance (MRMR)')),
                                   h4('The MRMR method is based on mutual information as a measure of the relevancy and redundancy, where the redundancy of a selected feature subset is an aggregate mutual information measure'),
                                   h4('between each pair of features in the selected feature and the relevancy to a class is an aggregate mutual information measure between each feature with respect to the class.[4]'),
                                   hr(),
                                   h5('[4] C. Ding and H. Peng. Minimum redundancy feature selection from microarray gene expression data. Journal of Bioinformatics and Computational Biology, 3(2):185–205, 2005.'),
                                   uiOutput("url.MRMR"),
                                   h1(strong('Classification algorithm')),
                                   hr(),
                                   h3(strong('Random Forest')),
                                   h4('Random forest is an ensemble of decision trees, where each tree is built on a different bagging sample of the original data set.[5] For each split, a subset of variables is selected randomly '), 
                                   h4('and the one is selected that allows to achieve the highest Gini coefficient for the resulting leaves. Random Forest works well on data sets with a small number of objects, '),
                                   h4('has few tunable parameters that do not relate directly to the data, and very rarely fails.'),
                                   hr(),
                                   h5('[5]  L Breiman. Random forests. Machine Learning, 45:5-32, 2001.'),
                                   uiOutput("url.RF"),
                                   h1(strong('Model evaluation metrics')),
                                   hr(),
                                   h3(strong('The area under receiver operator curve (AUC)')),
                                   h4('AUC - ROC curve is a performance measurement for the classification problems at various threshold settings. ROC is a probability curve and AUC represents the degree or measure of separability.'),
                                   hr(),
                                   h3(strong('Accuracy (ACC)')),
                                   h4('The accuracy of a machine learning classification algorithm is one way to measure how often the algorithm classifies a data point correctly.'),
                                   h4('This evaluation metric is the number of correctly predicted data points out of all the data points. The ACC values range from 0 (no correct decision/prediction) to 1 (perfect decision/prediction)'),
                                   h4('The accuracy is expressed by:'),
                                   img(src="ACC.png", align = "center",height='120px',width='950px'),
                                   hr(),
                                   h3(strong('Matthews correlation coefficient (MCC)')),
                                   h4('The Matthews correlation coefficient or phi coefficient is used in machine learning as a measure of the quality of binary classifications. The MCC value between -1 and +1.'),
                                   h4('A coefficient of +1 represents a perfect prediction, 0 an average random prediction and -1 an inverse prediction. The MCC is expressed by formula:'),
                                   img(src="MCC.png", align = "center",height='60px',width='950px'),
                                   h1(strong('Feature stability measure')),
                                   hr(),
                                   h3(strong('The Lustgarten adjusted stability measure (ASM)')),
                                   h4('The  total stability of feature selection method is  measured  as  the  average  of  the  pairwise  similarity  for  all  pairs'),
                                   h4('of the most informative feature subsets from n runs of a model. For this, the ASM metric was used, that is described by the following formula:'),
                                   img(src="ASM.png", align = "center",height='150px',width='950px'),
                                   hr(),
                                   h5('[6] J.L. Lustgarten, V. Gopalakrishnan, and S. Visweswaran, Measuring stability of feature selection in biomedical datasets, AMIA ... Annual Symposium proceedings, AMIA Symposium vol. 2009 406-10, 14 Nov. 2009'),
                                   uiOutput("url.Lustgarten"),
                                   hr()
                          ),
                          tabPanel('Tutorial', 
                                   h1(strong('YouTube: EnsembleFS tutorial')),
                                   hr(),
                                   HTML('<iframe width="800" height="400" src="https://www.youtube.com/embed/e5oHiaigA68" frameborder="0" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen"></iframe>'),
                                   h1(strong('A short workflow')),
                                   hr(),
                                   
                                   h3(strong('The Ensemble feature selection process')),
                                   h4('The process of the selection of the most informative biomarkers includes the following steps:'),
                                   h4('1. Navigate to the FEATURE SELECTION tab.'),
                                   h4('2. Load csv/txt file (separator = ";", decimal = ",") and enter the column number for the binary decision variable'),
                                   h4('3. Choose the following parameters:'),
                                   h4('- Feature selection methods: MRMR, MCFS, MDFS-1D'),
                                   h4('- Multitest correction: fdr'),
                                   h4('- MRMR parameter: number of significant features: 100'),
                                   h4('- MCFS parameter: cut-off method: k-means'),
                                   h4('- Correlation coefficient: 0.75'),
                                   h4('- Validation methods: 3-fold cross-validation'),
                                   h4('- Number of repetitions: 10'),
                                   h4('4. Press RUN FEATURE SELECTION'),
                                   h4('5. Navigate to RANKING LIST tab for the set of most informative biomarkers'),
                                   h4('6. Navigate to FS STABILITY tab for stability calculation results'),
                                   h4('7. Navigate to MODEL ACCURACY tab for the model building results'),
                                   h4('8. Navigate to PLOT to visualize stability results and model building results'),
                                   h4('9. Navigate to DOWNLOAD ZIP tab to download all results as one archive'),
                                   
                                   h3(strong('Searching biological information about biomarkers')),
                                   h4('The process of aggregating information about the most informative biomarkers includes the following steps:'),
                                   h4('1. Navigate to GENE INFORMATION tab for biological information on the top biomarkers'),
                                   h4('2. Number of relevant biomarkers:'),
                                   h4('- Top N variables with FS filter: 100'),
                                   h4('- Combination of a set of biomarkers: union'),
                                   h4('- Data bases: all'),
                                   h4('3. Press GET ANALYSIS')
                                   
                          ),
                          tabPanel('Example', 
                                   h1(strong('Study experimental dataset')),
                                   hr(),
                                   h5('The RNA-sequencing data of tumor-adjacent normal tissues of lung adenocarcinoma cancer patients.'),
                                   uiOutput("data.sample.set"),
                                   h5('The preprocessing of data involved standard steps for RNA-Seq data. The log2 transformation was performed.'), 
                                   h5('Features with zero and near-zero (1%) variance across patients were removed.'),
                          )
                        )),
               tabPanel(h4('Licence Information'), icon = icon("file-alt", "fa-2x"),
                        h2('Availability and requirements:'),
                        hr(),
                        h4(strong('Project name:')),
                        h4('EnsembleFS: ensemble feature selection methods for analysis of molecular data'),
                        hr(),
                        h4(strong('Source code:')),
                        h4('https://github.com/biocsuwb/EnsembleFS'),
                        hr(),
                        h4(strong('Operation system(s):')),
                        h4('Web based, Platform independent'),
                        hr(),
                        h4(strong('Programming language:')),
                        h4('R, R SHINY'),
                        hr(),
                        h4(strong('Other requirements:')),
                        h4('Modern Browser'),
                        hr(),
                        h4(strong('License:')),
                        h4('BSD License'),
                        hr(),
                        h4(strong('Any restrictions to use by non-academics:')),
                        h4('None.'),
                        hr(),
                        h4(strong('Source code is available from GitHub (https://github.com/biocsuwb/EnsembleFS) and is free open source software under an MIT license.'))
               )
               
                ))

server <- function(session, input, output){
    
  options(shiny.maxRequestSize=200*1024^2)
  
  store.result <- reactiveValues(gene.top = list(),
                                 list.gene.analysis = list(),
                                 result.gene.info = data.frame())
  
  url_app = a("here", href = "https://github.com/biocsuwb/EnsembleFS")
  url_package = a("here", href = "https://github.com/biocsuwb/EnsembleFS-pacakge")
  email= a("here", href = 'http://matinf.uwb.edu.pl/pl/wydzial/kadra/pracownikii.php?cID=180')
  url.art.EnsembleFS = a('10.1007/978-3-030-50420-5_39', href = "https://link.springer.com/chapter/10.1007%2F978-3-030-50420-5_39")
  url.data.sample = a('https://www.cancer.gov/tcga', href = "https://www.cancer.gov/tcga")
  url.art.MDFS = a('10.1016/j.ins.2020.03.024', href = "https://www.sciencedirect.com/science/article/abs/pii/S0020025520302048?via%3Dihub")
  url.art.MCFS = a('10.1093/bioinformatics/btm486', href = 'https://academic.oup.com/bioinformatics/article/24/1/110/204931?login=true')
  url.art.MRMR =  a('10.1142/S0219720005001004' , href ="https://www.worldscientific.com/doi/abs/10.1142/S0219720005001004") 
  url.art.Lustgarten = a('PMCID: PMC2815476', href = "https://pubmed.ncbi.nlm.nih.gov/20351889/")
  url.art.RF = a('10.1023/A:1010933404324', href = "https://link.springer.com/article/10.1023/A:1010933404324")
  output$data.sample.set = renderUI({
    tagList("The dataset from The Cancer Genome Atlas TCGA", url.data.sample)
  })
  
  output$home.email <- renderUI({
    tagList(h4("Write to the help desk  ", email))
  })
  
  output$home.url_app <- renderUI({
    tagList(h4("The source code EnsembleFS is available here ", url_app))
  })
  
  output$home.url_package <- renderUI({
    tagList(h4("The R package EnsembleFS is available ", url_package))
  })
  
  output$url.article = renderUI({
    tagList(h4("DOI:", url.art))
  })
  
  
  output$url.MCFS = renderUI({
    tagList("DOI:", url.art.MCFS)
  })
  
  output$url.MDFS = renderUI({
    tagList("DOI:", url.art.MDFS)
  })
  
  output$url.MRMR = renderUI({
    tagList("DOI:", url.art.MRMR)
  })
  
  output$url.Lustgarten = renderUI({
    tagList("DOI:", url.art.Lustgarten)
  })
  
  output$url.RF = renderUI({
    tagList("DOI:", url.art.RF)
  })
  
  output$url.art.EnsembleFS = renderUI({
    tagList("DOI:", url.art.EnsembleFS)
  })
  
  output$home.url <- renderUI({
    tagList("The R package of EnsembleFS is available ", url_package)
  })
  #####################  
  output$info.load.main.data <- renderText({"Please select a data set"})
    

  data <- reactive({
    if(is.null(input$file1) && input$data.default == TRUE){
      read.csv2('www/LUAD_test_dataset_500samples.csv', check.names = FALSE)
    }else if(is.null(input$file1)){
        return(NULL)
    }else{ 
      data <- read.csv2(input$file1$datapath, check.names = FALSE)
      output$info.load.main.data <- renderText({""})
      return(data)  
      }
    })
    

   output$contents <- renderDataTable({
        data()[,1:9]
   })
 
    
    output$run <- renderUI({
        if(!is.null(data()) && length(input$methods) != 0 && length(input$cv)){
            actionButton("Run", "Run Feature Selection")
        }
    })
    
    observeEvent(input$Run, {
    withProgress(message = 'Feature Selection in progress. Please wait ...', {
    data <- data()
    nums <- unlist(lapply(data[,-input$num], is.numeric))
    if(input$num > ncol(data)){
        showModal(modalDialog(
            title = "Error in column decision",
            "The number of the decision column is greater than the number of columns!"
        ))
    }
    else if(FALSE %in% nums){
        showModal(modalDialog(
            title = "Error in data",
            "Columns in data must be of type numeric!"
        ))
    }
    
    else if(input$nvar > ncol(data[,-input$num]) && 'fs.mrmr' %in% input$methods){
        showModal(modalDialog(
            title = "Error in Number features for MRMR",
            "Error in Number features for MRMR methods, should not exceed the number of features"
        ))
    }
    
    else if(!all(data[,input$num] == 0 | data[,input$num] == 1)){
        showModal(modalDialog(
            title = "Error in column decision",
            "Decision must be binary!"
        ))
    }
    else if(all(data[,input$num] == 0) || all(data[,input$num] == 1)){
        showModal(modalDialog(
            title = "Error in column decision",
            "Both classes have to be represented!"
        ))
    }
    else{
    start_time <- Sys.time()

    res <-  ensembleFS(x = data[,-input$num],
                       y = data[,input$num],
                       methods = input$methods,
                       method.cv = input$cv,
                       params.cv = list(niter = input$niter, test.size = 0.3, k =3),
                       level.cor = input$level.cor,
                       params = list(adjust = input$adjust,
                                     feature.number = input$nvar,
                                     alpha = 0.05,
                                     use.cuda = TRUE,
                                     cutoff.method = c(input$cutoff)),
                       asm = c(input$methods),
                       model = c(input$methods)) 
 
  
    end_time <- Sys.time()
    print(end_time - start_time)
        
    result_full <-reactive({
        res
    }) 
        
    result_ranking <-reactive({
        res$ranking.feature
    }) 
    
    result_asm <- reactive({
        res$stability
    })
    
    ###add to gene infromation 
    store.result$gene.top <- res$selected.feature
  
    result_model <- reactive({
        res$model
    })
    
    info_app <- reactive({
        info <- list(methods = input$methods,
                     p.adjust = input$adjust,
                     level.corelation = input$level.cor,
                     validation = input$cv,
                     gene.info = input$type.info,
                     number.repeats = input$niter)
    })
    
    
    ## PLOT ASM
    plotly.reult.asm = reactive({
       data = result_asm()
         
       ggplot(data, aes(x= N, y = stability.asm, group= method, color= method)) +
       geom_line() +
       geom_point() +  
       scale_x_continuous(breaks= c(seq(0,100, by = 10))) + 
       labs(title= "ASM similarity measure between feature subsets vs top N variables.", y="ASM", x = "N") + 
       theme_light()+
       theme(legend.position = "bottom")
   })
    
    ## PLOT ACC
    plotly.reult.acc = reactive({
      data = result_model()
      
      ggplot(data, aes(x= N, y = mean.acc, group= method, color= method)) +
        geom_line() +
        geom_point() +  
        scale_x_continuous(breaks= c(seq(0,100, by = 10))) + 
        labs(title= "The accuracy vs top N variables.", y="ACC", x = "N") + 
        theme_light()+
        theme(legend.position = "bottom")
    })
    
    ## PLOT AUC
    plotly.reult.auc = reactive({
      data = result_model()
      
      ggplot(data, aes(x= N, y = mean.auc, group= method, color= method)) +
        geom_line() +
        geom_point() +  
        scale_x_continuous(breaks= c(seq(0,100, by = 10))) + 
        labs(title= "Area under the ROC curve vs top N variables.", y="AUC", x = "N") + 
        theme_light() +
        theme(legend.position = "bottom")
    })
    
    ## PLOT MCC
    plotly.reult.mcc = reactive({
      data = result_model()
      
      ggplot(data, aes(x= N, y = mean.mcc, group= method, color= method)) +
        geom_line() +
        geom_point() +  
        scale_x_continuous(breaks= c(seq(0,100, by = 10))) + 
        labs(title= "Matthews Correlation Coefficient vs top N variables.", y="MCC", x = "N") + 
        theme_light()+
        theme(legend.position = "bottom")
    })
    
    
    output$ranking <- renderDataTable(result_ranking())
    
    output$stability <- renderDataTable(result_asm()) 
    
    output$model <- renderDataTable(result_model())
  
    output$plot.stab <- renderPlotly(plotly.reult.asm())
  
    output$plot.model.acc <- renderPlotly(plotly.reult.acc())
    
    output$plot.model.auc <- renderPlotly(plotly.reult.auc())
    
    output$plot.model.mcc <- renderPlotly(plotly.reult.mcc())
    
    #render download button
    output$downloadAsm <- renderUI({
        if(!is.null(result_asm())){
            downloadButton('download_asm', 'Download CSV')
        }
    })
    
    output$downloadModel <- renderUI({
        if(!is.null(result_model())){
            downloadButton('download_model', 'Download CSV')
        }
    })
    
    output$downloadPlot <- renderUI({
        if(!is.null(plotly.reult.asm())){
            downloadButton('download_plot', 'Download PDF')
        }
    })
    
    ####
    output$downloadRanking <- renderUI({
        if(!is.null(result_ranking())){
            downloadButton('download_ranking', 'Download CSV')
        }
    })
    
     output$downloadZip <- renderUI({
         if(!is.null(result_full())){
             downloadButton('download_zip', 'Download Zip')
         }
     })
    
    ###
    #download handler
    output$download_asm <- downloadHandler(
        filename = function() {
            paste('stability', ".csv", sep = "")
        },
        content = function(file) {
            write.csv(result_asm(), file, row.names = FALSE)
        }
    )
    
    
    output$download_ranking <- downloadHandler(
        filename = function() {
            paste('ranking', ".csv", sep = "")
        },
        content = function(file) {
            write.csv(result_ranking(), file, row.names = FALSE)
        }
    )
    
    
    output$download_model <- downloadHandler(
        filename = function() {
            paste('model', ".csv", sep = "")
        },
        content = function(file) {
            write.csv(result_model(), file, row.names = FALSE)
        }
    )
    
    output$download_info <- downloadHandler(
        filename = function() {
            paste(result_info(), ".csv", sep = "")
        },
        content = function(file) {
            write.csv(result_info(), file, row.names = FALSE)
        }
    )
    
    
    output$download_plot = downloadHandler(
        filename = 'result.pdf',
        content = function(file) {
            pdf(file)
            
            arrangeGrob(print(plotly.reult.asm()),
                        print(plotly.reult.acc()), 
                        print(plotly.reult.auc()),
                        print(plotly.reult.mcc(), ncol = 4))  
            dev.off()
        })
    
    
    
    output$download_zip <- downloadHandler(
        filename = function() {
            paste("output", "zip", sep=".")
        },
        content = function(filename) {
            tmpdir <- tempdir()
            setwd(tempdir())
            print(tempdir())

            fs <- c("info.txt","ranking.csv", "stability.csv", "model.csv",
                    "utest.txt", "mcfs.txt", "mdfs1d.txt", "mdfs2d.txt", "mrmr.txt", "full_result.RData", "result.pdf")

            writeLines(deparse(info_app()), "info.txt")

            write.csv(result_ranking(), file = "ranking.csv")
            write.csv(result_asm(), file = "stability.csv")
            write.csv(result_model(), file = "model.csv")

            writeLines(deparse(fs.utest), "utest.txt")
            writeLines(deparse(fs.mcfs), "mcfs.txt")
            writeLines(deparse(fs.mdfs.1D), "mdfs1d.txt")
            writeLines(deparse(fs.mdfs.2D), "mdfs2d.txt")
            writeLines(deparse(fs.mrmr), "mrmr.txt")

            save(res, file = "full_result.RData")

            pdf('result.pdf')
            arrangeGrob(print(plotly.reult.asm()),
                        print(plotly.reult.acc()), 
                        print(plotly.reult.auc()),
                        print(plotly.reult.mcc(), ncol = 4))  
            dev.off()

            print (fs)
            zip(zipfile=filename, files=fs)
        },
        contentType = "application/zip"
    )
    }    
    })
    })
    ############GENE INFORMATION#################################################
    output$info.load.data <- renderText({
      if(length(store.result$gene.top) == 0){
        paste(" Load data from the FEATURE SELECTION  (relevant biomarkers)")
      }else{
        paste("Data from the FEATURE SELECTION (relevant biomarkers) was load")
      }
    })
    
    
    output$checkbox.condition.method <- renderUI({
      if(length(store.result$gene.top) != 0){
        radioButtons("condition.methods",
                     label = h4("Combination of a set of biomarkers"), 
                     choices = list("Union" = 'union',
                                    "Intersection" = 'intersect'),
                     selected = 1)
      }
    })
    
    output$get.information <- renderUI({
      if(length(store.result$gene.top) != 0 && length(input$condition.methods) != 0){
        actionButton("get.analysis", "GET ANALYSIS")
      }
    })
    
    ###################HELER FUNCTION###################
    funGProfiler = function(rel.var){
      df = list()
      for (i in 1:length(rel.var)){
        if(try(is.null(gost(query = rel.var[i]))))next
        all.info = gost(query = rel.var[i])
        df.res = all.info$result
        df[[i]] = data.frame(term = rel.var[i], source = df.res$source, term.ID = df.res$term_id, term.name = df.res$term_name)
      }
      df = do.call(rbind,df)
      return(df)
    }
    ###
    ranking.var <- function(list.imp.var, level.freq, n){
      result <- c()
      name.method <- c()
      for(name in names(list.imp.var)){
        var.one.method <- list.imp.var[[name]]
        var.imp <- list()
        for(i in var.one.method){
          if(length(i$name) < n) {var.imp <- append(var.imp, i$name)}
          else{var.imp <- append(var.imp, i$name[1:n])}
        }
        var.frequency <- as.data.frame(table(unlist(var.imp)))
        colnames(var.frequency) <- c('name', 'frequency')
        var.frequency <- var.frequency[order(var.frequency$frequency, decreasing = TRUE),]
        var.frequency <- var.frequency[var.frequency$frequency >= level.freq,]
        name.method <- append(name.method, substring(name, 4))
        result <- append(result, list(var.frequency$name))
      }
      names(result) <- name.method
      return(result)
    }
    #####################################################
    
    venn.for.methods <- eventReactive(input$geneNumber, {
      if(length(store.result$gene.top) != 0){
        if(input$cv == 'rsampling') level.freq = round(input$niter / 2)
        if(input$cv == 'kfoldcv') level.freq = round((3 * input$niter) / 2)
        var.venn.methods <- ranking.var(store.result$gene.top, level.freq, input$geneNumber)
        store.result$list.gene.analysis <- var.venn.methods
        venn(var.venn.methods, ilabels = TRUE, zcolor = "style", size = 25, cexil = 5, cexsn = 5, box = FALSE)
        #venn(venn00, ilabels = TRUE, zcolor = "style", plotsize = 15, ilcs = 1.2, sncs = 1.2)
      }
    })
    
    output$graph.venn.methods <- renderPlot({
      venn.for.methods()
    })
    
    observeEvent(input$get.analysis, {
      withProgress(message = 'Get analysis in progress. Please wait ...', {
        start_time <- Sys.time()
        
        gene.for.analysis <- store.result$list.gene.analysis
        if(input$condition.methods == 'intersect'){
          var.imp <- Reduce(intersect, gene.for.analysis)
        }
        else if(input$condition.methods == 'union'){
          var.imp <- unique(unlist(gene.for.analysis))}
        result <-  funGProfiler(var.imp)
        store.result$result.gene.info <- result
        name.source <- unique(result$source)
        list.var.source <- c()
        for(i in name.source){
          var <- result[result$source == i, 1]
          list.var.source <- append(list.var.source, list(var))
        }
        names(list.var.source) <- name.source
        end_time <- Sys.time()
        print(end_time - start_time)
        output$graph.venn.result <- renderPlot({
          venn(list.var.source, ilabels = TRUE, zcolor = "style", size = 25, cexil = 5, cexsn = 5, box = FALSE)
        })
      })})
    
    
    select.information.GProfiler = eventReactive(input$typeBase, {
      if(nrow(store.result$result.gene.info) != 0){
        all.info = store.result$result.gene.info
        if (input$typeBase == 1) df = subset(all.info, source == 'GO:MF')
        if (input$typeBase == 2) df = subset(all.info, source == 'GO:CC')
        if (input$typeBase == 3) df = subset(all.info, source == 'GO:BP')
        if (input$typeBase == 4) df = subset(all.info, source == 'KEGG')
        if (input$typeBase == 5) df = subset(all.info, source == 'REAC')
        if (input$typeBase == 6) df = subset(all.info, source == 'WP')
        if (input$typeBase == 7) df = subset(all.info, source == 'TF')
        if (input$typeBase == 8) df = subset(all.info, source == 'MIRNA')
        if (input$typeBase == 9) df = subset(all.info, source == 'HPA')
        if (input$typeBase == 10) df = subset(all.info, source == 'CORUM')
        if (input$typeBase == 11) df = subset(all.info, source == 'HP')
        if (input$typeBase == 12) df = all.info}
      return(df)
      
    })
    
    output$information = 
      renderDataTable(
        select.information.GProfiler(),
        options = list(pageLength = 5))
    
    
    output$save.information <- renderUI({
      if(!is.null(select.information.GProfiler())){
        downloadButton('save.data', 'Save')
      }
    })
    
    output$save.data <- downloadHandler(
      filename = function() {
        paste('gene_inforamtion', ".csv", sep = "")
      },
      content = function(file) {
        write.csv(select.information.GProfiler(), file, row.names = FALSE)
      }
    )
    
    
}

shinyApp(ui = ui, server = server)