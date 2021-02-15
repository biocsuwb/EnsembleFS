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
ui <- fluidPage(theme = shinytheme("sandstone"), 
    navbarPage("WebEFS: ensemble feature selection methods for analysis of molecular data", 
               tabPanel('Home',
                        mainPanel(
                            h2('Welcome to WebEFS-tools '),
                            hr(),
                            div(
                                h4("WebEFS is tool that allows the user to: "),
                                h5('- filter the most informative biomarkers from molecular data generated from high-throughput microarray experiments
                  that could be a new diagnostic/prognostic markers or therapeutic targets;'),
                                h5('- establish the selected parameters for predictive models, such as the best number of most informative variables/biomarkers;'),
                                h5('- examine the impact of correlation between informative features on the predictive power of predictive model;'),
                                h5('- evaluate stability of biomarker sets and performance of predictive models;'),
                                h5('- find information about gene collection (gene ontology, pathways, tissue specificity, miRNA targets,
                  regulatory motif, protein complexes, disease phenotypes) in several biological databases.'),
                                h5('It can be applied to two-class problems. WebEFS based on the several fil-ter  feature  selection  algorithms,  such  as  the  U-test,  the  Monte  Carlo  FeatureSelection  (MCFS),
                                the  MultiDimensional  Feature  Selection  (MDFS)  and  theMinimum Redundancy Maximum Relevance (MRMR) for discovering the mostimportant  biomarkers  and  used  the  machine  learning  algorithms  to  evaluatequality  of  the  set  of  variables.
                                Predictive models are built using the Random Forest algorithm.'),
                                h5('The information about each of the biomarkers was obtained from the biological databases, namely Gene Ontology (molecular function, cellular component, biological process),
                                Kyoto Encyclopedia of Genes and Genomes (pathways), Reactome (pathways), WikiPathways (pathways), Transfac (regulatory motif),
                                miRNA targets (miRTarBase), Human  Protein  Atlas (tissue specificity),  CORUM  protein  complexes, and Human Phenotype Ontology (human  disease phenotypes).'),
                            hr(),
                            h2('Overview'),
                            img(src="Procedure2.png", align = "center",height='700px',width='1100px'),
                            h5('For details on used notation, please refer to the Help -> Terminology'),
                            h5('For examples on how to use WebBFS-tools, please refer to the Help -> Tutorial'),
                            h5('For description of sample dataset, please refer to the Help -> Sample Data'),
                            hr(),
                            h4('Contact'),
                            uiOutput("home.email"),
                            uiOutput("home.url"),
                            hr(),
                            h5('Developed by Pavel Hrablis, Aneta Polewko-Klim'),
                            h5('Institute of Computer Science, University of Bialystok, Bialystok, Poland'),
                            h5('Computational Center, University of Bialystok, Bialystok, Poland ')
                        ),
               )),
               tabPanel('FEATURE SELECTION',
                        sidebarPanel(
                            fileInput('file1', 'Load file',
                                      accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                            
                            
                            numericInput("num",
                                         label = h4("Column number of decision variable"),
                                         value = 1, min = 1),
                            
                            checkboxGroupInput("methods",
                                               label = h4("Feature selection methods (FS)"), 
                                               choices = list("U-test" = 'fs.utest',
                                                              "MRMR" = 'fs.mrmr',
                                                              "MCFS-ID" = 'fs.mcfs',
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
                                                            "SGoF" = "SGoF",
                                                            "None" = "none"), 
                                             selected = 'fdr')
                            ),
                            
                            conditionalPanel(
                                condition = "input.methods.includes('fs.mrmr')",
                                numericInput("nvar", label = h4("Number of relevant variables"), value = 10),
                                helpText("Note: hyperparameter of MRMR,
                                         setup no more than the number of all variables")
                            ),
                            
                            sliderInput("level.cor", label = h4("Correlation coefficient"), min = 0, max = 1, value = 0.75),
                            
                            radioButtons("cv", label = h4("Validation methods"),
                                         choices = list("3-fold cross-validation" = 'kfoldcv', "random sample (test set 30%)" = 'rsampling'), 
                                         selected = 'kfoldcv', inline = TRUE),
                            
                            sliderInput("niter", label = h4("Number iteration"), min = 1, max = 30, value = 10),
                            
                            uiOutput("run"),

                        ),
                        mainPanel(
                            
                            tabsetPanel(
                                type = "tabs",
                                tabPanel("Data",  dataTableOutput('contents')),
                                tabPanel("Ranking List", dataTableOutput("ranking"), uiOutput('downloadRanking')),
                                tabPanel("FS Stability", tableOutput("stability"), uiOutput('downloadAsm')),
                                tabPanel("Model Accuracy", tableOutput("model"), uiOutput('downloadModel')),
                                tabPanel("Plots",
                                         plotOutput("plot.stab"),
                                         plotOutput("plot.model.acc"),
                                         plotOutput("plot.model.auc"),
                                         plotOutput("plot.model.mcc"),
                                         br(),
                                         uiOutput('downloadPlot')),
                                tabPanel("Download Zip", br(), uiOutput('downloadZip'))
                            ))),
               
               tabPanel('GENE INFORMATION',
                        sidebarPanel(
                            textOutput('info.load.data'),        
                            hr(), 
                            h4("Parameters"),
                            hr(),
                            sliderInput("geneNumber",
                                        label = "Number of top N biomarkers",
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
               tabPanel('HELP',
                        tabsetPanel(
                            type = "tabs",
                            tabPanel('Terminology',
                                     h2('Feature selection algorithm'),
                                     h4('U-test'),
                                     h5('The U-test is a nonparametric statistical test that assigns probability to the hypothesis that two samples corresponding'),
                                     h5('to two decision classes (e.g. normal and tumor tissue) are drawn from populations with the same average value.'),
                                     h4('MDFS'),
                                     h5('MDFS method measures decrease of the information entropy of the decision variable due to knowledge of k-dimensional tuples of variables'),
                                     h5('and measures influence of each variable in the tuple. It performs an exhaustive search over all possible k-tuples and assigns to each'),
                                     h5('variable a maximal information gain due to given variable that was achieved in any of the k-tuple that included this variable. '),
                                     h5('MDFS-1D is one-dimensional version of this algorithm. MDFS-2D is two-dimensional version of this algorithm.'),
                                     uiOutput("url.MDFS"),
                                     h4('MCFS-IG'),
                                     h5('Monte Carlo Feature Selection and Interdependency Discovery is a Monte Carlo method-based tool for feature selection.'),
                                     uiOutput("url.MCFS"),
                                     h4('MRMR'),
                                     h5('Minimum redundancy maximum relevance feature selection.'),
                                     uiOutput("url.MRMR"),
                                     h2('Classifier'),
                                     h4('Random Forest'),
                                     h5('Random forest is an ensemble of decision trees, where each tree is built on a different bagging sample of the original data set.'), 
                                     h5('For each split, a subset of variables is selected randomly and the one is selected that allows to achieve the highest Gini coefficient'),
                                     h5('for the resulting leaves. Random Forest works well on data sets with a small number of objects, has few tunable parameters that do not'),
                                     h5('relate directly to the data, and very rarely fails.'),
                                     uiOutput("url.RF"),
                                     h2('Measuring the quality of models'),
                                     h4('AUC'),
                                     h5('AUC - ROC curve is a performance measurement for the classification problems at various threshold settings.'),
                                     h5('ROC is a probability curve and AUC represents the degree or measure of separability.'),
                                     h4('ACC'),
                                     h5('The accuracy of a machine learning classification algorithm is one way to measure how often the algorithm classifies a data point correctly.'),
                                     h5('Accuracy is the number of correctly predicted data points out of all the data points.'),
                                     h4('MCC'),
                                     h5('Matthews correlation coefficient or phi coefficient is used in machine learning as a measure of the quality of binary classifications.'),
                                     h2('Measuring the stability of feature selection'),
                                     h5('The  total  stability  of  feature selection method  is  measured  as  the  average  of  the  pairwise  similarity  for  all  pairs'),
                                     h5('of the most informative feature subsets from n runs of a model.'),
                                     h4('ASM'),
                                     uiOutput("url.Lustgarten")
                            ),
                            tabPanel('Tutorial', 
                                     h2('YouTube: WebEFS tutorial'),
                                     HTML('<iframe width="560" height="315" src="https://www.youtube.com/embed/e5oHiaigA68" frameborder="0" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen"></iframe>'),
                                     hr(),
                                     h2('Example'),
                                     h4('Feature selection process'),
                                     h5('The process of the selection of the most informative biomarkers includes the following steps:'),
                                     h5('1. Navigate to the FEATURE SELECTION tab.'),
                                     h5('2. Load csv/txt file (separator = ";", decimal = ",") and enter the column number for the binary decision variable'),
                                     h5('3. Choose the following parameters:'),
                                     h5('- Feature selection methods: U-test, MRMR'),
                                     h5('- Multitest correction: fdr'),
                                     h5('- Number of relevant variables: 100'),
                                     h5('- Correlation coefficient: 0.75'),
                                     h5('- Validation methods: 3-fold cross-validation'),
                                     h5('- Number of repetitions: 10'),
                                     h5('4. Press RUN FEATURE SELECTION'),
                                     h5('5. Navigate to RANKING LIST for the set of most informative biomarkers'),
                                     h5('6. Navigate to FS STABILITY for stability calculation results'),
                                     h5('7. Navigate to MODEL ACCURACY for the for model building results'),
                                     h5('8. Navigate to PLOT for the to visualize stability results and build a model'),
                                     h5('9. Navigate to DOWNLOAD ZIP for the to download all results in one archive'),
                                     h4('Searching information about genes'),
                                     h5('The process of aggregating information about the most informative biomarkers includes the following steps:'),
                                     h5('1. Navigate to GENE INFORMATION for biological information on top-100 biomarkers'),
                                     h5('2. Choose the following parameters:'),
                                     h5('- Number genes: 100'),
                                     h5('- Biomarkers set: union'),
                                     h5('- Data bases: all'),
                                     h5('3. Press GET ANALYSIS')
                                    
                            ),
                            tabPanel('Sample Data', 
                                     h2('Study dataset'),
                                     hr(),
                                     h4('Experimental data'),
                                     h5('The RNA-sequencing data of tumor-adjacent normal tissues of lung adenocarcinoma cancer patients.'),
                                     uiOutput("data.sample.set"),
                                     h4('Preprocessing of data'),
                                     h5('The preprocesing of data involved standard steps for RNA-Seq data. First the log2 transformation was performed.'), 
                                     h5('Then features with zero and near zero (1%) variance across patients were removed.'), 
                                     h5('The number of probes was randomly limited to 2000 DEGs.'),
                                     h5('The tumor tissue samples were randomly limited to 100 samples'),
                                     h5('The number of normal tissue samples is equal to 59 samples')
                                     
                                     
                            ),
                            tabPanel('WebEFS Flow Chart',
                                     img(src="FullProtocole.png", align = "center",height='600px',width='1200px'),
                                     h5('Fig1. Pipeline of the procedure to select the potential diagnostic/prognostic molecular markers.'),
                                     
                                     hr(),
                                     h4('For more details on used procedure for building predictive models, please refer to: '),
                                     h5('A. Polewko-Klim, W.R. Rudnicki. Analysis of Ensemble Feature Selection for Correlated High-Dimensional RNA-Seq Cancer Data.'),
                                     h5('In: Krzhizhanovskaya V. et al. (eds) Comp.Scien. ICCS 2020. Lecture Notes in Computer Science 12139 (2020), 525-538, Springer, Cham.'),
                                     uiOutput("url.art.benchmark")
                            )
                            )),
               tabPanel('Licence Information',
                        h4('Availability and requirements:'),
                        hr(),
                        h5('Project name: WebEFS: ensemble feature selection methods for analysis of molecular data'),
                        hr(),
                        h5('Source code: https://github.com/biocsuwb/WebEFS'),
                        hr(),
                        h5('Operation system(s): Web based, Platform independent'),
                        hr(),
                        h5('Programming language: R, R SHINY'),
                        hr(),
                        h5('Other requirements: Modern Browser'),
                        hr(),
                        h5('License: BSD License'),
                        hr(),
                        h5('Any restrictions to use by non-academics: None.'),
                        hr(),
                        h5('It is available from GitHub (https://github.com/biocsuwb/WebEFS) and is free open source software under an MIT license.'))
    
))

server <- function(session, input, output) {
    
    store.result <- reactiveValues(gene.top = list(),
                                   list.gene.analysis = list(),
                                   result.gene.info = data.frame())
    
    
    url = a("here", href = "https://github.com/biocsuwb/ensembleFS")
    email= a("here", href = 'http://matinf.uwb.edu.pl/pl/wydzial/kadra/pracownikii.php?cID=180')
    url.art.benchmark = a('DOI: 10.1007/978-3-030-50420-5_39', href = "https://link.springer.com/chapter/10.1007%2F978-3-030-50420-5_39")
    url.data.sample = a('https://www.cancer.gov/tcga',href = "https://www.cancer.gov/tcga")
    url.art.MDFS = a('DOI: 10.1016/j.ins.2020.03.024',href = "https://www.sciencedirect.com/science/article/abs/pii/S0020025520302048?via%3Dihub")
    url.art.MCFS = a('DOI: 10.1093/bioinformatics/btm486', href = 'https://academic.oup.com/bioinformatics/article/24/1/110/204931?login=true')
    url.art.MRMR =  a('DOI: 10.1142/S0219720005001004' , href ="https://www.worldscientific.com/doi/abs/10.1142/S0219720005001004") 
    url.art.Lustgarten = a('PMCID: PMC2815476', href = "https://pubmed.ncbi.nlm.nih.gov/20351889/")
    url.art.RF = a('DOI: 10.1023/A:1010933404324', href = "https://link.springer.com/article/10.1023/A:1010933404324")
    output$data.sample.set = renderUI({
        tagList("The dataset from The Cancer Genome Atlas TCGA", url.data.sample)
    })
    
    output$home.email <- renderUI({
        tagList("Write to the help desk  ", email)
    })
    
    output$home.url <- renderUI({
        tagList("The R package of WebEFS-tools is available ", url)
    })
    
    output$url.article = renderUI({
        tagList("For more information, see", url.art)
    })
    
    
    output$url.MCFS = renderUI({
        tagList("For more information, see", url.art.MCFS)
    })
    
    output$url.MDFS = renderUI({
        tagList("For more information, see", url.art.MDFS)
    })
    
    output$url.MRMR = renderUI({
        tagList("For more information, see", url.art.MRMR)
    })
    
    output$url.Lustgarten = renderUI({
        tagList("Lustgarten stability measure. For more information, see", url.art.Lustgarten)
    })
    
    output$url.RF = renderUI({
        tagList("For more information, see", url.art.RF)
    })
    
  ####
    
    output$home.url <- renderUI({
        tagList("The R package of WebBGS-tools is available ", url)
    })
    
    data <- reactive({
        validate(
            need(input$file1 != "", "Please select a data set")
        )
        inFile <- input$file1
        if (is.null(inFile)){return(NULL)}
        data <-  read.csv2(inFile$datapath)
        return(data)
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
    res <-  ensembleFS(x = data[,-input$num],
                       y = data[,input$num],
                       methods = input$methods,
                       method.cv = input$cv,
                       params.cv = list(niter = input$niter, test.size = 0.3, k =3),
                       level.cor = input$level.cor,
                       params = list(adjust = input$adjust, feature.number = input$nvar, alpha = 0.05),
                       asm = c(input$methods),
                       model = c(input$methods)) 
 
  
    
        
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
    
    plotinput <- function() {
        list(asm = graph.result(res$stability, 'stability'),
             acc =  graph.result(res$model, 'acc'),
             auc = graph.result(res$model, 'auc'),
             mcc = graph.result(res$model, 'mcc'))
    }
    
    result_plot <- reactive({
        plotinput()
    })
    
    output$ranking <- renderDataTable(result_ranking())
    
    output$stability <- renderTable(result_asm())
    
    output$model <- renderTable(result_model())
    
    
    output$plot.stab <- renderPlot(result_plot()['asm'])
    output$plot.model.acc <- renderPlot(result_plot()['acc'])
    output$plot.model.auc <- renderPlot(result_plot()['auc'])
    output$plot.model.mcc <- renderPlot(result_plot()['mcc'])
    
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
        if(!is.null(result_plot())){
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
            
            arrangeGrob(print(result_plot()[['asm']]),
                        print(result_plot()[['acc']]), 
                        print(result_plot()[['auc']]),
                        print(result_plot()[['mcc']]), ncol = 4)  
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
            arrangeGrob(print(result_plot()[['asm']]),
                        print(result_plot()[['acc']]),
                        print(result_plot()[['auc']]),
                        print(result_plot()[['mcc']]), ncol = 4)
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
            paste("Data from the FEATURE SELECTION will be loaded here")
        }else{
            paste("Data from the FEATURE SELECTION was load")
        }
    })
    
    
    output$checkbox.condition.method <- renderUI({
        if(length(store.result$gene.top) != 0){
            radioButtons("condition.methods",
                               label = h4("Biomarkers set"), 
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
      }
    })
    
    output$graph.venn.methods <- renderPlot({
      venn.for.methods()
    })

    observeEvent(input$get.analysis, {
        withProgress(message = 'Get analysis in progress. Please wait ...', {
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

        output$information = renderDataTable(
        select.information.GProfiler(),
        options = list(pageLength = 5)
    )
        
        
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