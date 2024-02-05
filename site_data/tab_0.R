tab0 <- tabPanel("Overview",
                 
                 fluidRow(
                   column(6,
                          wellPanel(
                            #QC plots
                            #create metrics
                            wellPanel(strong("Welcome to scRNA-Explorer"), align = "center")
                          )       
                   ),
                   
                   column(3,
                          
                   ),
                   column(3,
                          #helpText("Post-filtering")
                   )
                 ),
                 HTML("<br/>"),
                 
                 fluidRow(
                   column(6,
                          helpText("Main concept"),
                          wellPanel(
                           "scRNA explorer pipeline allows users to interrogate in an 
                          interactive manner scRNA-sequencing data sets to explore via gene
                        expression correlations possible function(s) of a gene of interest." 
                          )       
                 )),
                 
                 fluidRow(
                   column(6,
                          helpText("Significance"),
                          wellPanel(
                            "In addition to identifying expression
correlations with the gene of interest (our “bait”) we evaluate the expressions of the
putatively correlated genes in the cell population that expresses the gene of interest versus
the cell population that does not express our “bait”. We thus investigate the possibility of an
orchestrated effect by the co-expressed genes. If the expression of a putatively correlated
gene is comparable in the two populations, we deem this putative correlation of limited
biological significance. In this way, we keep only gene correlations to our bait that may hold
biological insights of a given process." 
                          )       
                   )),
                 fluidRow(
                   column(6,
                          helpText("Analysis components"),
                          wellPanel(
                            "To move towards this direction we coupled gene expression correlation with a review of the quality of cells sequenced
                            and gene enrichment on correlated genes that emerge from correlation analysis." 
                          )       
                   )),
                 fluidRow(
                   column(6,
                          helpText("Starting point: Data input/preprocessing"),
                          wellPanel(
                            "At the very start someone has to upload the results of a scRNA-seq experiment (or use a pre-loaded test case provided in the tab called 'Data input/preprocessing'.
                            After reading corresponding data, metadata regarding features present in the dataset are added and a minimal filtering of cells takes place." 
                          )       
                   )),
                 fluidRow(
                   column(6,
                          helpText("Cell quality assessment: Quality Control Plots"),
                          wellPanel(
                            "Afterwards quality control of cells takes place in tab 'Quality Control Plots'. Within this tab various metrics are plotted, 
                            giving the opportunity to filter out cells by specifying thresholds for these metrics. Upon user request, a filtering occurs 
                            and plots are re-drawn. In plots visualizing numerical values, range of sliders is initialized with minimum and maximum values." 
                          )
                   )
                 ),
                 fluidRow(
                   column(6,
                          helpText("Cell type clustering: Explore the presence of different cell types"),
                          wellPanel(
                            "Prior to gene correlation someone can review different cell types present in the assay." 
                          )
                   )
                 ),
                 fluidRow(
                   column(6,
                          helpText("Find related genes: Gene correlation"),
                          wellPanel(
                            "Having filtered and unfiltered (i.e. initial filtering that took place in the pre-processing phase) data, as well as cell type clusters, someone can run a correlation analysis 
                            for a desired gene among all single cells or cell type clusters in the assay. Resulting correlated genes to the desired gene are shown in tables. We generate three 
                            tables that correspond to all correlated, positive correlated and negative correlated genes" 
                          )
                   )
                 ),
                 fluidRow(
                   column(6,
                          helpText("Explore biological background: Gene enrichment analysis"),
                          wellPanel(
                            "Finally, someone can explore enriched biological phenomena to which correlated genes are involved (and thus the desired gene might be involved) through gene enrichment 
                            on correlated genes, which takes place in 'Gene enrichment analysis' tab." 
                          )
                   )
                 ),
                 fluidRow(
                   column(6,
                        
                          wellPanel(
                            "Most data objects generated by this pipeline are available for download upon user request (raw and filtered data, tables)" 
                          )
                   )
                 ),
                 
)#end of tabPanel
