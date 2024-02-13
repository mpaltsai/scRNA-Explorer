#### scRNA-Seq: An end-user online tool for single cell RNA-seq data analysis featuring gene correlation and data filtering

##### Highlights
* Single cell RNA-Seq analysis utilises feature selection, to classify and characterise cell clusters withing a heterogeneous cell population. To our knowledge, there is no public pipeline to interrogate this plethora of data to identify by co-expression correlation the putative unknown or additional function of a specific gene.
* Utilizing i. Filtering out uninformative cells in an interactive manner via a web interface and ii. Gene correlation analysis coupled with an extra step of evaluating the biological importance of these correlations, we developed a pipeline to address the above question.
* scRNA-Explorer pipeline allows users to interrogate in an interactive manner scRNA-sequencing data sets to explore via gene expression correlations possible function(s) of a gene of interest

##### Workflow components
1. File input(s) upload, preprocessing and initial filtering of scRNA-Seq data
2. i. Visualisation of various metrics characterising the depth and effectiveness of sequencing through plots  
   ii. Filtering out non-informative cells with user defined cutoffs
3. Cell type clustering accompanied with plots related to functional annotation of cell clusters and expression of marker genes. Available for both filtered and unfiltered (though slightly filtered) data
4. Gene correlation  of a user defined gene against all genes present:
   i. Across all cells in the assay, with both filtered and unfiltered data
   ii. Within cell type clusters
5. Gene enrichment to functionally characterise correlated or anti-correlated genes to the selected gene

Please consider the following:
* This pipeline is suitable for a step-by-step usage as it depends on reactive components which update along the way.
* Experimental designs and sequencing platforms are highly diverse and may present limitations during the analysis. Error messages will inform you for almost any failure. If this is not true for your case, please contact us at ibaltsavia@gmail.com iliopj@med.uoc.gr
* scRNA-Seq is available at http://bioinformatics.med.uoc.gr/shinyapps/app/scrnaexplorer and as Docker image from DockerHub with:  docker pull mpaltsai/scrnaexplorer
  
