#' initial version of compound browser over pharmacoDb cells
#' @import shiny
#' @import ontoProc
#' @note simple shiny app demonstrating coverage of pharmacoDb compounds by CHEBI
#' @return only used for side effect of running shiny app
#' @examples
#' if (!requireNamespace("shiny")) stop("install shiny to use compoundsByCell")
#' if (interactive()) print(compoundsByCell())
#' @export
compoundsByCell = function() {
data(cell_lines_v1, package="pogos")
data(compounds_v1, package="pogos")
data(datasets_v1, package="pogos")
message("acquiring CHEBI")
if (!exists("cc")) cc = getChebiOnto()
message("done")
cellLineVec = cell_lines_v1[,"name"]
cellLineCodes = cell_lines_v1[,"id"]
compoundVec = compounds_v1[,"name"]
compoundCodes = compounds_v1[,"id"]
datasetVec = datasets_v1[,"name"]
datasetCodes = datasets_v1[,"id"]
names(datasetCodes) = datasetVec
names(cellLineCodes) = cellLineVec
names(compoundCodes) = compoundVec

ui = fluidPage(
 sidebarLayout(
  sidebarPanel(
   helpText("This app enumerates compounds that have been tested against selected cell lines in selected experiments."),
   selectInput("cell", "cell line", choices=cellLineVec, selected="MCF7"),
   #selectInput("compound", "compound", choices=compoundVec, selected=compoundVec[1]),
   selectInput("dataset", "dataset", choices=datasetVec, selected="CCLE")
   ),
  mainPanel( 
      textOutput("cell"), textOutput("dataset"),
      dataTableOutput("nres")
   )
  )
 )

server = function(input, output) {
 output$cell = renderText( input$cell )
 output$compound = renderText( input$compound )
 output$dataset = renderText( input$dataset )
 output$nres = renderDataTable({
    xx = GET(
       sprintf("https://pharmacodb.pmgenomics.ca/api/v1/intersections/2/%d/%d?indent=true",
       cellLineCodes[input$cell], datasetCodes[input$dataset] ))
    ans = fromJSON(readBin(xx$content, what="character"))
    validate(need(length(ans)>0, "cell line not tested in selected experiment"))
    if (length(ans)>0) {
      alld = sapply(ans, "[[", "drug")
      dn = as.character(alld[2,])
      chebn = tolower(cc$name)
      checo = names(cc$name)
      chepar = vapply(cc$parents[cc$id], function(x)x[1], character(1))
      cheparn = as.character(cc$name[chepar])
      lk = match(tolower(dn), chebn, nomatch=NA)
      chema = checo[lk]
      cheparn = cheparn[lk]
uwrap = function(id)
   sprintf("<A href='http://www.ebi.ac.uk/chebi/chebiOntology.do?chebiId=%s' target='_blank'>%s</A>", id, id)
      wrid = a(href=uwrap(chema), chema)
      nn = data.frame(drug=dn, chebi=uwrap(chema), parent=cheparn)
      } 
    else nn = data.frame(n=as.character(length(ans)))
    nn
    }, escape=FALSE)
}

print(shinyApp(ui, server))
}
