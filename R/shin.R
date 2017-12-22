#' initial version of compound browser over pharmacoDb cells
#' @note simple shiny app demonstrating coverage of pharmacoDb compounds by CHEBI
#' @export
compoundsByCell = function() {
library(pogos)
library(shiny)
data(cell_lines_v1)
data(compounds_v1)
data(datasets_v1)
library(ontoProc)
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
    if (length(ans)>0) {
      alld = sapply(ans, "[[", "drug")
      dn = as.character(alld[2,])
      chebn = tolower(cc$name)
      checo = names(cc$name)
      lk = match(tolower(dn), chebn, nomatch=NA)
      chema = checo[lk]
uwrap = function(id)
   sprintf("<A href='http://www.ebi.ac.uk/chebi/chebiOntology.do?chebiId=%s'>%s</A>", id, id)
      wrid = a(href=uwrap(chema), chema)
      nn = data.frame(drug=dn, chebi=chema, chbn=chebn[lk], url=uwrap(chema))
      } 
    else nn = data.frame(n=as.character(length(ans)))
    nn
    }, escape=FALSE)
}

print(shinyApp(ui, server))
}
