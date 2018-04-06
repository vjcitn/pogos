#onto_plot(oto, st@ontoTags, fontsize=50)
#    r = sapply(x@traces, function(x) f(x@DRProfiles[[1]]@responses))
#    data.frame(Cell_line_primary_name = x@cell_lines, resp = r, 
#        drug = x@drug, dataset = x@dataset)

#' trace extractor
#' @param x instance of DRTraceSet
#' @importFrom methods is
#' @examples
#' iri = iriCCLE()
#' str(traces(iri)[[1]])
#' @export
traces = function(x) {
 stopifnot(is(x, "DRTraceSet"))
 slot(x, "traces")
}

# #' profiles extractor
# #' @param x instance of DRTraceSet
# #' @examples
# #' iri = iriCCLE()
# #' str(DRProfiles(iri)[[1]])
# #' @export
# DRProfiles = function(x) {
#  stopifnot(is(x, "DRTraceSet"))
#  slot(x, "DRProfiles")
# }

#' DRProfSet is a class for managing dose-response information about 
#' cell lines from a pharmacogenomics dataset
#' @import methods
#' @importFrom graphics lines
#' @rdname DRProfSet
#' @aliases DRProfSet-class DRProfSet
#' @examples
#' if (interactive()) trs = DRTraceSet() else trs = iriCCLE()
#' ps = traces(trs)[[1]]
#' ps
#' getDrugs(ps)
#' @export
setClass("DRProfile", representation(cell_line = "character", 
    drug = "character", drug_code = "numeric", cell_line_code = "numeric", 
    dataset = "character", doses = "numeric", 
    responses = "numeric"))
setClass("DRProfSet", representation(cell_line = "character", 
    dataset = "character", DRProfiles = "list"))
setMethod("show", "DRProfSet", function(object) {
    cat(sprintf("DRProfSet with %d dose-response profiles for cell line %s\n dataset: %s\n", 
        length(object@DRProfiles), object@cell_line, object@dataset))
})

#' getDrugs extracts drug list
#' @rdname DRProfSet
#' @param x instance of DRProfSet
#' @return getDrugs: character vector
#' @aliases getDrugs,DRProfSet-method
#' @export
setGeneric("getDrugs", function(x) standardGeneric("getDrugs"))
setMethod("getDrugs", "DRProfSet", 
    function(x) vapply(x@DRProfiles, function(x) x@drug, character(1)))

#' subscripting on DRProfSet extracts a profile for a single drug 
#' whose name constitutes the index
#' @param x instance of DRProfSet
#' @param i character(1) drug name
#' @param j not used
#' @param \dots not used
#' @param drop logical(1) not used
#' @return a DRProfSet instance restricted to experiments involving the selected drug
#' @export
setMethod("[", c("DRProfSet", "character"), 
    function(x, i, j, ..., drop = TRUE) {
      ind = match(i, getDrugs(x))
      if (length(ind) == 0) 
          stop(sprintf("index drug (%s) not found\n", i))
      if (length(ind) > 1) {
        message("multiple matches found, using first")  # FIXME
        ind = ind[1]
      }
      x@DRProfiles = x@DRProfiles[ind]
      x
    })


DRtraces = function(drp, ylab = "response", ...) {
    doses = lapply(drp@DRProfiles, function(x) x@doses)
    resps = lapply(drp@DRProfiles, function(x) x@responses)
    xlim = range(unlist(doses))
    ylim = range(unlist(resps))
    plot(doses[[1]], resps[[1]], xlim = c(0.95, 1.05) * xlim, ylim = c(0.95, 1.05) * 
        ylim, ylab = ylab,  
        main = paste(drp@dataset, drp@cell_line, sep = " "), xlab = "dose", 
        type = "l", ...)
    if (length(doses) > 1) 
        uu = lapply(2:length(doses), function(x) lines(doses[[x]], resps[[x]], col = x, 
            lty = x))
    invisible(NULL)
}

lkc = function(cell_line = "MCF7", dataset = "CCLE") {
    data(cell_lines_v1)
    data(datasets_v1)
    data(compounds_v1)
    cellLineVec = cell_lines_v1[, "name"]
    cellLineCodes = cell_lines_v1[, "id"]
    compoundVec = compounds_v1[, "name"]
    compoundCodes = compounds_v1[, "id"]
    datasetVec = datasets_v1[, "name"]
    datasetCodes = datasets_v1[, "id"]
    names(datasetCodes) = datasetVec
    names(cellLineCodes) = cellLineVec
    names(compoundCodes) = compoundVec
    xx = GET(sprintf("https://pharmacodb.pmgenomics.ca/api/v1/intersections/2/%d/%d?indent=true", 
        cellLineCodes[cell_line], datasetCodes[dataset]))
    ans = fromJSON(readBin(xx$content, what = "character"))
    ans = lapply(ans, function(comb) {
        new("DRProfile", cell_line = cell_line, dataset = dataset, 
            drug = comb$compound$name, 
            drug_code = comb$compound$id, 
            cell_line_code = comb$cell_line$id, 
            doses = vapply(comb$dose_responses, 
                function(x) x$dose, numeric(1)), 
            responses = vapply(comb$dose_responses, 
                function(x) x$response, numeric(1)))
    })
    new("DRProfSet", cell_line = cell_line, dataset = dataset, DRProfiles = ans)
}

#' DRTraceSet class manages dose-response information for a single cell line, multiple drugs
#' @rdname DRTraceSet
#' @export
setClass("DRTraceSet", representation(cell_lines = "character", drug = "character", 
    dataset = "character", traces = "list"))
setMethod("show", "DRTraceSet", function(object) {
    cat(sprintf("DRTraceSet for %d cell lines, drug %s, dataset %s\n", 
        length(object@cell_lines), object@drug, object@dataset))
})
.traceDF = function(drts) {
    ds = lapply(drts@traces, function(x) x@DRProfiles[[1]]@doses)
    rs = lapply(drts@traces, function(x) x@DRProfiles[[1]]@responses)
    cls = vapply(drts@traces, function(x) x@cell_line, character(1))
    ns = vapply(ds, length, numeric(1))
    cls = rep(cls, ns)
    data.frame(dose = unlist(ds), response = unlist(rs), 
               cell_line = cls, dataset = drts@dataset, 
               drug = drts@drug)
}

#' @rdname DRTraceSet
#' @import ggplot2
#' @param x for plot: instance of DRTraceSet
#' @param y for plot: not used
#' @param \dots not used
#' @export
setMethod("plot", c("DRTraceSet", "missing"), function(x, y, ...) {
    df = .traceDF(x)  # defines cell_line, dose, response
    ggplot(df, aes(x = log(dose), y = response, colour = cell_line)) + geom_line()
})

#' DRTraceSet constructor for multiple cell lines, single drug, single dataset
#' @param cell_lines character vector of cell line names, must be found 
#' in `cell_lines_v1` data of pogos package
#' @param drug character(1) drug name in `compounds_v1`
#' @param dataset character(1) dataset known to pharmacodb.pmgenomics.ca
#' @note Will query pharmacodb for relevant dose-response information
#' @return instance of DRTraceSet
#' @examples
#' DRTraceSet()
#' @export
DRTraceSet = function(cell_lines = c("SK-ES-1", "TC-71", "MHH-ES-1", "HCC-56", "SK-HEP-1"), 
    drug = "Irinotecan", dataset = "CCLE") {
    if (length(drug) != 1) 
        stop("currently handles only one drug at a time")
    st1 = lapply(cell_lines, function(x) try(lkc(x, dataset), silent = TRUE))
    errs = vapply(st1, function(x) inherits(x, "try-error") || length(x@DRProfiles) == 
        0, logical(1))
    if (any(errs)) {
        if (all(errs)) 
            stop("no data for requested cell lines")
        message("no data for some requested cell lines:")
        message(paste(" ", cell_lines[which(errs)], sep = " "))
        st1 = st1[-which(errs)]
    }
    kp = seq_len(length(cell_lines))
    if (any(errs)) 
        kp = kp[-which(errs)]
    new("DRTraceSet", cell_lines = cell_lines[kp], drug = drug, dataset = dataset, 
        traces = lapply(st1, function(x) x[drug]))
}


#' DRProfSet manages all data from a given cell line from
#' a pharmacogenomics source
#' @param cell_line character(1) cell line name, entries in cell_lines_v1
#' @param dataset character(1) resource name, entries in datasets_v1
#' @return instance of DRProfSet
#' @examples
#' if (interactive()) DRProfSet()
#' @export
DRProfSet = function(cell_line = "MCF7", dataset = "CCLE") {
  try(lkc(cell_line=cell_line, dataset=dataset))
}
