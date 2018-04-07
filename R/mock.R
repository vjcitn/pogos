#' obtain an example trace set stored locally, for irinotecan and selected cell lines
#' @return an instance of DRTraceSet
#' @examples
#' iri = iriCCLE()
#' iri
#' plot(iri)
#' @export
iriCCLE = function() .load_mock("iriCCLE")

.load_mock = function (stub) 
{
    get(load(system.file(sprintf("mocks/%s.rda", stub), package = "pogos")))
}

