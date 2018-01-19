
library(pogos)

context("testTargets")
test_that("test suite knows how many functions are exported", {
    exported = ls("package:pogos")
    expect_true(length(exported) == 7)
})

context("testInterface")
test_that("query resolves", {
    qout = rxdbQuery_v1("cell_lines")  # yields 30; append '?all=true' to retrieve all
    expect_true(all(names(qout[[1]]) == c("id", "name")))
    expect_true(length(sapply(qout, function(x) x[[2]])) == 30)
    expect_error(DRTraceSet("SNU-719"))
    expect_true(is(DRTraceSet(c("MCF7", "SNU-719")), "DRTraceSet"))
})


