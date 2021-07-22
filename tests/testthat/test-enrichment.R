context("ORA")

test_that("get_interaction_from_ORA fails on invalid imput", {
    query <- c("IL15", "CDHR5", "TGFA", 'C4B')
    expect_error(get_interaction_from_ORA(query, sources = "qowiudjh"))
    expect_error(get_interaction_from_ORA(query, sources = c("GO", "KEGG"))) 
})

test_that("get_interaction_from_ORA works ", {
     query <- c("IL15", "CDHR5", "TGFA", 'C4B')
     expect_is(get_interaction_from_ORA(query, sources = "GO"), "igraph")
     expect_is(get_interaction_from_ORA(query, sources = c("KEGG")), "igraph")
     
     
     query <- list("All" = c("IL15", "CDHR5", "TGFA", 'C4B'),
                   "c1" = c("IL15", "CDHR5", "TGFA"))
     expect_is(get_interaction_from_ORA(query, sources = "GO"), "list.igraph")
     query <- list(c("IL15", "CDHR5", "TGFA", 'C4B'), c("IL15", "CDHR5", "TGFA"))
     expect_is(get_interaction_from_ORA(query, sources = "GO"), "list.igraph")
     
})

test_that("get_ORA works", {
    query <- list("All" = c("IL15", "CDHR5", "TGFA", 'C4B'),
                  "c1" = c("IL15", "CDHR5", "TGFA"))
    expect_is(get_ORA(query), "data.frame")
    query <- c("IL15")
    expect_is(get_ORA(query), "data.frame")
})

test_that("get_go_info works", {
    expect_is(get_go_info("GO:0044216"), "data.frame")
})
