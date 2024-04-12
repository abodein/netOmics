context("get_interaction_from_database")

db.data.frame <- as.data.frame(list(from = c("A", "B", "E", "A"), 
                                    to = c("B", "C", "F", "C")))
db.graph <- igraph::graph_from_data_frame(db.data.frame, directed = FALSE)

X <- c("A", "F")
X.list <- list(X)
names(X.list) <- "2"

#get_interaction_from_database(X, db = NULL, type = "db", user.ego = FALSE)

test_that("get_interaction_from_database fails on invalid input", {
    # if X fails, X recieve a default value
    # same for type and user.ego
    expect_error(get_interaction_from_database(X = NULL), "'db' must be an igraph or data.frame object")
    expect_error(get_interaction_from_database(X = NULL, db = ""), "'db' must be an igraph or data.frame object")
    expect_error(get_interaction_from_database(X = NULL, db = data.frame("a" = c(1,2,3), 'b'= c(2,3,4))), "'db' must contains the columns 'from' and 'to'")
    
})

test_that("get_interaction_from_database works", {
    expect_message(get_interaction_from_database(X = NULL, db = db.data.frame), "X is NULL, returning an empty graph")
    expect_is(get_interaction_from_database(X = X, db = db.data.frame, type = "db", user.ego = FALSE), "interaction.igraph")
    expect_is(get_interaction_from_database(X = X, db = db.graph, type = "db", user.ego = FALSE), "interaction.igraph")
    
    expect_is(get_interaction_from_database(X = X.list, db = db.graph, type = "db", user.ego = FALSE), "list.igraph")
    expect_is(get_interaction_from_database(X = X.list, db = db.data.frame, type = "db", user.ego = FALSE), "list.igraph")
    
    expect_message(get_interaction_from_database(X = list("A" = "Z"), db = db.graph, type = "db", user.ego = FALSE), "no shared elements between X and db, return empty graph", fixed = TRUE)
    expect_is(get_interaction_from_database(X = X.list, db= db.graph, type = "db", user.ego = TRUE), "list.igraph")
    
    expect_message(get_interaction_from_database(X = list("A" = "Z"), db = db.data.frame, type = "db", user.ego = FALSE), "no shared elements between X and db, return empty graph", fixed = TRUE)
    expect_is(get_interaction_from_database(X = X.list, db = db.data.frame, type = "db", user.ego = TRUE), "list.igraph")
    
    expect_message(get_interaction_from_database(X = list("Z"), db = db.data.frame, type = "db", user.ego = FALSE), "no shared elements between X and db, return empty graph", fixed = TRUE)
    
})
