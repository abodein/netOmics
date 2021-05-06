context("get_interaction_from_database")

db.data.frame <- as.data.frame(list(from = c("A", "B", "E", "A"), 
                                    to = c("B", "C", "F", "C")))
db.graph <- igraph::graph_from_data_frame(db.data.frame, directed = FALSE)

X <- c("A", "F")
X.list <- list(X)

#get_interaction_from_database(X, db = NULL, type = "db", user.ego = FALSE)

test_that("get_interaction_from_database fails on invalid input", {
    # if X fails, X recieve a default value
    # same for type and user.ego
    expect_error(get_interaction_from_database(X = NULL), "'db' must be an igraph or data.frame object")
    expect_error(get_interaction_from_database(X = NULL, db = ""), "'db' must be an igraph or data.frame object")
})

test_that("get_interaction_from_database works", {
    expect_message(get_interaction_from_database(X = NULL, db = db.data.frame), "X is NULL, returning an empty graph")
    expect_is(get_interaction_from_database(X = X, db = db.data.frame, type = "db", user.ego = FALSE), "interaction.igraph")
    expect_is(get_interaction_from_database(X = X, db = db.graph, type = "db", user.ego = FALSE), "interaction.igraph")
    
    expect_is(get_interaction_from_database(X = X.list, db = db.graph, type = "db", user.ego = FALSE), "list.igraph")
    
})
