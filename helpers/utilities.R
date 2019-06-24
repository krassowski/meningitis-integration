remove_leading_X = function(x) {
    x <- gsub("X", "", x)
    x
}

space_to_dot = function(x) {
    x <- gsub(" ", ".", x)
    x
}

dot_to_space = function(x) {
    x <- gsub("\\.", " ", x)
    x
}
