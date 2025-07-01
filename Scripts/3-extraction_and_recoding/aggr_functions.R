get_mean <- function(recoded, colname_prefix) {
    return(apply(recoded[,startsWith(colnames(recoded), colname_prefix)], 1, 
        function(x) if_else(
            sum(!is.na(x)) == 0, 
            NA, 
            mean(x, rm.na = TRUE)
        )
    ))
}

get_valid_flg <- function(recoded, colname_prefix, valid_values) {
    return(apply(
        recoded[,startsWith(colnames(recoded), colname_prefix)], 1,
        function(x) {
            for (val in as.array(x)) {
                if (val %in% valid_values)
                    return(TRUE)
            }
            return(FALSE)
        }
    ))
}

drop_col <- function(keep_col_mask, colname_prefix) {
    return(keep_col_mask & !startsWith(names(keep_col_mask), colname_prefix))
}