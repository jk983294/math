library("Tushare")

read_tushare_token <- function(path_ = "~/.ssh/tushare.token") {
    readLines(path_)
}

get_fut_basic <- function(api_, exch) {
    fut_basic <- api_(api_name = "fut_basic", exchange = exch, fut_type = 1)
    file_name <- paste("/tmp/fut_basic", exch, "csv", sep = ".")
    write.csv(fut_basic, file_name)
    message("get_fut_basic to file: ", file_name)
    
    fut_basic_cont <- api_(api_name = "fut_basic", exchange = exch, fut_type = 2)
    file_name <- paste("/tmp/fut_basic_cont", exch, "csv", sep = ".")
    write.csv(fut_basic_cont, file_name)
    message("get_fut_basic_cont to file: ", file_name)
}

get_fut_daily <- function(api_, symbol = "I1807.DCE", start_date_ = "20180101") {
    fut_daily <- api_(api_name = "fut_daily", ts_code = symbol, start_date = start_date_)
    file_name <- paste("/tmp/fut", symbol, "csv", sep = ".")
    write.csv(fut_basic, file_name)
    message("get_fut_daily to file: ", file_name)
}

get_stk_basic <- function(api_) {
    stock_basic <- api(api_name = "stock_basic")
    file_name <- paste("/tmp/stk_basic", "csv", sep = ".")
    write.csv(stock_basic, file_name)
    message("get_stk_basic to file: ", file_name)
}

get_stk_index_daily <- function(api_, symbol = "399300.SZ", start_date_ = "20180101") {
    index_data = api_(api_name = "index_daily", ts_code = symbol, start_date = start_date_)
    file_name <- paste("/tmp/index", symbol, "csv", sep = ".")
    write.csv(index_data, file_name)
    message("get_stk_index_daily to file: ", file_name)
}

get_stk_daily <- function(bar_, symbol = "000001.SZ", start_date_ = "20180101", fq = "hfq") {
    stk_data_adj <- bar_(ts_code = symbol, start_date = start_date_, adj = fq)
    file_name <- paste("/tmp/stk", symbol, fq, "csv", sep = ".")
    write.csv(stk_data_adj, file_name)
    message("get_stk_daily to file: ", file_name)
}
