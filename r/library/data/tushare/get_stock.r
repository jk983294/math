source("library/data/tushare/tushare_func.r")

api <- Tushare::pro_api(token = read_tushare_token())
bar <- Tushare::pro_bar(token = read_tushare_token())

get_stk_basic(api)
get_stk_index_daily(api, "399300.SZ")
get_stk_daily(bar, "000001.SZ")
