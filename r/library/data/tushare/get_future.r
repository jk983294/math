source("library/data/tushare/tushare_func.r")

api <- Tushare::pro_api(token = read_tushare_token())

get_fut_basic(api, "CFFEX")
get_fut_basic(api, "DCE")
get_fut_basic(api, "CZCE")
get_fut_basic(api, "SHFE")
get_fut_basic(api, "INE")

get_fut_daily(api, "I1807.DCE")
