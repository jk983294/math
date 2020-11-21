library(data.table)

# 1天10个bar, 两个ii
index <- 1:10
ticks <- 929000 + (1000 * index)
ticks <- rep(ticks, 4)
cp <- c(1.01^(1:10), (1.01^10) * (1.02^(1:10)), 1.02^(1:20))
pre_cp <- c(rep(1, 10), rep(1.01^10, 10), rep(1, 10), rep(1.02^10, 10))
ii <- c(rep(1, 20), rep(2, 20))
capt <- c(rep(0.5, 20), rep(2, 20))
beta <- c(rep(0.1, 20), rep(0.2, 20))
dates <- rep(c(rep(1, 10), rep(2, 10)), 2)
UNIVERSE <- rep(1, 40)

dt <- data.table(ukey = ii, DataDate = dates, ticktime = ticks, close = cp, vwap = cp, 
    pre_close = pre_cp, capt = capt, beta = beta, UNIVERSE = UNIVERSE)

# fill=pre_close[[1L]] 表示以组内第一个值来填充
dt[, `:=`(ret, close/shift(close, fill = pre_close[[1L]]) - 1), by = .(ukey, DataDate)]  # 每个bar收益
# open to now
dt[, `:=`(ret_on, close/close[[1L]] - 1), by = .(ukey, DataDate)]
# pre_close to now
dt[, `:=`(ret_pcn, close/pre_close[[1L]] - 1), by = .(ukey, DataDate)]
# pre_close to open
dt[, `:=`(ret_pco, close[[1L]]/pre_close[[1L]] - 1), by = .(ukey, DataDate)]
# today close to next day open
dt[ticktime == 930000L, `:=`(ret_co, shift(ret_pco, 1L, type = "lead", fill = 0)), 
    by = ukey]
dt[, `:=`(ret_co, first(ret_co)), by = .(ukey, DataDate)]
# open to close
dt[, `:=`(ret_oc, last(ret_on)), by = .(ukey, DataDate)]
# now to close
dt[, `:=`(ret_nc, last(close)/close - 1), by = .(ukey, DataDate)]
# now to next day open
dt[, `:=`(ret_no, (1 + ret_nc) * (1 + ret_co) - 1), by = .(ukey, DataDate)]
# previous day now to previous day close
dt[, `:=`(ret_pnc, shift(ret_nc, fill = 0)), by = .(ukey, ticktime)]
# previous day now to today open
dt[, `:=`(ret_pno, (1 + ret_pnc) * (1 + ret_pco) - 1)]
# vret: vwap to vwap ret
dt[, `:=`(vret, vwap/shift(vwap, fill = pre_close[[1L]]) - 1), by = .(ukey, DataDate)]

ret_cols <- grep("ret", names(dt), value = TRUE)
mret_cols <- paste0("m", ret_cols)
umret_cols <- paste0("u", mret_cols)

x <- 1:2
cap <- 2:3
weighted.mean(x, cap, na.rm = TRUE)  # 1.6 = 1 * (2/5) + 2 * (3/5)

dt[, `:=`((mret_cols), lapply(.SD, function(x) x - beta * weighted.mean(x, capt, 
    na.rm = TRUE))), by = .(DataDate, ticktime), .SDcols = ret_cols]
dt[, `:=`((umret_cols), lapply(.SD, function(x) x - weighted.mean(x, capt, na.rm = TRUE))), 
    by = .(DataDate, ticktime, UNIVERSE), .SDcols = mret_cols]
