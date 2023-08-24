(my_dates <- as.Date(c("2025-10-05", "2022-01-17", "2023-06-28", "2024-12-22")))
(out_dt <- data.table::fifelse(my_dates < as.Date("2024-01-01"), my_dates + 10, my_dates))
