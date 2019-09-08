df <- data.frame(cbind(matrix(rnorm(30, 1), ncol = 5), c(1, seq(5))))
write.csv(df, "/tmp/r_output.csv")
write.csv(df, "/tmp/no_row_index_output.csv", row.names = FALSE)


write.table(df, "/tmp/r_table.csv", sep = ",", quote = TRUE, na = "NA")
