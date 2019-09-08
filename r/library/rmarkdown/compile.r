# knit direct on *.Rmd will calculated in a seprate session, while manually call render, you can use current session's data
library(rmarkdown)
render("library/rmarkdown/demo", "pdf_document")
