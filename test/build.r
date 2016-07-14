library(devtools)

setwd("~/git/")
build('jkmeans')
install.packages("jkmeans_1.0.tar.gz", repos = NULL, type = "source")

require("jkmeans")



