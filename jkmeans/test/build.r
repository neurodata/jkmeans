library(devtools)

setwd("~/git/jkmeans")
build('jkmeans')
install.packages("jkmeans_1.0.tar.gz", repos = NULL, type = "source")

require("jkmeans")



