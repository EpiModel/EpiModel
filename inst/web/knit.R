
## knit/purl EpiModel web tutorials

library(knitr)

file.remove("WebTutorials/BasicDCMs.R")
purl("WebTutorials/BasicDCMs.Rmd", output = "WebTutorials/BasicDCMs.R")

file.remove("WebTutorials/BasicICMs.R")
purl("WebTutorials/BasicICMs.Rmd", output = "WebTutorials/BasicICMs.R")

file.remove("WebTutorials/BasicNet.R")
purl("WebTutorials/BasicNet.Rmd", output = "WebTutorials/BasicNet.R")

file.remove("WebTutorials/NewDCMs.R")
purl("WebTutorials/NewDCMs.Rmd", output = "WebTutorials/NewDCMs.R")

file.remove("WebTutorials/NewNet.R")
purl("WebTutorials/NewNet.Rmd", output = "WebTutorials/NewNet.R")
