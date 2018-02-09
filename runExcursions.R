library(R.matlab)
library(excursions)

outputPath = commandArgs(TRUE)[1]
method = commandArgs(TRUE)[2]
subjStr = commandArgs(TRUE)[3]
u = as.numeric(commandArgs(TRUE)[4])
alpha = 1 - as.numeric(commandArgs(TRUE)[5])

mat = readMat(paste(c(outputPath,'sub',subjStr,'/',method,'/excursionsInput',method,'.mat'),collapse=''))

X = mat$wVeck
MCMCexc = excursions.mc(X,alpha=alpha,type='>',u=u)

writeMat(paste(c(outputPath,'sub',subjStr,'/',method,'/excursionsResults',method,'.mat'),collapse=''),excuFuncMCMC=MCMCexc$F,u=u)

