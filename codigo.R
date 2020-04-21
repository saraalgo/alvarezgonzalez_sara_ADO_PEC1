#código a emplear durante la PEC1

setwd("C:/Users/Sara/Desktop/alvarezgonzalez_sara_ADO_PEC1")
library(affy)

#1.
archivos <- ReadAffy()
feno <- archivos@phenoData