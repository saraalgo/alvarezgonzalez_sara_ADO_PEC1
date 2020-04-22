#código a emplear durante la PEC1

setwd("C:/Users/Sara/Desktop/alvarezgonzalez_sara_ADO_PEC1")
library(affy)

archivos <- ReadAffy(compress = TRUE)

#1.
feno <- archivos@phenoData
nombres=c(rep("Control CA3",7),rep("Learning Activated CA3",7),rep("Control CA1",6),rep("Learning Activated CA1",6),rep("Control DG",7),rep("Learning Activated DG",7))
feno@data[,1]=nombres
feno@data
feno$sample

fenotipo<-file.path("matrix.txt")
fenotipo <- as.matrix(read.table(file = fenotipo, header = T, sep = "\t", row.names = 1, as.is = T))
dim(fenotipo)
#las filas serán los genes y sus expresiones
print(rownames(fenotipo), max = 10)
#las columnas serán los anteriores nombres de los archivos, que de hecho podemos cambiar como hicimos en la otra manera de hacerlo anterior
colnames(fenotipo)
sapply(fenotipo, FUN = class)
str(fenotipo)
colnames(fenotipo)=nombres
colnames(fenotipo)


