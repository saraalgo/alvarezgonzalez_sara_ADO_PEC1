#código a emplear durante la PEC1

setwd("C:/Users/Sara/Desktop/alvarezgonzalez_sara_ADO_PEC1")
library(affy)

archivos <- ReadAffy(compress = TRUE)

#1.
feno <- archivos@phenoData

nombres=c(rep("CA3_cont",7),rep("CA3_act",7),rep("CA1_cont",6),rep("CA1_act",6),rep("DG_cont",7),rep("DG_act",7))
nombres<-paste0(nombres,"_")
CA3_control<-paste0(nombres[1:7],1:7)
CA3_exp<-paste0(nombres[8:14],1:7)
CA1_control<-paste0(nombres[15:20],1:6)
CA1_exp<-paste0(nombres[21:26],1:6)
DG_control<-paste0(nombres[27:33],1:7)
DG_exp<-paste0(nombres[34:40],1:7)
nombres <- c(CA3_control,CA3_exp,CA1_control,CA1_exp,DG_control,DG_exp)

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


#2.
#boxplot
colores = c(rep("red",14),rep("blue",12),rep("yellow",14))
boxplot(archivos,las=2, cex.axis=0.7, names=nombres, col=colores)
#histograma
hist(archivos,col=colores)
legend("topright", c("CA3", "CA1", "DG"), fill=c("red","blue","yellow"))
#clúster jerárquico
clust.euclid.average <- hclust(dist(t(exprs(archivos))),method="average")
plclust(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)

calidad <- qc(archivos)
plot(calidad)


#3.
normalizacion <- rma(archivos)
normalizacion
class(normalizacion)
matriz <- exprs(normalizacion)
head(matriz)

dim(matriz)

colnames(matriz) <- nombres
head(matriz)


#4.
#boxplot
par(mfrow=c(1,2))
boxplot(normalizacion,las=2, cex.axis=0.7,names=nombres, col=colores,ylab="Luminiscencia",main="Datos normalizados")
boxplot(archivos,las=2, cex.axis=0.7, names=nombres, col=colores,ylab="Luminiscencia",main="Datos crudos")
#histograma
par(mfrow=c(1,2))
hist(normalizacion,col=colores,main="Datos normalizados")
legend("topright", c("CA3", "CA1", "DG"), fill=c("red","blue","yellow"))
hist(archivos,col=colores,main="Datos crudos")
legend("topright", c("CA3", "CA1", "DG"), fill=c("red","blue","yellow"))