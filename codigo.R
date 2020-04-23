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


#6 y 7.
feno@data[,2] <- c(rep("CA3",14),rep("CA1",12),rep("DG",14))
feno@data[,3] <- c(rep("Control",7),rep("Activado",7),rep("Control",6),rep("Activado",6),rep("Control",7),rep("Activado",7))
colnames(feno@data)[2] <- "area"
colnames(feno@data)[3] <- "estado"
head(feno@data,2)
area <- factor(feno$area)
estado <- factor(feno$estado)

grupo <- interaction(area,estado)
plotMDS(matriz,col=as.numeric(grupo))
modelo <- model.matrix(~0+grupo)
ajuste <- voom(matriz,modelo,plot=T)

ajuste2 <- lmFit(ajuste,modelo)
head(coef(ajuste2))

contraste_CA3 <- makeContrasts(grupoCA3.Control - grupoCA3.Activado, levels = colnames(coef(ajuste2)))
contraste_CA3
gen_CA3 <- contrasts.fit(ajuste2,contraste_CA3)
gen_CA3 <- eBayes(gen_CA3)
table_CA3 <- topTable(gen_CA3,sort.by = "P",n=Inf)
head(table_CA3,5)
length(which(table_CA3$adj.P.Val < 0.05))
table_CA3$Genes <- rownames(table_CA3)
table_CA3 <- table_CA3[,c("Genes",names(table_CA3)[1:6])]
write.table(table_CA3,file="CA3_Control_vs_Activado.txt",row.names = F,sep = "\t",quote = F)

contraste_CA1 <- makeContrasts(grupoCA1.Control - grupoCA1.Activado, levels = colnames(coef(ajuste2)))
contraste_CA1
gen_CA1 <- contrasts.fit(ajuste2,contraste_CA1)
gen_CA1 <- eBayes(gen_CA1)
table_CA1 <- topTable(gen_CA1,sort.by = "P",n=Inf)
head(table_CA1,5)
length(which(table_CA1$adj.P.Val < 0.5))
table_CA1$Genes <- rownames(table_CA1)
table_CA1 <- table_CA1[,c("Genes",names(table_CA1)[1:6])]
write.table(table_CA1,file="CA1_Control_vs_Activado.txt",row.names = F,sep = "\t",quote = F)

contraste_DG <- makeContrasts(grupoDG.Control - grupoDG.Activado, levels = colnames(coef(ajuste2)))
contraste_DG
gen_DG <- contrasts.fit(ajuste2,contraste_DG)
gen_DG <- eBayes(gen_DG)
table_DG <- topTable(gen_DG,sort.by = "P",n=Inf)
head(table_DG,5)
length(which(table_DG$adj.P.Val < 0.05))
table_DG$Genes <- rownames(table_DG)
table_DG <- table_DG[,c("Genes",names(table_DG)[1:6])]
write.table(table_DG,file="DG_Control_vs_Activado.txt",row.names = F,sep = "\t",quote = F)

#8.

length(which(table_CA3$adj.P.Val < 0.5 & table_CA1$adj.P.Val < 0.5))
length(which(table_CA1$adj.P.Val < 0.5 & table_DG$adj.P.Val < 0.5))
length(which(table_CA3$adj.P.Val < 0.05 & table_DG$adj.P.Val < 0.05))
length(which(table_CA3$adj.P.Val < 0.5 & table_CA1$adj.P.Val < 0.5 & table_DG$adj.P.Val < 0.5))