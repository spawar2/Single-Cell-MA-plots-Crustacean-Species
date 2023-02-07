# MA Plots for different species, Shrikant Pawar, revised 03/29/2021

# Read in data file
setwd("/Users/yalegenomecenter/Desktop/Mihika Single Cell Data/Mihika_Plos_2/Deseq-Volcano-Plosone")
data <- read.csv("PclaLF_evg_RSEM vs PclaBr_evg_RSEM.csv", header=TRUE)
y <- data[,c(3, 4)]

#####
plot(rowMeans(log2(abs(y))), log2(abs(y[, 1])-y[, 2]), cex=0.1, col = "grey", xlab="Mean of Normalized Counts (Log2)", ylab="Fold Change (Log2)", xlim=c(-15,15),
ylim=c(-20,20))
#####

plot(y$baseMean, y$log2FoldChange, xlim=c(0,50000), cex=0.1, col = "grey", xlab="Mean of Normalized Counts", ylab="Fold Change (Log2)")
abline(h=1.5, col="red")
abline(h=-1.5, col="green")

# Read in annotation id file
tag <- read.delim("/Users/yalegenomecenter/Desktop/Re__Volcano_plots_for_GPCR_manuscript/Pcla_LGRs.txt", header=FALSE)

df <- data.frame(matrix(NA,20,8))

for(i in 1:nrow(tag)) {
    row <- tag[i,]
    print(row)
    row1 <- data[which(data$countData...1.== row),]
    print(row1)
    df[i,] <- row1
}

x <- df[,c(3,4)]
x <- x[complete.cases(x),]


points(rowMeans(log2(abs(x))), log2(abs(x[, 1]))-abs(x[, 2]), col="brown", cex=0.3)

--------

tag <- read.delim("/Users/yalegenomecenter/Desktop/Re__Volcano_plots_for_GPCR_manuscript/Pcla_Neuropeptide.txt", header=FALSE)

df <- data.frame(matrix(NA,20,8))

for(i in 1:nrow(tag)) {
    row <- tag[i,]
    print(row)
    row1 <- data[which(data$countData...1.== row),]
    print(row1)
    df[i,] <- row1
}

x <- df[,c(3,4)]
x <- x[complete.cases(x),]

points(rowMeans(log2(abs(x))), log2(abs(x[, 1]))-abs(x[, 2]), col="dark green", cex=0.3)

---------

tag <- read.delim("/Users/yalegenomecenter/Desktop/Re__Volcano_plots_for_GPCR_manuscript/Pcla_Opsin.txt", header=FALSE)

df <- data.frame(matrix(NA,20,8))

for(i in 1:nrow(tag)) {
    row <- tag[i,]
    print(row)
    row1 <- data[which(data$countData...1.== row),]
    print(row1)
    df[i,] <- row1
}

x <- df[,c(3,4)]
x <- x[complete.cases(x),]

points(rowMeans(log2(abs(x))), log2(abs(x[, 1]))-abs(x[, 2]), col="orange", cex=0.3)


---------

tag <- read.delim("/Users/yalegenomecenter/Desktop/Re__Volcano_plots_for_GPCR_manuscript/Pcla_Orphans.txt", header=FALSE)

df <- data.frame(matrix(NA,20,8))

for(i in 1:nrow(tag)) {
    row <- tag[i,]
    print(row)
    row1 <- data[which(data$countData...1.== row),]
    print(row1)
    df[i,] <- row1
}

x <- df[,c(3,4)]
x <- x[complete.cases(x),]

points(rowMeans(log2(abs(x))), log2(abs(x[, 1]))-abs(x[, 2]), col="purple", cex=0.3)

---------

tag <- read.delim("/Users/yalegenomecenter/Desktop/Re__Volcano_plots_for_GPCR_manuscript/Pcla_SmallMolecule.txt", header=FALSE)

df <- data.frame(matrix(NA,20,8))

for(i in 1:nrow(tag)) {
    row <- tag[i,]
    print(row)
    row1 <- data[which(data$countData...1.== row),]
    print(row1)
    df[i,] <- row1
}

x <- df[,c(3,4)]
x <- x[complete.cases(x),]

points(rowMeans(log2(abs(x))), log2(abs(x[, 1]))-abs(x[, 2]), col="dark blue", cex=0.3)

-------------

legend("topright", inset=.02, title="Subclasses", c("Opsins", "Small-molecule receptors", "Neuropeptide Receptors", "LGRs", "Orphans"), col=c("orange", "dark blue", "dark green", "brown", "purple"), cex=0.8, lty=1)



