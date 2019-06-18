#load in rlog data from Aaron carlin


my_dat <- data.frame(read.delim("hDC_AR_rlog.txt"))
my_dat2 <- data.frame(my_dat[,9:17], row.names = my_dat[,1])

colnames(my_dat2) <- c( "hDC7_mock", "hDC7_ZIKVneg", "hDC7_ZIKVpos",
                         "hDC8_mock", "hDC8_ZIKVneg", "hDC8_ZIKVpos",
                         "hDC10_mock", "hDC10_ZIKVneg", "hDC10_ZIKVpos")

sample_dat <- data.frame(my_dat2[,2:3] - my_dat2$hDC7_mock,
                         my_dat2[,5:6] - my_dat2$hDC8_mock,
                         my_dat2[,8:9] - my_dat2$hDC10_mock)

pca <- prcomp(t(sample_dat))

plot(pca$x[,1], pca$x[,2])
text(pca$x[,1], pca$x[,2], colnames(sample_dat))
summary(pca)

load_val <- sort((abs(pca$rotation[,1])), decreasing = TRUE)

plot(load_val)
abline(h=0.013)

#sample_dat_clust <- hclust(dist(scale(sample_dat)))

#plot(sample_dat_clust)







#filter the my_dat2 matrix
#example from stack overflow of the logic
#d<-d[!(d$A=="B" & d$E==0),]


#application

all_true_fun <- function(datframe, x) {
  all_true <- (abs(datframe$hDC7_ZIKVneg) < x & abs(datframe$hDC10_ZIKVpos) < x
               & abs(datframe$hDC8_ZIKVneg) < x & abs(datframe$hDC10_ZIKVneg) < x
               & abs(datframe$hDC7_ZIKVpos) < x & abs(datframe$hDC8_ZIKVpos) < x)
} 
sample_dat2 <- sample_dat

sample_dat2 <- sample_dat2[!(all_true_fun(sample_dat2, 4)),]


sample_dat_clust2 <- hclust(dist(scale(sample_dat2)))

plot(sample_dat_clust2)


#create a matrix of upregulated and down regulated genes

#using the filtered data

sample_dat3 <- sample_dat2

sample_dat3 <- data.frame(sample_dat2$hDC7_ZIKVpos - sample_dat2$hDC7_ZIKVneg,
                          sample_dat2$hDC8_ZIKVpos - sample_dat2$hDC8_ZIKVneg,
                          sample_dat2$hDC10_ZIKVpos - sample_dat2$hDC10_ZIKVneg, row.names = row.names(sample_dat2))

# positive values is upregulation, negative value is downregulation

pos_count <- length(sample_dat3[sample_dat3 >0,1])

neg_count <- length(sample_dat3[sample_dat3 <0,1])


count_matrix <- c(pos_count, neg_count)
#pplot of ZIKV+ vs ZIKV-
barplot(count_matrix)



#creating a heatmap + dendrogram


##clustering

#load necessary libraries
library("ggplot2")
library("ggdendro")
library("reshape2")
library("grid")
scale.sample_dat2 <- scale(sample_dat2)
matrix.sample_dat2 <- as.matrix(scale.sample_dat2)
#rownames(matrix.sample_dat2) <- row.names(sample_dat2)
sample_dat2_dendro <- as.dendrogram(hclust(d = dist(x = matrix.sample_dat2)))

dendro.plot <- ggdendrogram(data = sample_dat2_dendro, rotate = TRUE)

print(dendro.plot)


dendro.plot <- dendro.plot + theme(axis.text.y = element_text(size = 6))

print(dendro.plot)


##Heatmap preparation


sample.long <- melt(scale.sample_dat2)


heatmap.plot <- ggplot(data = sample.long, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2() +
  theme(axis.text.y = element_text(size = 6))


# Preview the heatmap
print(heatmap.plot)





#putting it together


grid.newpage()
print(heatmap.plot, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot, vp = viewport(x = 0.90, y = 0.445, width = 0.2, height = 1.0))




#Matching the dendrogram to the heatmap


sample.order <- order.dendrogram(sample_dat2_dendro)

sample.long$Var1 <- factor(x = sample.long$Var1,
                           levels = row.names(sample_dat2)[sample.order],
                           ordered = TRUE)

heatmap.plot2 <- ggplot(data = sample.long, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2() +
  theme(axis.text.y = element_text(size = 6))

print(heatmap.plot2)

grid.newpage()
print(heatmap.plot2, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot, vp = viewport(x = 0.90, y = 0.445, width = 0.2, height = 1.0))





#reposition the legend and make the figure publication ready

heatmap.plot2 <- ggplot(data = sample.long, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2() +
  theme(axis.text.y = element_text(size = 6),
        legend.position = "top")

print(heatmap.plot2)



#more finer details

#change the tips to align and remove the double name

colnames(sample.long)[3] <- c("log2")
heatmap.plot2 <- ggplot(data = sample.long, aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = log2)) +
  scale_fill_gradient2() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top",
        text = element_text(size = 8))

#text = element_text(size=20),
#axis.text.x = element_text(angle=90, hjust=1)

print(heatmap.plot2)



#and finally together

grid.newpage()
print(heatmap.plot2, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.05))
print(dendro.plot, vp = viewport(x = 0.90, y = 0.415, width = 0.2, height = 0.96))



