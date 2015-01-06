WS.orthoMCL.out <- read.delim("~/R-projects/orthoMCL_analysis/WS-orthoMCL-out.txt", header=FALSE, quote="", stringsAsFactors=FALSE)
clusters <- WS.orthoMCL.out[1]
temp <- sapply(clusters, function(x) gsub("[:,/,/(,/)]"," ",x))
temp2 <- read.table(text = temp, sep = " ", stringsAsFactors = FALSE)
#get column names
cn <- temp2[1,][seq(from = 4, to = length(temp2[1,]) - 1, by = 2)]

temp3 <- data.frame(temp2[seq(from = 5, to = length(temp2[1,]) - 1, by = 2)],stringsAsFactors = FALSE)
colnames(temp3) <- cn
ortho_tbl<- cbind(temp2[2:3],temp3)
#colnames(ortho_tbl[1:2]) <- c("numSpecies","numGenes")

##PCA
ortho_pca <- prcomp(ortho_tbl[3:13])
plot(ortho_pca$rotation)
#whale_shark rotations looks very different to the other genomes.  Probably a data quality issue