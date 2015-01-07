library("dplyr", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")

WS.orthoMCL.out <- read.delim("~/R-projects/orthoMCL_analysis/WS-orthoMCL-out.txt", header=FALSE, quote="", stringsAsFactors=FALSE)
clusters <- WS.orthoMCL.out[1]
genes <- WS.orthoMCL.out[-1]

temp <- sapply(clusters, function(x) gsub("[:,/,/(,/)]"," ",x))
temp2 <- read.table(text = temp, sep = " ", stringsAsFactors = FALSE)
#get column names
cn <- temp2[1,][seq(from = 4, to = length(temp2[1,]) - 1, by = 2)]

temp3 <- data.frame(temp2[seq(from = 5, to = length(temp2[1,]) - 1, by = 2)],stringsAsFactors = FALSE)
colnames(temp3) <- cn
ortho_tbl<- cbind(temp2[1:3],temp3)
#colnames(ortho_tbl[1:2]) <- c("numSpecies","numGenes")

##PCA
ortho_pca <- prcomp(ortho_tbl[1:1000,4:13])
plot(ortho_pca$rotation)
#whale_shark & lamprey rotations looks very different to the other genomes.  Probably a data quality issue

#USe dplyr to easliy filter the larger table!!
#numerb of genes shared with elephant_shark
ws_es <- filter(ortho_tbl,whale_ > 0,eleph_ > 0, V2 == 2) #num = 216
not_in_ws <- filter(ortho_tbl,whale_ == 0, V2 == 10) #865
not_in_es <- filter(ortho_tbl,eleph_ == 0, V2 == 10) #108
ws_lamp <- filter(ortho_tbl,whale_ > 0,lamp_ > 0, V2 == 2) #num = 104
ws_zebra <- filter(ortho_tbl,whale_ > 0,zebra_ > 0, V2 == 2) #num = 38
ws_coel <- filter(ortho_tbl,whale_ > 0,coel_ > 0, V2 == 2) #num = 165

#number of core genses
core <- filter(ortho_tbl,V2 == 11) #1846
per_core <- filter(ortho_tbl,V2 == 11, V3 == 11)#155

#gene lists
ws_es_genes <- genes[unlist(ws_es[1] +1 ),][2]
