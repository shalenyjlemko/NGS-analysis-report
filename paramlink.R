setwd('C:/Users/denis/OneDrive/Desktop/Uni/NGS/Project/TP_DATA_PROG/I.a.Paramlink')
library(paramlink)
fam = read.table('fam.txt')

fam[1:5, 1:10]

x = linkdat(fam)

summary(x)

# generating family tree figure
plot(x, marker=1)

library(kinship2)
library(igraph)
library(ggraph)
library(tidyverse)
library(dplyr)


id <- fam[ , 2]
dadid <- fam[ , 3]
momid <- fam[ , 4]
sex <- fam[ , 5]
affected <- fam[ , 6]

genotypes <- fam[ , 7:32]

mark1_gt <- paste(fam[ ,7], fam[, 8], sep = "/")

ped <- pedigree(id = id, dadid = dadid, momid = momid, sex = sex, affected = affected)

colors <- ifelse(fam$V6 == 2, "lightgreen", ifelse(fam$V6 == 1, "black", "grey"))

plot(ped, symbolsize = 1.3, cex = 0.2, col = colors, 
     packed = TRUE,
     align = TRUE,
     width = 30,
     branch = 0.9,
     #subregion = TRUE,
     pconnect = 0,
     packed.spacing = 2,
     main = "Pedigree with 5 generations"
     )

text(ped$plot$x, ped$plot$y, labels = mark1_gt, cex = 0.6, col = "blue")
# trying for no result



plot(x, marker = 1, symbolsize = 1.1, col = colors, cex = 0.4, pconnect = 0, branch = 0.2, width = 50, main = "Pedigree with Genotypes")




# Save plot as PNG with higher resolution and larger dimensions
png("pedigree_plot.png", width = 3840, height = 2160, res = 500)  # Customize width, height, and resolution
plot(x, marker = 1, symbolsize = 1.2, col = colors, cex = 0.5, pconnect = 0, branch = 0.2, width = 50, main = "Pedigree with Genotypes")

dev.off()  # Close the device to save the file

# running lod score

xdom = setModel(x, model=1, penetrances = c(0.00001, 1, 1), dfreq = 0.00001)

result_dom = lod(xdom, theta=c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5)) 
result__dom_df = as.data.frame(result_dom)
write.csv(result__dom_df, "lod_results_dom.csv", row.names = TRUE)

result_dom2 = lod(xdom, theta=c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5)) 
result__dom2_df = as.data.frame(result_dom2)
write.csv(result__dom2_df, "lod_results_dom2.csv", row.names = TRUE)


zmax_minus_11 <- 6.6738 

lod_max <- lod(xdom, theta = 0)

theta_range <- seq(0.06, 0.2, by = 0.001) 
i_range <- seq(1, 10, by = 1)
for (i in i_range) {
    zmax_minus_11 <- lod_max[, i] - 1  # Calculate zmax - 1 for the marker
    for (theta in theta_range) {
        result_dom_ex <- lod(xdom, theta = theta)
        
        # Check if the LOD score drops below zmax - 1
        if (result_dom_ex[, i] < zmax_minus_11) {
            previous_theta <- theta - 0.001  # Backtrack to the previous theta
            previous_result <- lod(xdom, theta = previous_theta)
            
            print(paste("Marker", i, "Lod-score:", previous_result[, i], "theta = ", previous_theta))
            break # Stop searching for this marker once the condition is met
        }
    }
}
lod(xdom, marker=c(5,7,8,12), theta='max') 

xdom5=modifyMarker(xdom,marker = 5, afreq = c(0.1, 0.1, 0.1, 0.7))
lod(xdom5, marker=5, theta=c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5))

xrec=setModel(x, model=2, penetrances=c(0.00001,0.00001, 1), dfreq=0.00001)
result_rec = lod(xrec)

result_rec_df = as.data.frame(result_rec)
write.csv(result_rec_df, "lod_results_rec.csv", row.names = TRUE)

result_rec2 = lod(xrec, theta=c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5))
result_rec2_df = as.data.frame(result_rec2)
write.csv(result_rec2_df, "lod_results_rec2.csv", row.names = TRUE)
