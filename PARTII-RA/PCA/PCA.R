
library(tidyverse)

## Genetic distances between individuals
# actual command:
system("/Users/ngoc/plink/plink --bfile /Users/ngoc/All/GENIOMHE1/subjects/NGS/TP_DATA_PROG/II.a.Plink/output160125/lndMkMAFHWE --distance-matrix --out PCAdata")

# pretty command:
system("plink --bfile lndMkMAFHWE --distance-matrix --out PCAdata")


## Load data
dist_populations<-read.table("PCAdata.mdist",header=F)

## Getting countries of origin
country <- data.frame(country=read.table("PCAdata.mdist.id")[,1])
country$country <- substr(country$country, 1, nchar(country$country) - 3)
country$country <- ifelse(country$country == "HCB", "CHB", country$country)

## Perform PCA using the cmdscale function
mds_populations <- cmdscale(dist_populations, eig = TRUE, k = 5)

## Extract the eigen vectors
eigenvec_populations <- cbind(country,mds_populations$points)

## Proportion of variation captured by each eigen vector
eigen_percent <- round(((mds_populations$eig)/sum(mds_populations$eig))*100,2)

ggplot(data = eigenvec_populations) +
  geom_point(mapping = aes(x = `1`, y = `2`, color = country), show.legend = T, size=2) + 
  labs(x = paste0("Principal component 1 (",eigen_percent[1]," %)"),
       y = paste0("Principal component 2 (",eigen_percent[2]," %)")) + 
  scale_color_manual(values = c("red2", "blue2")) + 
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
    axis.line.x= element_blank(), 
    axis.line.y= element_blank(), 
    panel.grid.major = element_line(color="grey97"),  
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12, color = "black"), 
    axis.title = element_text(size = 12, color = "black"))
