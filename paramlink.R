setwd('C:/Users/denis/OneDrive/Desktop/Uni/NGS/Project/TP_DATA_PROG/TP_DATA_PROG/I.a.Paramlink')
library(paramlink)
fam = read.table('fam.txt')

fam[1:5, 1:10]

x = linkdat(fam)

summary(x)

plot(x, marker=1)

xdom = setModel(x, model=1, penetrances = c(0.00001, 1, 1), dfreq = 0.00001)

result = lod(xdom, theta=c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5)) 

result_df = as.data.frame(result)

write.csv(result_df, "lod_results.csv", row.names = TRUE)

zmax_minus_11 <- 6.6738 

theta_range <- seq(0.06, 0.08, by = 0.001) 

for (theta in theta_range) {
    # Compute the LOD score for the current theta
    result <- lod(xdom, theta = theta)
    
    # Check if the first column matches the zmax - 1 value
    if (result[, 1] > zmax_minus_11) {
        print(paste("Theta at zmax - 1 is:", theta))
        print(result)
        break
    }
}

lod(xdom, marker=c(5,7,8,12), theta='max') 

xdom5=modifyMarker(xdom,marker = 5, afreq = c(0.1, 0.1, 0.1, 0.7))
lod(xdom5, marker=5, theta=c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5))

xrec=setModel(x, model=2, penetrances=c(0.00001,0.00001, 1), dfreq=0.00001)
result_rec = lod(xrec)

result_rec_df = as.data.frame(result_rec)
write.csv(result_rec_df, "lod_results_rec.csv", row.names = TRUE)

result_rec2 = lod(xrec, theta=c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5))
result_rec2_df = as.data.frame(result_rec2)
write.csv(result_rec2_df, "lod_results_rec2.csv", row.names = TRUE)
