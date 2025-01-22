setwd('C:/Users/denis/OneDrive/Desktop/Uni/NGS/Project/TP_DATA_PROG/I.a.Paramlink')
library(paramlink)
fam = read.table('fam.txt')


x = linkdat(fam)

xdom = setModel(x, model=1, penetrances = c(0.00001, 1, 1), dfreq = 0.00001)

lod_max <- lod(xdom, theta = 0) # Get the maximum LOD scores at theta = 0

theta_range <- seq(0.05, 0.2, by = 0.001)
i_range <- seq(1, 10, by = 1)
# Precompute LOD scores for all theta values
theta_lod_scores <- t(sapply(theta_range, function(theta) lod(xdom, theta = theta)))

# Loop for markers 1–10
for (i in i_range) {
    zmax_minus_11 <- lod_max[, i] - 1  # Calculate zmax - 1 for the marker
    if (lod_max[, i] <= 3) { 
        # Inconclusive case: add NA for lower and upper bounds
        ci_bounds <- rbind(ci_bounds, data.frame(Marker = i, Lower = NA, Upper = NA, Zmax = lod_max[, i], Zmax_Theta = 0))
        next
    }
    
    for (j in seq_along(theta_range)) {
        result_dom_ex <- theta_lod_scores[j, ]
        
        # Check if the LOD score drops below zmax - 1 or 3
        if (result_dom_ex[i] < zmax_minus_11 || result_dom_ex[i] < 3) {
            previous_theta <- theta_range[j] - 0.001  # Backtrack to the previous theta
            ci_bounds <- rbind(ci_bounds, data.frame(Marker = i, Lower = 0, Upper = previous_theta, Zmax = lod_max[, i], Zmax_Theta = 0))
            break
        }
    }
}

# Additional loop for markers 11–13
i_range_2 <- seq(11, 13, by = 1)
theta_range <- seq(0.02, 0.2, by = 0.001)

for (i in i_range_2) {
    # Find the maximum LOD score and its corresponding theta
    lod_scores <- theta_lod_scores[, i]
    zmax <- max(lod_scores)  # Find the maximum LOD score for this marker
    zmax_theta <- theta_range[which.max(lod_scores)]  # Find the theta at which Zmax occurs
    
    if (zmax < 3) {
        # Inconclusive case: add NA for lower and upper bounds
        ci_bounds <- rbind(ci_bounds, data.frame(Marker = i, Lower = NA, Upper = NA, Zmax = zmax, Zmax_Theta = zmax_theta))
        next
    }
    
    # Find lower bound
    lower_theta <- NA  # Initialize lower bound
    for (j in seq_along(theta_range)) {
        result_dom_ex <- theta_lod_scores[j, ]
        if (result_dom_ex[i] > zmax - 1 || result_dom_ex[i] > 3) {
            lower_theta <- theta_range[j] - 0.001
            break
        }
    }
    
    # Find upper bound
    upper_theta <- NA  # Initialize upper bound
    for (j in rev(seq_along(theta_range))) { # Search in reverse for the upper bound
        result_dom_ex <- theta_lod_scores[j, ]
        if (result_dom_ex[i] > zmax - 1 || result_dom_ex[i] > 3) {
            upper_theta <- theta_range[j] + 0.001
            break
        }
    }
    
    # Add both bounds and Zmax information to ci_bounds
    ci_bounds <- rbind(ci_bounds, data.frame(Marker = i, Lower = lower_theta, Upper = upper_theta, Zmax = zmax, Zmax_Theta = zmax_theta))
}

# Print CI bounds
print(ci_bounds)
