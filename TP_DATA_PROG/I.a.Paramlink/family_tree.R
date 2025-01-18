library(paramlink)
library(kinship2)
library(igraph)
library(ggraph)
library(tidyverse)
library(dplyr)

setwd('C:/Users/denis/OneDrive/Desktop/Uni/NGS/Project/TP_DATA_PROG/I.a.Paramlink')

fam = read.table('fam.txt')

x = linkdat(fam)


fam <- fam %>%
    rename(
        famID  = V1,
        id     = V2,
        dadID  = V3,
        momID  = V4,
        sex    = V5,
        status = V6
        
    )

pedigree <- as.data.frame(x$pedigree)
colnames(pedigree) <- c("id", "dadID", "momID", "sex", "status")

# Add generation information
pedigree <- pedigree %>%
    mutate(
        # Identify founders (no parents in the data)
        founder = (dadID == 0 & momID == 0),
        
        # Map sex to labels
        sex_label = ifelse(sex == 1, "Male", "Female"),
        
        # Map status to labels
        status_label = case_when(
            status == 2 ~ "Affected",
            status == 1 ~ "Unaffected",
            TRUE        ~ "Unknown"
        ),
        
        # Initialize generation for founders
        generation = ifelse(founder, 0, NA_integer_)
    )

# Assign generations using a BFS approach
changed <- TRUE
while (changed) {
    changed <- FALSE
    for (i in seq_len(nrow(pedigree))) {
        if (is.na(pedigree$generation[i])) {
            dad_gen <- pedigree$generation[pedigree$id == pedigree$dadID[i]]
            mom_gen <- pedigree$generation[pedigree$id == pedigree$momID[i]]
            new_gen <- min(dad_gen, mom_gen, na.rm = TRUE) + 1
            if (!is.infinite(new_gen) && is.na(pedigree$generation[i])) {
                pedigree$generation[i] <- new_gen
                changed <- TRUE
            }
        }
    }
}

# Step 3: Ensure couples share the same generation
# Identify couples (pairs who have children together)
couples <- pedigree %>%
    filter(dadID != 0, momID != 0) %>%
    select(dadID, momID) %>%
    distinct()

for (i in seq_len(nrow(couples))) {
    dad <- couples$dadID[i]
    mom <- couples$momID[i]
    dad_gen <- pedigree$generation[pedigree$id == dad]
    mom_gen <- pedigree$generation[pedigree$id == mom]
    # If one has a known generation, propagate it to the other
    if (!is.na(dad_gen) && is.na(mom_gen)) {
        pedigree$generation[pedigree$id == mom] <- dad_gen
    } else if (is.na(dad_gen) && !is.na(mom_gen)) {
        pedigree$generation[pedigree$id == dad] <- mom_gen
    } else if (!is.na(dad_gen) && !is.na(mom_gen) && dad_gen != mom_gen) {
        # If generations differ, adjust the higher one
        gen_to_update <- max(dad_gen, mom_gen)
        pedigree$generation[pedigree$id == dad] <- gen_to_update
        pedigree$generation[pedigree$id == mom] <- gen_to_update
    }
}



edges_parent_child <- bind_rows(
    pedigree %>%
        filter(dadID != 0) %>%
        transmute(from = dadID, to = id, relationship = "parent"),
    
    pedigree %>%
        filter(momID != 0) %>%
        transmute(from = momID, to = id, relationship = "parent")
)

edges_couple <- pedigree %>%
    filter(dadID != 0, momID != 0) %>%
    select(dadID, momID) %>%
    distinct() %>%
    transmute(
        from = pmin(dadID, momID),
        to = pmax(dadID, momID),
        relationship = "couple"
    ) %>%
    distinct()

edges <- bind_rows(edges_parent_child, edges_couple)

# Create graph
g <- graph_from_data_frame(
    d = edges,
    directed = TRUE,  # Use directed for parent-child relationships
    vertices = pedigree %>%
        select(id, sex_label, status_label, generation)
)

ggraph(g, layout = 'sugiyama') +  # 'sugiyama' layout is good for hierarchies
    geom_edge_link(aes(color = relationship), arrow = arrow(length = unit(2, 'mm'))) +
    geom_node_point(aes(color = status_label, shape = sex_label), size = 5) +
    geom_node_text(aes(label = id), vjust = -1.2) +
    scale_shape_manual(values = c("Male" = 15, "Female" = 16)) +
    scale_color_manual(values = c("Affected" = "red", "Unaffected" = "green", "Unknown" = "grey")) +
    theme_void() 
