# Simple Network Visualization
# Create a basic network plot with ggplot2 instead of ggraph

# Create output directory if it doesn't exist
dir.create("inst/example_plots/images", recursive = TRUE, showWarnings = FALSE)

# Load required libraries
library(microFGT)
library(ggplot2)
library(dplyr)
library(igraph)

# Generate a simple correlation matrix
set.seed(42)
n <- 10
cor_mat <- matrix(runif(n*n, -1, 1), nrow=n, ncol=n)
diag(cor_mat) <- 1
cor_mat[lower.tri(cor_mat)] <- t(cor_mat)[lower.tri(cor_mat)]

# Set node names
node_names <- paste0("Taxon_", 1:n)
rownames(cor_mat) <- node_names
colnames(cor_mat) <- node_names

# Create node dataframe
nodes <- data.frame(
  id = node_names,
  type = sample(c("Firmicutes", "Bacteroidetes", "Actinobacteria"), n, replace = TRUE),
  abundance = runif(n, 0, 1)
)

# Create edges dataframe
edges <- data.frame(
  from = character(),
  to = character(), 
  weight = numeric(),
  type = character()
)

# Filter significant correlations
threshold <- 0.6
for (i in 1:n) {
  for (j in i:n) {
    if (i != j && abs(cor_mat[i,j]) >= threshold) {
      edges <- rbind(edges, data.frame(
        from = node_names[i],
        to = node_names[j],
        weight = abs(cor_mat[i,j]),
        type = ifelse(cor_mat[i,j] > 0, "positive", "negative")
      ))
    }
  }
}

# Create the network with igraph
g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)

# Get layout
set.seed(42)
layout <- layout_with_fr(g)
layout_df <- data.frame(x = layout[,1], y = layout[,2])
layout_df$id <- nodes$id

# Prepare node data for plotting
node_data <- nodes
node_data$x <- layout_df$x
node_data$y <- layout_df$y
node_data$size <- 3 + 5 * nodes$abundance

# Prepare edge data for plotting
edge_data <- edges
edge_data$x1 <- layout_df$x[match(edge_data$from, layout_df$id)]
edge_data$y1 <- layout_df$y[match(edge_data$from, layout_df$id)]
edge_data$x2 <- layout_df$x[match(edge_data$to, layout_df$id)]
edge_data$y2 <- layout_df$y[match(edge_data$to, layout_df$id)]
edge_data$width <- 0.5 + 2.5 * edge_data$weight

# Create network plot
network_plot <- ggplot() +
  # Add edges
  geom_segment(data = edge_data, 
               aes(x = x1, y = y1, xend = x2, yend = y2,
                   color = type, linewidth = width),
               alpha = 0.7) +
  # Add nodes
  geom_point(data = node_data,
             aes(x = x, y = y, size = size, fill = type),
             shape = 21, color = "white") +
  # Add labels
  geom_text(data = node_data,
            aes(x = x, y = y, label = id),
            size = 3, vjust = -1.5) +
  # Scales
  scale_color_manual(values = c("positive" = "forestgreen", "negative" = "firebrick")) +
  scale_fill_brewer(palette = "Set1") +
  scale_size_identity() +
  scale_linewidth_identity() +
  # Theme
  theme_void() +
  labs(
    title = "Taxa Co-occurrence Network",
    subtitle = "Edge color indicates positive (green) or negative (red) correlation",
    color = "Correlation",
    fill = "Phylum"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  coord_fixed()

# Save the network plot
ggsave("inst/example_plots/images/taxa_network_simple.png", network_plot, width = 10, height = 8, dpi = 300)

# Print success message
cat("Simple network plot created successfully!\n")
cat("Plot saved to inst/example_plots/images/\n")