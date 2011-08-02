library(paloma)

# Read data set
data(yeast)

# Get exceptional motifs
out <- getExceptional( yeast$graph, 2, 4,
                       yeast$classes, yeast$pi, 0.001, directed=TRUE )

# Show the exceptional motifs
plot( out )
