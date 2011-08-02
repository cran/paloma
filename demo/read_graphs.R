library(paloma)

# Get the path of training files
path  <- system.file("networks",package="paloma")
fname <- paste( path, '/', 'example.txt', sep='' )

# Reading a text file with the following format
#  src -> dst
edges <- as.matrix( read.table( fname ) )

# Get over-represented motif according to an Erdos-Renyi model
out <- getExceptional( edges, 3, 5, pvalue=0.05 )

readline("continue :")

# print over-represented motifs
out

readline("continue :")

# The same with 'igraph' package
cat("This part of the demonstration need igraph package)")
library(igraph)

G     <- read.graph( fname )
# Remove redondant edges and loops
G     <- simplify( G, remove.multiple = TRUE, remove.loops = TRUE )
edges <- get.edgelist( G )

# Get over-represented motif according to an Erdos-Renyi model
out   <- getExceptional( edges, 3, 5, pvalue=0.05 )

# plot over-represented motifs
plot( out )

