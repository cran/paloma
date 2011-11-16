pkgname <- "paloma"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('paloma')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("computeModelMean")
### * computeModelMean

flush(stderr()); flush(stdout())

### Name: computModelMean
### Title: Computes the expected occurrence number of a motif
### Aliases: computeModelMean
### Keywords: stats models cluster

### ** Examples

## Load the data set
data(example)

# Get over-represented motifs under a statistical
# block model with two classes
out <- getExceptional( example$graph, 3, 5, example$classes,
          example$pi, pvalue=0.001 )


# Get the first over-represented motif
m  <- out[[1]]$motifs$adjacency[[1]]

# Compute the expected counts of the motif occurrences
computeModelMean( example$classes, example$pi, m)




cleanEx()
nameEx("countAllMotifs")
### * countAllMotifs

flush(stderr()); flush(stdout())

### Name: countAllMotifs
### Title: Counts the number of all motifs having k nodes in a given graph
### Aliases: countAllMotifs
### Keywords: graphs

### ** Examples

## Load the data set
data(example)

# Count all motifs with 3 nodes in a graph
out <- countAllMotifs( example$graph, 3, directed=TRUE )

# Number of discovery motifs
motif_nbr <- length( out$adjacency)

# Getting the fifth motif
motif <- matrix( out$adjacency[[5]], 3, 3  )

# Number of occurrences of this motif
nb <- out$count[5]



cleanEx()
nameEx("data-example")
### * data-example

flush(stderr()); flush(stdout())

### Name: example
### Title: Graph example for package tests
### Aliases: example
### Keywords: datasets

### ** Examples

## load the data set
data(example)

out <- getExceptional( example$graph, 3, 4,
                       example$classes, example$pi, 0.001 )



cleanEx()
nameEx("data-yeast")
### * data-yeast

flush(stderr()); flush(stdout())

### Name: yeast
### Title: Graph example for package tests
### Aliases: yeast
### Keywords: datasets

### ** Examples

## load the data set
data(yeast)

# Get exceptional motifs
out <- getExceptional( yeast$graph, 3, 4,
                       yeast$classes, yeast$pi, 0.001, directed=TRUE )

# Show the exceptional motifs
plot( out )




cleanEx()
nameEx("getAllExtensions")
### * getAllExtensions

flush(stderr()); flush(stdout())

### Name: getAllExtensions
### Title: Searchs all occurrences and extensions in a network.
### Aliases: getAllExtensions
### Keywords: graphs

### ** Examples


## Load the data set
data(example)

# Get all occurrences and extensions of motifs
# containing 3 nodes in a graph.
out <- getAllExtensions( example$graph, 3, directed=TRUE )

# Number of nodes
out[[1]]$k

# Considering the second motif
out[[1]]$motifs$adjacency[[2]]
motif <- out[[1]]$motifs$m.mp[[2]]

# There is 3 deletion class for this motifs
# Select the third one
# D class (node list)
motif[[3]]$DelClass  # deletion class (node list)

# Second occurrence in this class (k-1 nodes)
motif[[3]]$Occur[[2]]

# Possible extensions of this occurrence
# in the same deletion class
motif[[3]]$Ext[[2]]



cleanEx()
nameEx("getCanonic")
### * getCanonic

flush(stderr()); flush(stdout())

### Name: getCanonic
### Title: Compute the canonical adjacency matrix
### Aliases: getCanonic
### Keywords: graphs

### ** Examples

# Define a 3 node motif
motifV <-  matrix( c(0, 1, 1, 1, 0, 0, 1, 0, 0), 3, 3) 

# Compute canonical form
cf <- getCanonic( motifV )




cleanEx()
nameEx("getExceptional")
### * getExceptional

flush(stderr()); flush(stdout())

### Name: getExceptional
### Title: Seeks for exceptional motifs under a statistical block model.
### Aliases: getExceptional
### Keywords: graphs models cluster

### ** Examples


## Load the data set
data(example)

# Get couple (motif, deletion class) over-represented
# under a statistical block model with two classes
# The connectivity matrix is stored in "example$pi"
# and the mapping node/class is stored in "example$classes"
out <- getExceptional( example$graph, 3, 5, example$classes,
          example$pi, pvalue=0.001 )

# The same without the hierarchical filter
out <- getExceptional( example$graph, 3, 5, example$classes,
          example$pi, pvalue=0.001, h.filter=FALSE )

# Print exceptional motifs
out

# Without classes
out <- getExceptional( example$graph, 3, 5, pvalue=0.001)

#
# Use "mixer" R-package to cluster the graph
#
data(yeast)

# Clustering the graph (obtain the 'pi' matrix and the vertex classes)
res <- mixer( t( yeast$graph ), method="bayesian", qmin=2, qmax=10, directed=TRUE)

out <- getExceptional( res, 3, 5, pvalue=0.001)
out




cleanEx()
nameEx("getExtensions")
### * getExtensions

flush(stderr()); flush(stdout())

### Name: getExtensions
### Title: Search occurrences and extensions in a network for a given motif
### Aliases: getExtensions
### Keywords: graphs

### ** Examples


## Load the data set
data(example)

# Get couple (motif, deletion class) over-represented
# under a statistical block model with two classes
# The connectivity matrix is stored in "example$pi"
# and the mapping node/class is stored in "example$classes"
out <- getExceptional( example$graph, 3, 5, example$classes,
          example$pi, pvalue=0.001 )

# Print exceptionnal motifs
out


# Get all positions (with k-1 nodes) with the deletion class of node 0
# for the 11th exceptional motif of size k=5
u <- getExtensions( example$graph, out[[3]]$motifs$adjacency[[11]], 0 )

# Print occurrences
u




cleanEx()
nameEx("plot.exceptional.motif")
### * plot.exceptional.motif

flush(stderr()); flush(stdout())

### Name: plot.paloma.exceptional
### Title: Plots the over-represented motif list
### Aliases: plot.paloma.exceptional
### Keywords: plot

### ** Examples


## Load the data set
data(example)

# Get couple (motif, deletion class) over-represented
# under a statistical block model with two classes
# The  the connectivity matrix is stored in "example$pi"
# and the mapping node/class is stored in "example$classes"
out <- getExceptional( example$graph, 3, 5, example$classes,
          example$pi, pvalue=0.001 )

# Get the all the over-represented motifs (motif, deletion class)
plot( out )




cleanEx()
nameEx("print.exceptional.motif")
### * print.exceptional.motif

flush(stderr()); flush(stdout())

### Name: print.paloma.exceptional
### Title: Prints the over-represented motifs list
### Aliases: print.paloma.exceptional
### Keywords: print

### ** Examples


## Load the data set
data(example)

# Get couple (motif, deletion class) over-represented
# under a statistical block model with two classes
# The  the connectivity matrix is store in "example$pi"
# and the mapping node/class is store in "example$classes"
out <- getExceptional( example$graph, 3, 5, example$classes,
          example$pi, pvalue=0.001 )

# Get the all the over-represented motifs (motif, deletion class)
print( out )
# or
out

# Display  the couple (motif, deletion class) for the first motif
print( out, k=1, m=1 )

# The same with only the second deletion class
print( out, k=1, m=1, mp=2 )




cleanEx()
nameEx("print.occurrences")
### * print.occurrences

flush(stderr()); flush(stdout())

### Name: print.paloma.occurrences
### Title: Prints the motifs occurrences and their extensions
### Aliases: print.paloma.occurrences
### Keywords: print

### ** Examples

## Load data set
data( example )

# Get the all the extensions for a motif size 
occ <- getAllExtensions( example$graph, 3, directed=TRUE )

# Display extensions (occurrences)
print( occ )
# or
occ

# Display the occurrences and extensions of the first motif
# Two deletion classes are displayed 
print( occ, k=1, m=1 )

# Display the occurrences and extensions of the second deletion class of
# the first motif  
print( occ, k=1, m=1, mp=2 )




### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
