# jkmeans
an R package for jk-means


# input:
Y: N x p (single test)
 or N x p x batchN (multiple tests)

# output:

single test:
mu:  K x p matrix, each row is the estimate for the component mean
w:   p vector, each element is the estimate for the component weight (fixed to 1/J) at this point
zeta: N x K matrix, each row is the affinity for the data to each component
M:   N vector, MAP for clustering

multiple tests:
add one furhter dimension of size batchN to each of the arrays above