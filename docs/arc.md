1. Model-based clustering of multi-dimensional data (e.g. large dimension, curves, shapes) tends to produce clusters are similar to each other.

2. The regularization on controling the inter-cluster distance is theoretically elegant, but faces practical computational issue and is sensitive to over-regularization that pushes cluster centers too far.

3. There is an undiscovered connection between regularized clustering assignment and regularized inter-cluster distance (based on current theorm 1). The former does not have sensitivity issue  (probably a theory with an upper bound on the ratio based on inter-/intra-clustering distance) and is easy to apply.

theorem 1:
the complete membership leads to larger distance between the centers, when there is an overlap between centers

theorem 2 (conjecture):
the distance caused by the automatic repulsion is upper bounded by the distance between the centers learned from clustering without the overlapped part.

4.

sim1:
Univariate: two clusters, each with size n1 and n2, with m interfering noisy data in the middle. Setting n1=n2*2, keep m fixed at moderate and increase n2. Plot the misclassification rate, NMI & center estimates.

sim2:
Multivariate: repeat sim1, except each cluster is p dimensional, increase p as well (don't need new theory here, just to show the proposed model isn't just limited to toy setting)

5. 
application1:
eigenface
application2:
MNIST

Additional thought: the complete assignment may be an overly strong regularization, to mimic the inter-cluster distance regularization, there might be ways to adjust the regularization.
