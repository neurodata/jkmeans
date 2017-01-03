1. Model-based clustering of multi-dimensional data (e.g. large dimension, curves, shapes) tends to produce clusters are similar to each other.

2. The regularization on controling the inter-cluster distance is theoretically elegant, but faces practical computational issue and is sensitive to over-regularization that pushes cluster centers too far.

3. There is an undiscovered connection between regularized clustering assignment and regularized inter-cluster distance (based on current theorm 1). The former does not have sensitivity issue  (probably a theory with an upper bound on the ratio based on inter-/intra-clustering distance) and is easy to apply.

4. Simulation results are shown from simple univariate to complicated setting like curve clustering, compared against GMM, k-means.

5. Some applications like eigenface, MNIST.

Additional thought: the complete assignment may be an overly strong regularization, to mimic the inter-cluster distance regularization, there might be ways to adjust the regularization.
