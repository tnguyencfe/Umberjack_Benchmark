"""
one-off script to check why increasing ave robinson foulds distance between window tree and true tree
increases the accuracy of umberjack  in find_covar_accuracy_multisample.R
"""

# TestCase:  make a window where all the sequences are exactly like the full population.  (not reads, the full sequence)
#   - do umberjack on window
# Expected:  we expect that the window tree won't have the exact same topology due to greedy fasttree.
# But we expect the robinson foulds distance to be minimal.

# TestCase:  make a window where all the sequences are exactly like a slice of full population.
#   - do umberjack on window
# Expected:  we expect that the window tree won't have the exact same topology due to greedy fasttree.
#   But we expect the robinson foulds distance to be minimal, but more than the robinson foulds when window matched entire full population
#   The umberjack accuracy should decrease

# TestCase:  make a window where all the sequences are exactly like the full population.  (not reads, the full sequence)
#   - do umberjack on window
# Expected:  we expect that the window tree won't have the exact same topology due to greedy fasttree.
# But we expect the robinson foulds distance to be minimal.



