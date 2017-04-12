# DP_count
Privacy Preserving Medical Count Queries

This repository contains the code for our work on medical count queries (see the manuscript, to be added once available).

In particular, the code in countQuery.py implements two methods for releasing differentially private count queries: Vinterbo et al's method, and our modified method. These are implemented in the NaiveExp and SmartExp functions.

In both cases, the function takes in the parameters for the loss function (see the manuscript) and returns a noisy estimate of the count (the y parameter).

The Risk method estimates the expected Risk for a given choice of parameters.

The remaining functions were used to produce figures for the paper.
