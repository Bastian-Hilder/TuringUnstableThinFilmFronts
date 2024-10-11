# Fast-moving pattern interfaces close to a Turing instability in an asymptotic model for the three-dimensional Bénard–Marangoni problem

This repository contains the supplementary code files related to the article "Fast-moving pattern interfaces close to a Turing instability in an asymptotic model for the three-dimensional Bénard–Marangoni problem" by Bastian Hilder and Jonas Jansen. The preprint can be found at [arXiv:2410.02708](https://arxiv.org/abs/2410.02708).

### Abstract

We study the bifurcation of planar patterns and fast-moving pattern interfaces in an asymptotic long-wave model for the three-dimensional Bénard–Marangoni problem, which is close to a Turing instability. We derive the model from the full free-boundary Bénard–Marangoni problem for a thin liquid film on a heated substrate of low thermal conductivity via a lubrication approximation. This yields a quasilinear, fully coupled, mixed-order degenerate-parabolic system for the film height and temperature. As the Marangoni number $M$ increases beyond a critical value $M^*$, the pure conduction state destabilises via a Turing(–Hopf) instability. Close to this critical value, we formally derive a system of amplitude equations which govern the slow modulation dynamics of square or hexagonal patterns. Using center manifold theory, we then study the bifurcation of square and hexagonal planar patterns. Finally, we construct planar fast-moving modulating travelling front solutions that model the transition between two planar patterns. The proof uses a spatial dynamics formulation and a center manifold reduction to a finite-dimensional invariant manifold, where modulating fronts appear as heteroclinic orbits. These modulating fronts facilitate a possible mechanism for pattern formation, as previously observed in experiments.

### Supplementary Material

The Supplementary Material for the article containing the explicit expressions for the coefficients and the eigenvalues of the linearisation about the mixed modes is available [here](https://github.com/Bastian-Hilder/TuringUnstableThinFilmFronts/blob/4aee854da1a4afe367ea77802eb0ad60ea9047d2/supplement.pdf).

### Videos

Here are three examples of moving pattern interfaces constructed in the article. Please note that clicking the thumbnails below will redirect you to the videos on YouTube. First, an invasion of the pure conduction state by a hexagon pattern.

[![Modfront-H-to-T](https://github.com/Bastian-Hilder/TuringUnstableThinFilmFronts/blob/5409b90a7e84a7db7ed597b94de9a45d9700d70b/pictures/Modfront-H-to-T.jpg)](https://www.youtube.com/watch?v=i0xGN3Wy6gs)

Second, an invasion of down-hexagons by up-hexagons.

[![Modfront-up-H-to-down-H](https://github.com/Bastian-Hilder/TuringUnstableThinFilmFronts/blob/32ea84d40611dd648698a307bdb9f5180183409b/pictures/Modfront-up-H-to-down-H.jpg)](https://www.youtube.com/watch?v=_b2uqaCqw2M)

Third, a two-stage invasion process consisting of a primary invasion of the pure conduction state by roll waves followed by a secondary invasion of the roll waves by a hexagon pattern.

[![Modfront-H-to-R-to-T](https://github.com/Bastian-Hilder/TuringUnstableThinFilmFronts/blob/32ea84d40611dd648698a307bdb9f5180183409b/pictures/Modfront-H-to-R-to-T.jpg)](https://www.youtube.com/watch?v=il-8Szh2VOk)

The videos are generated with Matlab using the code provided in the folder "videos". More videos are available [in this YouTube playlist](https://www.youtube.com/playlist?list=PLc-ENVP1CPzfEf60dixeT_h6C5PVBvqLB).

### Coefficients

The coefficients for the formal amplitude equations (and the reduced equations on the center manifolds) are calculated using the Mathematica Notebook in the folder "Coefficient-Mathematica". The output of the notebook is also available in PDF and HTML formats.

Please note that the notebook can be read using the Wolfram Player, which is freely available [on the Wolfram website](https://www.wolfram.com/player/?source=nav).

### Coefficient plots

The plots of the coefficients are generated using Python using the files contained in the folder "Coefficients-Plots". The Python scripts read in the explicit expressions provided by the Mathematica Notebook, also available in "Coefficients-Plots".

### Phase plane

The phase plane diagrams for the heteroclinic orbits on the reduced system are also generated by the streamplot function provided by Mathematica. The corresponding Mathematica Notebooks for both the square and hexagonal lattice are available in the folder "Phase planes". A pdf version of the notebook outputs can also be found there.
