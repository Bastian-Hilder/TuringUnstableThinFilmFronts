# Fast-moving pattern interfaces close to a Turing instability in an asymptotic model for the three-dimensional Bénard–Marangoni problem

In this repository you find the supplementary code files related to the article "Fast-moving pattern interfaces close to a Turing instability in an asymptotic model for the three-dimensional Bénard–Marangoni problem" by Bastian Hilder and Jonas Jansen. The preprint can be found at [arXiv:???](https://arxiv.org/abs/??).

### Abstract

We study the bifurcation of planar patterns and fast-moving pattern interfaces in an asymptotic long-wave model for the three-dimensional Bénard–Marangoni problem, which is close to a Turing instability. We derive the model from the full free-boundary Bénard–Marangoni problem for a thin liquid film on a heated substrate of low thermal conductivity via a lubrication approximation. This yields a quasilinear, fully coupled, mixed-order degenerate-parabolic system for the film height and temperature. As the Marangoni number $M$ increases beyond a critical value $M^*$, the pure conduction state destabilises via a Turing(–Hopf) instability. Close to this critical value, we formally derive a system of amplitude equations which govern the slow modulation dynamics of square or hexagonal patterns. Using center manifold theory, we then study the bifurcation of square and hexagonal planar patterns. Finally, we construct planar fast-moving modulating travelling front solutions that model the transition between two planar patterns. The proof uses a spatial dynamics formulation and a center manifold reduction to a finite-dimensional invariant manifold, where modulating fronts appear as heteroclinic orbits. These modulating fronts facilitate a possible mechanism for pattern formation, as previously observed in experiments.

### Supplementary Material

The Supplementary Material for the article containing the explicit expressions for the coeffients and the eigenvalues of the linearisation about the mixed modes is available [here](https://github.com/Bastian-Hilder/TuringUnstableThinFilmFronts/blob/4aee854da1a4afe367ea77802eb0ad60ea9047d2/supplement.pdf) or in the top level of the repository above.

### Videos

Here are three examples of moving pattern interfaces constructed in the article. First, an invasion of the pure conduction state by a hexagon pattern.

https://github.com/user-attachments/assets/6b7c4264-52ee-4328-9e12-e7b753f2ac8e

Second, an invasion of down-hexagons by up-hexagons.

https://github.com/user-attachments/assets/8d0aef92-c6e7-47c2-bd6f-268478be127a

Third, a two-stage invasion process consisting of a primary invasion of the pure conduction state by roll waves followed by a secondary invasion of the roll waves by a hexagon pattern.

https://github.com/user-attachments/assets/22c11ca5-411e-451a-810e-c1a4782e5c01

The videos are generated with Matlab using the code provided in the folder "videos". There, also the (compressed) videos displayed above can be found.

### Coefficients

The coefficients for the formal amplitude equations (and the reduced equations on the center manifolds) are calculated using the Mathematica Notebook found in the folder "Coefficient-Mathematica". The output of the notebook is also available as pdf and html there.

Please also note that the notebook can be read using the Wolfram Player, which is freely available [here](https://www.wolfram.com/player/?source=nav).

### Coefficient plots

The plots of the coefficients are generated using Python using the files contained in the folder "Coefficients-Plots". The python scripts read in the explicit expressions provided by the Mathematica Notebook, which are also made available in "Coefficients-Plots".

### Phase plane

The phase plane diagrams for the heteroclinic orbits on the reduced system are also generated the streamplot-function provided by Mathematica. The corresponding Mathematica Notebooks for both the square and hexagonal lattice are available in the folder "Phase planes". There, also a pdf version of the outputs of the notebooks can be found.
