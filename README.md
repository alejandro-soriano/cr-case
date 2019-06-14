# cr_case
Creation of Python-based interactive charts on various topics of fundamental structural analysis. The concepts covered are some of those in the course AE2135-I Structural Analysis & Design of the BSc Aerospace Engineering at the Technical University of Delft.

The charts are intended to show interactively some of the ideas introduced in the course, serving as comparison between exact and approximate approaches and visualizing the influence of the shape in the geometrical properties and the elastic analysis. It is possible to vary geometries and shapes as well as to reach out the information transmitted by the chart in an intuitive manner.

The main Python library that is used in the interactive tool creation is Plotly (plot.ly) which allows for the implementation of sliders, selection boxes and fill boxes.

The current created charts are:

- Shear flow discretization

This chart shows the shear flow obtained from a static ideal analysis in various cross sections (I-beam, C-beam, T-beam,...) side-by-side with the shear flow obtained from a lumped area discretization. It allows to visualize the effect of discretizing the cross section in 1-dimensional elements as far as the shear flow is concerned. Variable cross-section dimensions may be introduced as well as various levels of discretization refinement.
