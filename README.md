# cr-case
Creation of Python-based interactive charts on various topics of fundamental structural analysis. The concepts covered are some of those in the course AE2135-I Structural Analysis & Design of the BSc Aerospace Engineering at the Technical University of Delft.

The charts are intended to show interactively some of the ideas introduced in the course, serving as comparison between exact and approximate approaches and visualizing the influence of the shape in the geometrical properties and the elastic analysis. It is possible to vary geometries and shapes as well as to reach out the information transmitted by the chart in an intuitive manner.

The main Python libraries that are used in the interactive tool creation are Plotly (https://plot.ly) and Dash (https://plot.ly/products/dash/) which allow for the implementation of sliders, dropdown menus and fill-in boxes.

## Shear Flow Discretization [flow]

This chart shows the shear flow obtained from a static ideal analysis in various cross sections (I-beam, T-beam, C-beam, U-beam, rectangular and open circular) side-by-side with the shear flow obtained from a lumped area discretization. It allows to visualize the effect of discretizing the cross section in 1-dimensional elements as far as the shear flow is concerned. Variable cross-section dimensions may be introduced as well as various levels of discretization refinement.

The shear flow is the most common way to conceptualize the mechanical response in a thin-walled section when a shear force is applied. The shear flow is obtained as the shear stress per given thickness which accordingly has the dimensions of force per unit of length. For a given structure the shear load is considered to be applied at the shear center, i.e. the point in space at which the shear force causes no torsional deformation of the section.

The shear flow through a section of thickness t is calculated using the following equation ([1],[2]):

![continuous equation](./resources/eq_cont.png)

where Vx and Vy are the shear forces perpendicular to the neutral axisx; Ix, Iy and Ixy are the second moments of area about the neutral axis; and x and y are the distances from the section coordinate to the centroid of the cross-section.

When discretizing the section a discretized version of the shear flow equation is used:


where Br are the equivalent contributing areas corresponding to each of the discretization points. These boom areas are calculated upon the force resultant equilibrium of the forces acting at each discretization point in the direction of interest:



## Usage

The apps are created using Dash which allows them to be run locally in "development" mode by default. However to share the Dash app it is necessary to "deploy" it to a server. In order to deploy the app to a server like Heroku the following instructions may be followed (https://dash.plot.ly/deployment).

## Development Environment

This app has been developed using Python 3.7.3.

## Bibliography

[1] Megson, T. H. G. 2013. Aircraft Structures for Engineering Students (version Fifth edition.). Fifth. Elsevier Aerospace Engineering Series. Oxford: Butterworth-Heinemann

[2] Aerospace Mechanics and Materials. Tu Delft OpenCourseWare. 7/12/2019. (https://ocw.tudelft.nl/courses/aerospace-mechanics-of-materials/)

