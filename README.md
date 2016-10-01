#kirkpatrick-point-location
This repo contains an implementation of the Kirkpatrick's point location algorithm. Given a set of *N* points, the algorithm preprocesses and triangulates the set in **O(N)** time. Henceforth, it supports locating query points *q* (i.e. finding which triangle from the triangulation *q* is contained in) in **O(log N)** time. 

##Implementation
The bulk of preprocessing required by Kirkpatrick's algorithm takes place in the constructor of the **Kirkpatrick** class (*kirkpatrick.py*), which creates the DAG that supports location queries. We store the initial triangulation of the bounding triangle and the *N* points in graph data structure that we created. The data structure, a class called **MyGraph** in *MyGraph.py*, is a modified adjacency matrix. At each possible edge [i, j], instead of storing a boolean indicating the existence of that edge, we store up to 2 faces that that edge is a part of. A face is represented by an instance of class **Piece**, which is also in *MyGraph.py*. If the edge doesn’t exist, then no faces are stored at [i, j]. Our graph structure supports most operations in constant time given the nature of the Kirkpatrick's algorithm. To compute the degree/neighbors of a vertex does require **O(N)** time however.

The rest of our implementation is fairly standard. For triangulation of the initial *N* points, the triangulation including the bounding triangle, as well as the triangulation of the holes during the DAG creation, we use the *SciPy* library’s triangulation method and the *Tri* library’s constrained triangulation method. Both methods perform Delaunay triangulation; for example, the *Tri* library uses triangle flips (if the Delaunay criterion don’t hold) to construct it’s a triangulation.

## Examples
To gain a better appreciation of how the algorithm works, it can be useful to view these [examples](https://github.com/brianonpig/kirkpatrick-point-location/blob/master/examples.pdf), which provides step-by-step illustrations of preprocessing and point locations. There are also some bonus timing benchmarks!

##Acknowledgements
In addition to using a handful of common Python libraries like *numpy*, *SciPy*, and *matplotlib*, the code also leverages the *tri* library for Delauney triangulation. Instructions for setup can be found [here](https://pypi.python.org/pypi/tri/0.3).