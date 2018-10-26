# LocalGyrificationIndex

## Introduction
The amount of cortical folding, or gyrification, is typically measured within local cortical regions covered by an equidistant geodesic or nearest neighborhood-ring kernel. However, without careful design, such a kernel can easily cover multiple sulcal and gyral regions that may not be functionally related. Furthermore, this can result in smoothing out details of cortical folding, which consequently blurs local gyrification measurements. Here, we propose a novel kernel shape to locally quantify cortical gyrification within sulcal and gyral regions.
### Input
* surface file (.vtk): triangular 3D mesh
* sulcal and gyral curves (.scurve and .gcurve - or .bary): outputs of <a href="https://github.com/ilwoolyu/CurveExtraction">CurveExtraction</a> [3]
* outer hull file (.vtk): output of <a href="https://github.com/ilwoolyu/klaplace">klaplace</a> [4]
### Output
* lgi file (.txt): local gyrification index per vertex
### Usage
After build and install the required packages, type:
```
script/lgi --i input.vtk
```
### Details
#### Sulcal/gyral curve extraction
The following command line will generate "output.scurve", "output.gcurve", and ".bary":<br />
```
CurveExtraction -i input.vtk -o output --sulcus --gyrus --bary --noVtk
```
Or if both pial and white surfaces are available, the following commands give better extraction results:<br />
```
CurveExtraction -i pial.vtk -o output --sulcus --bary --noVtk
CurveExtraction -i white.vtk -o output --gyrus --bary --noVtk
```
#### Outer hull creation
To create outer hull, an initial outer hull surface needs to be generated. The cortical surface is voxelized and the morphological operation is applied on it.<br />
To create a binary volume image of the input surface using <a href="https://www.mathworks.com/matlabcentral/fileexchange/27390-mesh-voxelisation">Mesh voxelisation</a>. From the volume, the outer hull can be obtained using <a href="https://www.mathworks.com/help/matlab/ref/isosurface.html">isosurface</a>. Both are implemented in MATLAB.<br />
```
OuterHull('input.vtk', 'outer_hull.vtk');
```
#### Outer hull correspondence
To find a Laplacian shape correspondence, the following command will give Laplacian trajectories:<br />
```
klaplace -dims 256 input.vtk outer_hull.vtk -surfaceCorrespondence outer_hull
```
The trjectories will be generated in "outer_hull_warpedMesh.vtp".
Let's trace the final destinations of the trajectories to obtain the outer hull:<br />
```
klaplace -conv outer_hull_warpedMesh.vtp outer_hull_corr.vtk
```
* Note 1: It would be useful if do some smoothing on "outer_hull_corr.vtk" since it's very rough mesh since isosurface does not provide smooth mesh.
* Note 2: Since the outputs of klaplace consume a huge disk space, it is recommended to delete all but "outer_hull_corr.vtk".<br />
#### Local gyrification index
The following command line gives local gyrification index per vertex in "output.lgi.map.316.txt":
```
Gyrification -i input.vtk -o output --outer outer_hull_corr.vtk -s output.scurve.bary -g output.gcurve.bary -m 316 --speed 0.2
```
Note barycentric curves can provide dense points along sulcal/gyral regions.
More technical details (theory, parameter choice, etc.) can be found in [1,2].<br />
* Note 1: If a population area is known, --poulationArea [area] will adjust the area size of "-m" with respect to the input surface area.
* Note 2: -t [area] will create different lgi measurements in a given interval of area; e.g., -t 100 -m 300 will give lgi at area of 100, 200, and 300 mm^2.
## Dependency
* <a href="https://github.com/ilwoolyu/MeshLib">MeshLib (general mesh processing)</a><br />
* <a href="https://github.com/Slicer/SlicerExecutionModel">SlicerExecutionModel (CLI)</a>

## Required Components
* <a href="https://github.com/ilwoolyu/CurveExtraction">CurveExtraction (sulcal/gyral curves)</a>
* <a href="https://github.com/ilwoolyu/klaplace">klaplace (outer hull correspondence)</a>
* <a href="https://www.mathworks.com/products/matlab.html">MATLAB (initial outer hull creation)</a>

## References
<ol>
<li>Lyu, I., Kim, S., Girault, J., Gilmore, J., Styner, M., <a href="https://doi.org/10.1016/j.media.2018.06.009">A Cortical Shape-Adaptive Approach to Local Gyrification Index</a>, <i>Medical Image Analysis</i>, 48, 244-258, 2018
<li>Lyu, I., Kim, S., Bullins, J., Gilmore, J., Styner, M., <a href="http://dx.doi.org/10.1007/978-3-319-66182-7_4">Novel Local Shape-Adaptive Gyrification Index with Application to Brain Development</a>, <i>Medical Image Computing and Computer Assisted Intervention (MICCAI) 2017</i>, LNCS10433, 31-39, 2017
<li>Lyu, I., Kim, S., Woodward, N., Styner, M., Landman, B., <a href="http://dx.doi.org/10.1109/TMI.2017.2787589">TRACE: A Topological Graph Representation for Automatic Sulcal Curve Extraction</a>, <i>IEEE Transactions on Medical Imaging</i>, 37(7), 1653-1663, 2018</li>
<li>Lee, J., Kim, S., Oguz, I., Styner, M., <a href="http://dx.doi.org/10.1117/12.2216420">Enhanced Cortical Thickness Measurements for Rodent Brains via Lagrangian-based RK4 Streamline Computation</a>, <i>SPIE Medical Imaging 2016</i>, SPIE9784, 97840B-1-97840B-10, 2016</li>
</ol>
