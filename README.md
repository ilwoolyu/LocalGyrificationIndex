# Local Gyrification Index

## Description
The amount of cortical folding, or gyrification, is typically measured within local cortical regions covered by an equidistant geodesic or nearest neighborhood-ring kernel. However, without careful design, such a kernel can easily cover multiple sulcal and gyral regions that may not be functionally related. Furthermore, this can result in smoothing out details of cortical folding, which consequently blurs local gyrification measurements. In this paper, we propose a novel kernel shape to locally quantify cortical gyrification within sulcal and gyral regions. We adapt wavefront propagation to generate a spatially varying kernel shape that encodes cortical folding patterns: neighboring gyral crowns, sulcal fundi, and sulcal banks. For this purpose, we perform anisotropic wavefront propagation that runs fast along gyral crowns and sulcal fundi by solving a static Hamilton–Jacobi partial differential equation. The resulting kernel adaptively elongates along gyral crowns and sulcal fundi, while keeping a uniform shape over flat regions like sulcal banks. We then measure local gyrification within the proposed spatially varying kernel.

![image](https://user-images.githubusercontent.com/9325798/47728629-b6e60600-dc2c-11e8-8138-094841ae5c46.png)

## Installation
You can download and compile the source code using <a href="https://cmake.org/">CMake</a>. Or you can pull <a href="https://hub.docker.com/r/ilwoolyu/cmorph/">docker image</a>:
```
$ docker pull ilwoolyu/cmorph
```
## Usage
### Input
* surface file (.vtk): triangular 3D mesh
### Intermediate input and output
* sulcal and gyral curves (.scurve and .gcurve - or .bary): outputs of <a href="https://github.com/ilwoolyu/CurveExtraction">CurveExtraction</a> [[3]](#ref3)
* outer hull file (.vtk): output of <a href="https://github.com/ilwoolyu/klaplace">klaplace</a> [[4]](#ref4)
### Output
* lgi file (.txt): local gyrification index per vertex
### Commands
After build and install the required packages, type:
```
$ script/lgi -i input.vtk
```
If you have both pial and white surfaces, the results will be more accurate. Type:
```
$ script/lgi -i pial.vtk --white white.vtk
```
To change kernel size
```
$ script/lgi --kernel <area mm^2>
```
If you have a known reference population area, the kernel size will be automatically adjusted by the ratio between the surface area and the reference area.
```
$ script/lgi --ref <area mm^2>
```
For example, if kernel size=300 mm^2, reference area=150000 mm^2, input surface area=100000 mm^2, the kernel size is adjusted to 200 mm^2. In our work [[1](#ref1),[2](#ref2)], we used reference area of 166000 mm^2 (pial surface) and 77100 mm^2 (cerebral hull) with kernel size of 316 mm^2. To disable the kernel size adjustment, set --ref 0 or ignore this argument (default: 0).
>**Note 1**: Please use ~~--ref **166000** (pial surface area)~~ --ref **77100** (cerebral hull area) with --kernel **316** for consistent quantification independent of individual surfce areas unless you know a specific kernel size for each individual or plan to use an absolute kernel size.

>**Note 2**: The reference area refers to cerebral hull rather than pial surface.

In Docker, you need a sudo acces. To run local gyrification, type:
```
$ docker run \
         -v <LOCAL_INPUT_PATH>:/INPUT/ \
         --rm ilwoolyu/cmorph:1.0 \
         lgi -i /INPUT/input.vtk
```
## Implementation Details
### Sulcal/gyral curve extraction
The following command line will generate "output.scurve", "output.gcurve", and ".bary":
```
$ CurveExtraction -i input.vtk -o output --sulcus --gyrus --bary --noVtk
```
Or if both pial and white surfaces are available, the following commands give better extraction results:
```
$ CurveExtraction -i pial.vtk -o output --sulcus --bary --noVtk
$ CurveExtraction -i white.vtk -o output --gyrus --bary --noVtk
```
See [CurveExtraction (sulcal/gyral curves)](https://github.com/ilwoolyu/CurveExtraction) for more options.
### Outer hull creation
An initial outer hull surface needs to be generated for outer hull creation. The cortical surface is voxelized and the morphological operation is applied on it.

To create a binary volume image of the input surface using <a href="https://www.mathworks.com/matlabcentral/fileexchange/27390-mesh-voxelisation">Mesh voxelisation</a>. From the volume, the outer hull can be obtained using <a href="https://www.mathworks.com/help/matlab/ref/isosurface.html">isosurface</a>. Both are implemented in MATLAB.
```
$ matlab OuterHull('input.vtk', 'outer_hull.vtk');
```
### Outer hull correspondence
To find a Laplacian shape correspondence, the following command will give Laplacian trajectories using [klaplace (outer hull correspondence)](https://github.com/ilwoolyu/klaplace):
```
$ klaplace -dims 128 input.vtk outer_hull.vtk -surfaceCorrespondence outer_hull
```
The trjectories will be generated in "outer_hull_warpedMesh.vtp".
Let's trace the final destinations of the trajectories to obtain the outer hull:
```
$ klaplace -conv outer_hull_warpedMesh.vtp outer_hull_corr.vtk
```
> **Note 1**: It would be useful if do some smoothing on "outer_hull_corr.vtk" since it's very rough mesh since isosurface does not provide smooth mesh.

> **Note 2**: ~~Since the outputs of klaplace consume a huge disk space, it is recommended to delete all but "outer_hull_corr.vtk".~~ This issue has been fixed.
### Local gyrification index
The following command line gives local gyrification index per vertex in "output.lgi.map.316.txt":
```
$ Gyrification \
               -i input.vtk \
               -o output \
               --outer outer_hull_corr.vtk \
               -s output.scurve.bary \
               -g output.gcurve.bary \
               -m 316 \
               --speed 0.2
```
Barycentric curves can provide dense points along sulcal/gyral regions.

To enable multi-thread support (OpenMP):
```
$ HSD --nThreads <# of threads>
```
More technical details (theory, parameter choice, etc.) can be found in [[1](#ref1),[2](#ref2)].
> **Note 1**: If a population area is known, --refHullArea [area] will adjust the area size of "-m" with respect to the input surface area. The use of --refHullArea is *recommended* particularly for neurodevelopmental studies.

> **Note 2**: -t [area] will create different lgi measurements in a given interval of area; e.g., -t 100 -m 300 will give lgi at area of 100, 200, and 300 mm^2.

## Dependency
* [MeshLib (general mesh processing)](https://github.com/ilwoolyu/MeshLib)
* [SlicerExecutionModel (CLI)](https://github.com/Slicer/SlicerExecutionModel)

## Required Components
* [CurveExtraction (sulcal/gyral curves)](https://github.com/ilwoolyu/CurveExtraction)
* [klaplace (outer hull correspondence)](https://github.com/ilwoolyu/klaplace)
* [MATLAB (initial outer hull creation)](https://www.mathworks.com/products/matlab.html)

## References
<ol>
<li><a id="ref1"></a>Lyu, I., Kim, S., Girault, J., Gilmore, J., Styner, M., <a href="https://doi.org/10.1016/j.media.2018.06.009">A Cortical Shape-Adaptive Approach to Local Gyrification Index</a>, <i>Medical Image Analysis</i>, 48, 244-258, 2018
<li><a id="ref2"></a>Lyu, I., Kim, S., Bullins, J., Gilmore, J., Styner, M., <a href="http://dx.doi.org/10.1007/978-3-319-66182-7_4">Novel Local Shape-Adaptive Gyrification Index with Application to Brain Development</a>, <i>Medical Image Computing and Computer Assisted Intervention (MICCAI) 2017</i>, LNCS10433, 31-39, 2017
<li><a id="ref3"></a>Lyu, I., Kim, S., Woodward, N., Styner, M., Landman, B., <a href="http://dx.doi.org/10.1109/TMI.2017.2787589">TRACE: A Topological Graph Representation for Automatic Sulcal Curve Extraction</a>, <i>IEEE Transactions on Medical Imaging</i>, 37(7), 1653-1663, 2018</li>
<li><a id="ref4"></a>Lee, J., Kim, S., Oguz, I., Styner, M., <a href="http://dx.doi.org/10.1117/12.2216420">Enhanced Cortical Thickness Measurements for Rodent Brains via Lagrangian-based RK4 Streamline Computation</a>, <i>SPIE Medical Imaging 2016</i>, SPIE9784, 97840B-1-97840B-10, 2016</li>
</ol>
