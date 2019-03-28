# sculpting-and-simulations-sample

<p align="center">
  <img src="https://github.com/fbsamples/sculpting-and-simulations-sample/releases/latest/download/deformation.gif"/>
</p>

This repo contains sample code that supplements the talk "Sculpting and Simulations with 6DoF Controllers" presented at VRDC 2019.

The slides from the talk (with video + speaker notes) can be downloaded at [SculptingAndSimulations.pptx](https://github.com/fbsamples/sculpting-and-simulations-sample/releases/latest/download/SculptingAndSimulations.pptx)

This sample code shows how to deform meshes using regularized Kelvinlets, as well as using affine transformations. It also includes several flavors of explicit ODE solvers. This code originated in Oculus Medium.

The code compiles in both C++ and GLSL, consisting of only header files.

There is a small test project in the `test` directory. This deforms a mesh eight different ways both as a test and also as an example of how to use this code.

For additional notes on the ODE solvers, see the included document [NotesOnODESolvers.pdf](NotesOnODESolvers.pdf?raw=true).

For information on how the Kelvinlets equations were calibrated, see the included document [KelvinletsCalibration.pdf](KelvinletsCalibration.pdf?raw=true).


## Example

First, you need to `#include "deformation.h"`. In C++, it's recommended to wrap the `#include` in a namespace:

```
using namespace deformation
{
    #include "../code/deformation.h"
};
```

(if this results in many `identifier not found` compiler errors about `sqrt`, `atan2`, `pow`, etc., then first `#include <cmath>`).

Next, you need to create a `deformation::Pose` object once each frame from the 6DoF controller's state:

```
deformation::Pose pose;
pose.translation = ovrPose.translation;
pose.orientation = ovrPose.orientation;
pose.scale = 1.0f;
pose.time = gApp.getElapsedTime();
```

Then, you build the deformation/Kelvinlets data once per frame on the CPU:

```
deformation::Motion motion = buildMotion(previousFramePose, thisFramePose);
deformation::Deformation deformation = buildDeformation(motion);
float stiffness = 1.0f;
float compressibility = 0.5f;
deformation::Kelvinlet kelvinlet = buildKelvinlet(deformation, stiffness, compressibility, radius);
```

Finally, you can deform the mesh. To deform the mesh each frame on the CPU, loop over your vertices and do a single RK4 integration:

```
for (uint i = 0; i < vertices.size(); i++)
{
    vertices[i] = deformation::IntegrateNonElastic_RungeKutta(vertices[i], kelvinlet.time, kelvinlet.time+kelvinlet.dt, kelvinlet);
}
```

Alternatively, to deform the mesh on the GPU, upload the Kelvinlet struct to the GPU and do a single RK4 integration:

```
#include "deformation.h"

vertexpos = IntegrateNonElastic_RungeKutta(vertexpos, kelvinlet.time, kelvinlet.time+kelvinlet.dt, kelvinlet);
```

For more in-depth examples, look in `test/test.cpp`. That test loads a mesh from disk, deforms it eight different ways, and writes each deformed mesh as an .OBJ file to `test/data/testresultN.obj`.

## Requirements
* This code has only been verified on Windows using Visual Studio 2017, but should run anywhere.

## Building
There is a test project that can be built with test/test.sln. Open that in Visual Studio 2017 and compile.

## How this sample code works
There are simple structures in `deformation.h` that contain the data needed for mesh deformation: `Pose`, `Motion`, `Deformation`, and `Kelvinlet`.

To use continuous deformation, where every frame the mesh is deformed by movement from the previous frame to now, you should use a single R4 integration per frame. Use one of these functions for that:

```
IntegrateNonElastic_RungeKutta(vec3 position, float tstart, float tend, Deformation deformation);
IntegrateKelvinlet_RungeKutta(vec3 position, float tstart, float tend, Kelvinlet kelvinlet);
```

For Medium move tool style deformation, where there are two keyframes, and the user freely move around the second keyframe, a single RK4 integration has far too much error. Instead, use one of the adaptive integrators. Use one of these functions for that:

```
IntegrateNonElastic_AdaptiveRK(vec3 position, float tstart, float tend, float maxerror, Deformation deformation);
IntegrateNonElastic_AdaptiveRKF45(vec3 position, float tstart, float tend, float maxerror, Deformation deformation);
IntegrateNonElastic_AdaptiveDP54(vec3 position, float tstart, float tend, float maxerror, Deformation deformation);
IntegrateNonElastic_AdaptiveBS32(vec3 position, float tstart, float tend, float maxerror, Deformation deformation);
IntegrateKelvinlet_AdaptiveRK(vec3 position, float tstart, float tend, float maxerror, Kelvinlet kelvinlet);
IntegrateKelvinlet_AdaptiveRKF45(vec3 position, float tstart, float tend, float maxerror, Kelvinlet kelvinlet);
IntegrateKelvinlet_AdaptiveDP54(vec3 position, float tstart, float tend, float maxerror, Kelvinlet kelvinlet);
IntegrateKelvinlet_AdaptiveBS32(vec3 position, float tstart, float tend, float maxerror, Kelvinlet kelvinlet);
```

There are also functions to blend two deformers. For continuous deformation using two deformers, use one of these functions:

```
IntegrateNonElasticTwoDeformers_RungeKutta(vec3 position, float tstart, float tend, Deformation deformation);
IntegrateKelvinletTwoDeformers_RungeKutta(vec3 position, float tstart, float tend, Kelvinlet kelvinlet);
```

For Medium move tool style deformation with two blended deformers, use one of these functions:
```
IntegrateNonElasticTwoDeformers_AdaptiveRK(vec3 position, float tstart, float tend, float maxerror, Deformation deformation0, Deformation deformation1);
IntegrateNonElasticTwoDeformers_AdaptiveRKF45(vec3 position, float tstart, float tend, float maxerror, Deformation deformation0, Deformation deformation1);
IntegrateNonElasticTwoDeformers_AdaptiveDP54(vec3 position, float tstart, float tend, float maxerror, Deformation deformation0, Deformation deformation1);
IntegrateNonElasticTwoDeformers_AdaptiveBS32(vec3 position, float tstart, float tend, float maxerror, Deformation deformation0, Deformation deformation1);
IntegrateKelvinletTwoDeformers_AdaptiveRK(vec3 position, float tstart, float tend, float maxerror, Kelvinlet kelvinlet0, Kelvinlet kelvinlet1);
IntegrateKelvinletTwoDeformers_AdaptiveRKF45(vec3 position, float tstart, float tend, float maxerror, Kelvinlet kelvinlet0, Kelvinlet kelvinlet1);
IntegrateKelvinletTwoDeformers_AdaptiveDP54(vec3 position, float tstart, float tend, float maxerror, Kelvinlet kelvinlet0, Kelvinlet kelvinlet1);
IntegrateKelvinletTwoDeformers_AdaptiveBS32(vec3 position, float tstart, float tend, float maxerror, Kelvinlet kelvinlet0, Kelvinlet kelvinlet1);
```

The different flavors of the `Adaptive*` functions have different tradeoffs in terms of performance. Medium uses AdaptiveBS32.

`maxerror` is very application specific. It is generally a good idea to set it to some small world space value. In
Medium, which uses units of meters, `maxerror` is set to 0.00013f, but that is scaled as you scale your sculpt up and
down. Larger values of maxerror are faster for the adaptive algorithms to compute, but return less accurate answers.

## Questions?

Email davidfarrell@oculus.com with any questions.

## License
sculpting-and-simulations-sample is BSD licensed, as found in the LICENSE file.
