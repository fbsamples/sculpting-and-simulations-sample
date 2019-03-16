# sculpting-and-simulations-sample
is sample code that supplements the Sculpting and Simulations with 6DoF Controllers presentation at VRDC 2019.

## Examples
This code compiles in both C++ and GLSL. It is a header file only project, so it is recommended to wrap the #include in
a namespace in C++.

On the C++ side:
```
using namespace deformation
{
#include "deformation.h"
};

... each frame, fill out the Deformation::Pose struct with the current state of the 6DoF controller

deformation::Motion motion = buildMotion(previousFramePose, thisFramePose);
deformation::Deformation deformation = buildDeformation(motion);
deformation::Material material;
material.stiffness = 1.0f;
material.compressiblity = 0.5f;
deformation::Kelvinlet kelvinlet = buildKelvinlet(deformation, material, radius);

float maxerror = 0.001f;

for (uint i = 0; i < vertices.size(); i++)
{
    vertices[i] = deformation::deformKelvinlet(vertices[i], kelvinlet, material, stroke.innerRadius, maxerror);
}
```

The above code will deform the vertices on the CPU.

To enable GPU deformation, upload the Kelvinlet, Material, radius, and maxerror variables to the GPU. Then do this
per-vertex:
```
#include "deformation.h"

... in the vertex shader ...
    KelvinletParameters parameters;
    parameters.kelvinlet = kelvinlet;
    parameters.material = material;
    parameters.radius = radius;
    position = Kelvinlets_AdaptiveBS32(position, kelvinlet.time, kelvinlet.time + kelvinlet.dt, parameters, maxerror);
```

For more in-depth examples, look in test/test/test.cpp and in code/deformers.h.


## Requirements
* Windows
* Visual Studio 2017

## Building
There is a test executable that can be built with test/test.sln. Open that in Visual Studio 2017 and compile.

## How this sample code works
The code is a header file only project. It compiles in both C++ and GLSL.

To use in either C++ or GLSL, just `#include "deformation.h"`.

## Questions?

Email davidfarrell@oculus.com with any questions.

## License
sculpting-and-simulations-sample is BSD licensed, as found in the LICENSE file.
