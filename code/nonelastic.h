// Copyright(c) Facebook, Inc. and its affiliates.
// All rights reserved.
// 
// This source code is licensed under the BSD - style license found in the
// LICENSE file in the root directory of this source tree.

#ifdef __cplusplus
#pragma once
#endif 

///////////////////////////////////////////////////////
// NonElastic (non-elastic) deformation
// The following code doesn't use Kelvinlets
///////////////////////////////////////////////////////

struct NonElasticTwoDeformerParameters
{
    Deformation deformation[2];
};

INLINE vec3
NonElasticEvaluateODE(float t, vec3 x, Deformation d)
{
    // advect the center of the move from the start position to end position
    float originLerp = t - d.time;
    vec3 originAdvected = d.origin + d.linearVelocity*originLerp;

    vec3 R = x - originAdvected;
    vec3 trans = d.linearVelocity;
    vec3 disp = d.displacementGradientTensor * R;
    return trans + disp;
}

INLINE vec3
NonElasticEvaluateODE_TwoDeformers(float t, vec3 x, NonElasticTwoDeformerParameters p)
{
    // advect the center of the move from the start position to end position
    float originLerp = t - p.deformation[0].time;
    vec3 originAdvected0 = p.deformation[0].origin + p.deformation[0].linearVelocity*originLerp;
    vec3 R0 = x - originAdvected0;
    vec3 trans0 = p.deformation[0].linearVelocity;
    vec3 disp0 = p.deformation[0].displacementGradientTensor * R0;
    vec3 u0 = trans0 + disp0;

    vec3 originAdvected1 = p.deformation[1].origin + p.deformation[1].linearVelocity*originLerp;
    vec3 R1 = x - originAdvected1;
    vec3 trans1 = p.deformation[1].linearVelocity;
    vec3 disp1 = p.deformation[1].displacementGradientTensor * R1;
    vec3 u1 = trans1 + disp1;

    return u0 + u1;
}

///////////////////////////////////////////////////////
// The following preprocessor code includes ODESolver multiple times
// to generate different variants of the solvers. This is necessary
// because GLSL doesn't support classes with functions, but we need
// to parameterize the solver with different evaluation functions.
///////////////////////////////////////////////////////

#define SCOPE(suffix) NonElastic##suffix
#define evaluate NonElasticEvaluateODE
#define Parameters Deformation
#include "ODESolvers.h"
#undef Parameters
#undef evaluate
#undef SCOPE

#define SCOPE(suffix) NonElasticTwoDeformerParameters##suffix
#define evaluate NonElasticEvaluateODE_TwoDeformers
#define Parameters NonElasticTwoDeformerParameters
#include "ODESolvers.h"
#undef Parameters
#undef evaluate
#undef SCOPE
