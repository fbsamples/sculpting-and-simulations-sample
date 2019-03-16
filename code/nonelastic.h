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
// The following code is an affine transform (consisting
// of translation, rotation, and uniform scale) in
// differential equation form
///////////////////////////////////////////////////////

INLINE vec3
NonElasticEvaluateODE(float t, vec3 x, Deformation deformer)
{
    // advect the center of the move
    float originLerp = t - deformer.time;
    vec3 originAdvected = deformer.origin + deformer.linearVelocity*originLerp;

    vec3 R = x - originAdvected;
    vec3 trans = deformer.linearVelocity;
    vec3 disp = deformer.displacementGradientTensor * R;
    return trans + disp;
}

INLINE vec3
NonElasticEvaluateODE_TwoDeformers(float t, vec3 x, Deformation deformer0, Deformation deformer1)
{
    // advect the center of the move
    float originLerp = t - deformer0.time;
    vec3 originAdvected0 = deformer0.origin + deformer0.linearVelocity*originLerp;

	vec3 R0 = x - originAdvected0;
    vec3 trans0 = deformer0.linearVelocity;
    vec3 disp0 = deformer0.displacementGradientTensor * R0;
    vec3 u0 = trans0 + disp0;

    vec3 originAdvected1 = deformer1.origin + deformer1.linearVelocity*originLerp;
    vec3 R1 = x - originAdvected1;
    vec3 trans1 = deformer1.linearVelocity;
    vec3 disp1 = deformer1.displacementGradientTensor * R1;
    vec3 u1 = trans1 + disp1;

    return u0 + u1;
}

///////////////////////////////////////////////////////
// The following preprocessor code includes ODESolver multiple times
// to generate different variants of the solvers. This is necessary
// because GLSL doesn't support classes with functions, but we need
// to parameterize the solver with different evaluation functions.
///////////////////////////////////////////////////////

#define SCOPE(suffix) IntegrateNonElastic##suffix
#define EVALUATE NonElasticEvaluateODE
#define PARAMETERLIST Deformation deformer
#define PARAMETERS deformer
#include "ODESolvers.h"
#undef PARAMETERS
#undef PARAMETERLIST
#undef EVALUATE
#undef SCOPE

#define SCOPE(suffix) IntegrateNonElasticTwoDeformers##suffix
#define EVALUATE NonElasticEvaluateODE_TwoDeformers
#define PARAMETERLIST Deformation deformer0, Deformation deformer1
#define PARAMETERS deformer0, deformer1
#include "ODESolvers.h"
#undef PARAMETERS
#undef PARAMETERLIST
#undef EVALUATE
#undef SCOPE
