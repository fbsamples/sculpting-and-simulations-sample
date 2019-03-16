// Copyright(c) Facebook, Inc. and its affiliates.
// All rights reserved.
// 
// This source code is licensed under the BSD - style license found in the
// LICENSE file in the root directory of this source tree.

#ifndef PI
#define PI 3.14159265f
#endif

///////////////////////////////////////////////////////
// Kelvinlets
///////////////////////////////////////////////////////

// Set this to 1 to use biscale Kelvinlets (for tighter falloff)
// Otherwise, set to 0 to use a single Kelvinlet
// See section 5 of Kelvinlets paper
#define BISCALE_FALLOFF 1

// The radius of the second biscale Kelvinlet is the first radius times this factor
// 1.1 is the same value in the Kelvinlets paper
#define BISCALE_RADIUS 1.1f

// The default stiffness. The calibration functions use this value so that 
// a Kelvinlet with a stiffness value of KELVINLETS_DEFAULT_STIFFNESS will
// cause the sculpt's deformation at the tool's tip to exactly follow the 
// tool's movement. The Kelvinlets evaluation function specifies the 
// material's stiffness value, which can be the same (to exactly match the 
// movement) or different (to simulate stiffness in the material).
// When evaluating the displacement vector, stiffness values between 1.0
// and 10.0 are useful. Stiffness values less than KELVINLETS_DEFAULT_STIFFNESS 
// cause a small movement in the Touch controller to generate a large deformation. 
// Values greater than KELVINLETS_DEFAULT_STIFFNESS cause the material to
// be more stiff.
// The actual default stiffness value doesn't matter, as the calibration function
// cancels it out, so we use 1.0 as a reasonable default.
// see equation 7 in Kelvinlets paper (parameterizing the brush tip displacement).
#define KELVINLETS_DEFAULT_STIFFNESS 1.0f

// Calibration factor functions
// These functions return a value that should be multiplied
// with the force vector and matrix before they are evaluated.
// The calibration is set up so that a stiffness value of
// KELVINLETS_DEFAULT_STIFFNESS causes a point at the move
// tool's tip to move exactly with the move tool.
// The calibration is linear in terms of the force vector/matrix,
// so this can be done on the CPU once a frame.
// See equation 7 in Kelvinlets paper for translation, 
// and equations 4/5/6 in the Kelvinlets supplemental material
// for twist/scale/pinch.

INLINE float
KTranslationCalibrationFactor(float radius, float compressibility)
{
    float a = 1 / (4 * PI * KELVINLETS_DEFAULT_STIFFNESS);
    float b = a / (4 * (1 - compressibility));
    float c = 2 / (3 * a - 2 * b);

#if BISCALE_FALLOFF
    float radius1 = radius * BISCALE_RADIUS;
    float biscale_falloff = 1 / (1 / radius - 1 / radius1);
    return c * biscale_falloff;
#else
    return c * radius;
#endif
}

INLINE float
KTwistCalibrationFactor(float radius, float compressibility)
{
    float a = 1 / (4 * PI * KELVINLETS_DEFAULT_STIFFNESS);
    float c = -2 / (5 * a);

#if BISCALE_FALLOFF
    float radius1 = radius * BISCALE_RADIUS;
    float biscale_falloff = 1 / (1 / (radius*radius*radius) - 1 / (radius1*radius1*radius1));
    return c * biscale_falloff;
#else
    float singlescale_falloff = radius * radius * radius;
    return c * singlescale_falloff;
#endif
}

INLINE float
KScaleCalibrationFactor(float radius, float compressibility)
{
    float a = 1 / (4 * PI * KELVINLETS_DEFAULT_STIFFNESS);
    float b = a / (4 * (1 - compressibility));
    float c = 2 / ((2 * b - a) * 5);

#if BISCALE_FALLOFF
    float radius1 = radius * BISCALE_RADIUS;
    float biscale_falloff = 1 / (1 / (radius*radius*radius) - 1 / (radius1*radius1*radius1));
    return c * biscale_falloff;
#else
    return c * singlescale_falloff;
#endif
}

INLINE float
KPinchCalibrationFactor(float radius, float compressibility)
{
    float a = 1 / (4 * PI * KELVINLETS_DEFAULT_STIFFNESS);
    float b = a / (4 * (1 - compressibility));
    float c = 2 * radius*radius*radius / (4 * b - 5 * a);

#if BISCALE_FALLOFF
    float radius1 = radius * BISCALE_RADIUS;
    float biscale_falloff = 1 / (1 / (radius*radius*radius) - 1 / (radius1*radius1*radius1));
    return c * biscale_falloff;
#else
    float singlescale_falloff = radius * radius * radius;
    return c * singlescale_falloff;
#endif
}

// Kelvinlets evaluation
struct KelvinletParameters
{
    Kelvinlet kelvinlet;
    Material material;
    float radius;
};

struct KelvinletTwoDeformerParameters
{
    Kelvinlet kelvinlet[2];
    Material material;
    float radius;
};

// Returns a displacement vector for the point R for a single Kelvinlet
// You probably want to call KTranslation(), which handles either single
// or biscale Kelvinlets based on the BISCALE_FALLOFF macro
// equation 6 in Kelvinlets paper
INLINE vec3 
KTranslationInner(vec3 R, vec3 loadForceVector, float radius, float stiffness, float compressibility)
{
    float a = 1 / (4 * PI * stiffness);
    float b = a / (4 * (1 - compressibility));

    float rnorm = length(R);
    float re = sqrt(rnorm*rnorm + radius * radius);

    float firstterm = (a - b) / re;

    vec3 Rsecondterm = b / (re*re*re) * R;
    mat3x3 secondterm = {
        vec3(R.x*Rsecondterm.x, R.x*Rsecondterm.y, R.x*Rsecondterm.z),
        vec3(R.y*Rsecondterm.x, R.y*Rsecondterm.y, R.y*Rsecondterm.z),
        vec3(R.z*Rsecondterm.x, R.z*Rsecondterm.y, R.z*Rsecondterm.z)
    };

    float thirdterm = (a * radius*radius) / (2 * re*re*re);

    vec3 displacement = firstterm*loadForceVector + secondterm*loadForceVector + thirdterm*loadForceVector;

    return displacement;
}

// Returns a displacement vector for the point R due to translation, using either
// a single or biscale Kelvinlet
// section 5 of Kelvinlets paper
INLINE vec3 
KTranslation(vec3 R, vec3 loadForceVector, float radius, float stiffness, float compressibility)
{
#if BISCALE_FALLOFF
    return KTranslationInner(R, loadForceVector, radius, stiffness, compressibility) - KTranslationInner(R, loadForceVector, radius*BISCALE_RADIUS, stiffness, compressibility);
#else
    return KTranslationInner(R, loadForceVector, radius, stiffness, compressibility);
#endif
}

// Returns the displacement vector for a skew-symmetric 3x3 rotation matrix
// You probably want to call KTwist(), which handles either single
// or biscale Kelvinlets based on the BISCALE_FALLOFF macro
// section 6 (extension to affine loads/twisting) of Kelvinlets paper
INLINE vec3
KTwistInner(vec3 R, mat3x3 loadForceMatrix, float radius, float stiffness, float compressibility)
{
    float a = 1 / (4 * PI * stiffness);

    float Rnorm = length(R);
    float re = sqrt(Rnorm*Rnorm + radius * radius);

    return -a * (1 / (re*re*re) + 3 * radius*radius / (2 * re*re*re*re*re)) * loadForceMatrix * R;
}

// Returns a displacement vector for a skew-symmetric 3x3 rotation matrix
// using either a single or biscale Kelvinlet
// section 5 of Kelvinlets paper
INLINE vec3
KTwist(vec3 R, mat3x3 loadForceMatrix, float radius, float stiffness, float compressibility)
{
#if BISCALE_FALLOFF
    return KTwistInner(R, loadForceMatrix, radius, stiffness, compressibility) - KTwistInner(R, loadForceMatrix, radius*BISCALE_RADIUS, stiffness, compressibility);
#else
    return KTwistInner(R, loadForceMatrix, radius, stiffness, compressibility);
#endif
}

// Returns the displacement vector for a uniform scale matrix
// You probably want to call KScale(), which handles either single
// or biscale Kelvinlets based on the BISCALE_FALLOFF macro
// section 6 (extension to affine loads/scaling) of Kelvinlets paper
INLINE vec3
KScaleInner(vec3 R, mat3x3 loadForceMatrix, float radius, float stiffness, float compressibility)
{
    float a = 1 / (4 * PI * stiffness);
    float b = a / (4 * (1 - compressibility));

    float Rnorm = length(R);
    float re = sqrt(Rnorm*Rnorm + radius * radius);

    return (2 * b - a) * (1 / (re*re*re) + 3 * radius*radius / (2 * re*re*re*re*re)) * loadForceMatrix * R;
}

// Returns a displacement vector for a uniform scale matrix (of the form Is where I=identity matrix, s=scalar)
// using either a single or biscale Kelvinlet
// section 5 of Kelvinlets paper
INLINE vec3
KScale(vec3 R, mat3x3 loadForceMatrix, float radius, float stiffness, float compressibility)
{
#if BISCALE_FALLOFF
    return KScaleInner(R, loadForceMatrix, radius, stiffness, compressibility) - KScaleInner(R, loadForceMatrix, radius*BISCALE_RADIUS, stiffness, compressibility);
#else
    return KScaleInner(R, loadForceMatrix, radius, stiffness, compressibility);
#endif
}

// Returns the displacement vector for a symmetric and traceless 3x3 matrix
// that represents pinch/bulge
// You probably want to call KPinch(), which handles either single
// or biscale Kelvinlets based on the BISCALE_FALLOFF macro
// section 6 (extension to affine loads/pinching) of the Kelvinlets paper
INLINE vec3
KPinchInner(vec3 R, mat3x3 loadForceMatrix, float radius, float stiffness, float compressibility)
{
    float a = 1 / (4 * PI * stiffness);
    float b = a / (4 * (1 - compressibility));

    float Rnorm = length(R);
    float re = sqrt(Rnorm*Rnorm + radius * radius);

    mat3x3 I = {
        vec3(1, 0, 0),
        vec3(0, 1, 0),
        vec3(0, 0, 1)
    };

    vec3 firstterm = (2 * b - a) / (re*re*re) * loadForceMatrix * R;

    vec3 secondterm = -3 / (2 * re*re*re*re*re) * (2 * b*dot(R, loadForceMatrix*R)*I + a * radius*radius*loadForceMatrix) * R;

    return firstterm + secondterm;
}

// Returns a displacement vector for a pinch matrix
// using either a single or biscale Kelvinlet
// section 5 of Kelvinlets paper
INLINE vec3
KPinch(vec3 R, mat3x3 loadForceMatrix, float radius, float stiffness, float compressibility)
{
#if BISCALE_FALLOFF
    return KPinchInner(R, loadForceMatrix, radius, stiffness, compressibility) - KPinchInner(R, loadForceMatrix, radius*BISCALE_RADIUS, stiffness, compressibility);
#else
    return KPinchInner(R, loadForceMatrix, radius, stiffness, compressibility);
#endif
}

// Returns the displacement vector for translation/rotation/scale of the position x at time t
// This advects the load origin so that a particle that starts at 
// loadOrigin is moved to loadEnd
INLINE vec3
KEvaluate(float t, vec3 x, KelvinletParameters p)
{
    // advect the center of the Kelvinlet from the start position to end position
    float originLerp = t - p.kelvinlet.time;
    vec3 loadOriginAdvected = p.kelvinlet.origin + p.kelvinlet.linearVelocity * originLerp;
    vec3 R = x - loadOriginAdvected;
    vec3 Kt = KTranslation(R, p.kelvinlet.forceVector, p.radius, p.material.stiffness, p.material.compressibility);
    vec3 Kr = KTwist(R, p.kelvinlet.twistForceMatrix, p.radius, p.material.stiffness, p.material.compressibility);
    vec3 Ks = KScale(R, p.kelvinlet.scaleForceMatrix, p.radius, p.material.stiffness, 0.0f);    // compressibility of 0.5 causes divide by 0 in Kelvinlets math, so for scale, we just set compressibility to 0.0f
    return Kt + Kr + Ks;
}

// Same as above, but with two deformers
INLINE vec3
KEvaluateTwoDeformers(float t, vec3 x, KelvinletTwoDeformerParameters p)
{
    // pose 0
    // advect the center of the Kelvinlet from the start position to end position
    float originLerp = t - p.kelvinlet[0].time;
    vec3 loadOriginAdvected0 = p.kelvinlet[0].origin + p.kelvinlet[0].linearVelocity * originLerp;
    vec3 R0 = x - loadOriginAdvected0;
    vec3 Kt0 = KTranslation(R0, p.kelvinlet[0].forceVector, p.radius, p.material.stiffness, p.material.compressibility);
    vec3 Kr0 = KTwist(R0, p.kelvinlet[0].twistForceMatrix, p.radius, p.material.stiffness, p.material.compressibility);
    vec3 Ks0 = KScale(R0, p.kelvinlet[0].scaleForceMatrix, p.radius, p.material.stiffness, 0.0f);    // compressibility of 0.5 causes divide by 0 in Kelvinlets math, so for scale, we just set compressibility to 0.0f

    // pose 1
    // advect the center of the Kelvinlet from the start position to end position
    vec3 loadOriginAdvected1 = p.kelvinlet[1].origin + p.kelvinlet[1].linearVelocity * originLerp;
    vec3 R1 = x - loadOriginAdvected1;
    vec3 Kt1 = KTranslation(R1, p.kelvinlet[1].forceVector, p.radius, p.material.stiffness, p.material.compressibility);
    vec3 Kr1 = KTwist(R1, p.kelvinlet[1].twistForceMatrix, p.radius, p.material.stiffness, p.material.compressibility);
    vec3 Ks1 = KScale(R1, p.kelvinlet[1].scaleForceMatrix, p.radius, p.material.stiffness, 0.0f);

    return Kt0 + Kr0 + Ks0 + Kt1 + Kr1 + Ks1;
}

// Returns the displacement vactor for pinch of x
INLINE vec3
KEvaluatePinch(float t, vec3 x, KelvinletParameters p)
{
    vec3 R = x - p.kelvinlet.origin;
    vec3 Kp = KPinch(R, p.kelvinlet.scaleForceMatrix, p.radius, p.material.stiffness, p.material.compressibility);
    return Kp;
}

INLINE vec3
KEvaluatePinchTwoDeformers(float t, vec3 x, KelvinletTwoDeformerParameters p)
{
    // pose 0
    vec3 R0 = x - p.kelvinlet[0].origin;
    vec3 Kp0 = KPinch(R0, p.kelvinlet[0].scaleForceMatrix, p.radius, p.material.stiffness, p.material.compressibility);

    // pose 1
    vec3 R1 = x - p.kelvinlet[1].origin;
    vec3 Kp1 = KPinch(R1, p.kelvinlet[1].scaleForceMatrix, p.radius, p.material.stiffness, p.material.compressibility);

    return Kp0 + Kp1;
}

///////////////////////////////////////////////////////
// The following preprocessor code includes ODESolver multiple times
// to generate different variants of the solvers. This is necessary
// because GLSL doesn't support classes with functions, but we need
// to parameterize the solver with different evaluation functions.
///////////////////////////////////////////////////////

#define SCOPE(suffix) Kelvinlets##suffix
#define evaluate KEvaluate
#define Parameters KelvinletParameters
#include "ODESolvers.h"
#undef Parameters
#undef evaluate
#undef SCOPE

#define SCOPE(suffix) KelvinletsTwoDeformers##suffix
#define evaluate KEvaluateTwoDeformers
#define Parameters KelvinletTwoDeformerParameters
#include "ODESolvers.h"
#undef Parameters
#undef evaluate
#undef SCOPE
