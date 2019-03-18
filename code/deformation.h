// Copyright(c) Facebook, Inc. and its affiliates.
// All rights reserved.
// 
// This source code is licensed under the BSD - style license found in the
// LICENSE file in the root directory of this source tree.

#ifdef __cplusplus
#pragma once
#define INLINE inline
INLINE float min(float a, float b) { return a < b ? a : b; };
INLINE float max(float a, float b) { return a > b ? a : b; };
INLINE int min(int a, int b) { return a < b ? a : b; };
INLINE int max(int a, int b) { return a > b ? a : b; };
#include "glslmathforcpp.h"
#else
// otherwise this assumes GLSL, which defines min()/max()
#define INLINE
struct quat { float r, i, j, k; };
#endif

struct Pose
{
    vec3  position;
    quat  orientation;
    float scale;
    float time;
};

struct Motion
{
    vec3  origin;
    vec3  linearVelocity;
    vec3  angularVelocity;
    float scaleFactor;
    float time;
    float dt;
};

struct Deformation
{
    vec3   origin;
    vec3   linearVelocity;
    mat3x3 rotationTensor;
    mat3x3 strainTensor;
    mat3x3 displacementGradientTensor;
    float  time;
    float  dt;
};

struct Kelvinlet
{
    vec3   origin;
    vec3   linearVelocity;
    vec3   forceVector;
    mat3x3 twistForceMatrix;
    mat3x3 scaleForceMatrix;
    float  time;
    float  dt;
	float  radius;
	float  stiffness;
	float  compressibility;
};

#include "kelvinlets.h"
#include "nonelastic.h"

// The code below is meant to be run on the CPU
// It does not need to be computed per-vertex
#ifdef __cplusplus

INLINE Motion buildMotion(Pose start, Pose end);

INLINE Deformation buildDeformation(Motion motion);

INLINE Kelvinlet buildKelvinlet(Deformation deformation, float stiffness, float compressibilty, float radius);

INLINE mat3x3 skewSymmetric(vec3 a)
{
    // This assumes mat3x3's constructor takes three column vectors
    // e.g. the zeroth column is (0, a.z, -a.y). When reading this,
    // you should transpose the vec3's and see them as columns.
    return mat3x3
    {
        vec3(   0,  a.z, -a.y),
        vec3(-a.z,    0,  a.x),
        vec3( a.y, -a.x,    0)
    };
}

INLINE mat3x3 identityMat3x3()
{
    return mat3x3
    {
        vec3(1, 0, 0),
        vec3(0, 1, 0),
        vec3(0, 0, 1)
    };
}

INLINE Motion buildMotion(Pose start, Pose end)
{
    Motion motion;

    quat rotation = end.orientation * inverse(start.orientation);
    vec3 axis;
    float angle;
    rotation.toAxisAngle(axis, angle);

    if (angle == 0.0f)
    {
        axis = vec3(1, 0, 0);
        angle = 0.0f;
    }

    vec3  t = end.position - start.position;
    vec3  r = axis * angle;
    float s = end.scale / start.scale;

    float dt = end.time - start.time;

    motion.origin = start.position;
    motion.linearVelocity = t / dt;
    motion.angularVelocity = r / dt;
    motion.scaleFactor = s;
    motion.time = start.time;
    motion.dt = dt;

    return motion;
}

INLINE Deformation buildDeformation(Motion motion)
{
    Deformation deformation;

    deformation.origin = motion.origin;
    deformation.linearVelocity = motion.linearVelocity;
    deformation.rotationTensor = skewSymmetric(motion.angularVelocity);
    deformation.strainTensor = identityMat3x3() * log(pow(motion.scaleFactor, 1/motion.dt));
    deformation.displacementGradientTensor = deformation.rotationTensor + deformation.strainTensor;
    deformation.time = motion.time;
    deformation.dt = motion.dt;
    
    return deformation;
}

INLINE Kelvinlet buildKelvinlet(Deformation deformation, float stiffness, float compressibility, float radius)
{
    Kelvinlet kelvinlet;

    // Compressibility of 0.5 causes divide by 0 in Kelvinlets equation for scale,
    // so use 0.0 (also done in KEvaluate() function)
    float scaleCompressibility = 0.0f;

    kelvinlet.origin = deformation.origin;
    kelvinlet.linearVelocity = deformation.linearVelocity;
    kelvinlet.forceVector = deformation.linearVelocity * KTranslationCalibrationFactor(radius, compressibility);
    kelvinlet.twistForceMatrix = deformation.rotationTensor * KTwistCalibrationFactor(radius, compressibility);
    kelvinlet.scaleForceMatrix = deformation.strainTensor * KScaleCalibrationFactor(radius, scaleCompressibility);
    kelvinlet.time = deformation.time;
    kelvinlet.dt = deformation.dt;
	kelvinlet.radius = radius;
	kelvinlet.stiffness = stiffness;
	kelvinlet.compressibility = compressibility;
    
    return kelvinlet;
}

#endif
