// Copyright(c) Facebook, Inc. and its affiliates.
// All rights reserved.
// 
// This source code is licensed under the BSD - style license found in the
// LICENSE file in the root directory of this source tree.

// This file contains a collection of ODE solvers
// compatible with C++ and GLSL code. Because this
// code is used by GLSL, we can't use a lot of nice
// things in C, like namespaces, classes, or virtual functions.
// We work around that by generating C style functions using the
// C preprocessor.
//
// To use this code, you must specify these macros, and
// then #include this file:
// SCOPE,
// EVALUATE, 
// PARAMETERLIST,
// PARAMETERS
//
// For example:
//
//  #define SCOPE(suffix) Kelvinlets##suffix
//  #define EVALUATE(...) KEvaluate(__VA_ARGS__)
//  #define PARAMETERLIST KelvinletParameters klvn, State state
//  #define PARAMETERS klvn, state
//  #include "ODESolver.h"
//
// The SCOPE macro must create a unique
// C identifier, given the suffix. This is used
// to provide a scope for function names and structs.
//
// The EVALUATE macro must be a function that
// returns a derivative. It has a declaration of:
// vec3 evaluator(float t, vec3 x, PARAMETERLIST)
//
// The PARAMETERLIST macro is the list of parameters
// used by the EVALUATE function, used in a function
// declaration.
//
// The PARAMETERS macro is the list of parameters
// used by the EVALUATE function, used to invoke the
// function.
//
// The ODE solver about the PARAMETERS, but the user
// needs a way to pass data to the evaluator function.
// For example, if your is f(x) = sin(x) + c + d^x, 
// then c and d would be parameters.

///////////////////////////////////////////////////////
// Tunable Constants
///////////////////////////////////////////////////////

// The adaptive integrators use this as the initial
// step size. This defaults to 10% of the distance
// between the start and end time. It's often faster
// to set this to 1.0, so that the initial step size
// takes you all the way to the end. After all,
// the adaptive integrator should detect that the
// error is too great, and use a smaller step size
// if necessary. Unforunately, an initial dt of 1.0
// does not work well for all problems. This number
// was tuned based on experiments in Medium.
// It should always in the range (0,1].
#define ADAPTIVE_INTEGRATOR_INITIAL_DT 0.1f

// The minimum dt that should be allowed in
// the adaptive integrators. If the adaptive dt
// is a very small value, then something went wrong
// (the problem is stiff?). We don't want to spiral
// into an infinite loop, so when this happens, 
// we use this dt instead of a smaller value.
#define ADAPTIVE_INTEGRATOR_MINIMUM_DT 0.001f

///////////////////////////////////////////////////////
// Integrator step functions
///////////////////////////////////////////////////////

INLINE vec3 
SCOPE(euler)(float t, float dt, vec3 x, PARAMETERLIST)
{
    return x + dt*EVALUATE(t, x, PARAMETERS);
}

INLINE vec3
SCOPE(rungekutta)(float t, float dt, vec3 x, PARAMETERLIST)
{
    float a1 = 0;
    float a2 = 1 / 2.0f;
    float a3 = 1 / 2.0f;
    float a4 = 1.0f;

    float b21 = 1 / 2.0f;
    float b31 = 0;
    float b32 = 1 / 2.0f;
    float b41 = 0;
    float b42 = 0;
    float b43 = 1;

    float c1 = 1 / 6.0f;
    float c2 = 1 / 3.0f;
    float c3 = 1 / 3.0f;
    float c4 = 1 / 6.0f;

    vec3 k1 = dt*EVALUATE(t + dt*a1, x, PARAMETERS);
    vec3 k2 = dt*EVALUATE(t + dt*a2, x + k1*b21, PARAMETERS);
    vec3 k3 = dt*EVALUATE(t + dt*a3, x + k1*b31 + k2*b32, PARAMETERS);
    vec3 k4 = dt*EVALUATE(t + dt*a4, x + k1*b41 + k2*b42 + k3*b43, PARAMETERS);

    return x + k1*c1 + k2*c2 + k3*c3 + k4*c4;
}

// Runge-Kutta-Fehlberg method
// This is an embedded method that computes a fourth and fifth order answer
// The original RK45 method computed a fourth order answer and a fifth order error estimator.
// We use 'local extrapolation', meaning we use the fifth order answer for our answer, and use
// the difference between the two terms for a fourth order accurate error.
struct SCOPE(RK45Result )
{
    vec3 fourthorder; 
    vec3 fifthorder;
};
INLINE SCOPE(RK45Result)
SCOPE(rungekuttafehlberg)(float t, float dt, vec3 x, PARAMETERLIST)
{
    // https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method
    float a1 = 0;
    float a2 = 1 / 4.0f;
    float a3 = 3 / 8.0f;
    float a4 = 12 / 13.0f;
    float a5 = 1.0f;
    float a6 = 1 / 2.0f;

    float b21 = 1 / 4.0f;
    float b31 = 3 / 32.0f;
    float b32 = 9 / 32.0f;
    float b41 = 1932 / 2197.0f;
    float b42 = -7200 / 2197.0f;
    float b43 = 7296 / 2197.0f;
    float b51 = 439 / 216.0f;
    float b52 = -8.0f;
    float b53 = 3680 / 513.0f;
    float b54 = -845 / 4104.0f;
    float b61 = -8 / 27.0f;
    float b62 = 2.0f;
    float b63 = -3544 / 2565.0f;
    float b64 = 1859 / 4104.0f;
    float b65 = -11 / 40.0f;

    float c1 = 16 / 135.0f;      // fifth order answer
    float c2 = 0.0f;
    float c3 = 6656 / 12825.0f;
    float c4 = 28561 / 56430.0f;
    float c5 = -9 / 50.0f;
    float c6 = 2 / 55.0f;

    float d1 = 25 / 216.0f;      // fourth order answer
    float d2 = 0;
    float d3 = 1408 / 2565.0f;
    float d4 = 2197 / 4104.0f;
    float d5 = -1 / 5.0f;
    float d6 = 0;

    vec3 k1 = dt*EVALUATE(t + dt*a1, x, PARAMETERS);
    vec3 k2 = dt*EVALUATE(t + dt*a2, x + k1*b21, PARAMETERS);
    vec3 k3 = dt*EVALUATE(t + dt*a3, x + k1*b31 + k2*b32, PARAMETERS);
    vec3 k4 = dt*EVALUATE(t + dt*a4, x + k1*b41 + k2*b42 + k3*b43, PARAMETERS);
    vec3 k5 = dt*EVALUATE(t + dt*a5, x + k1*b51 + k2*b52 + k3*b53 + k4*b54, PARAMETERS);
    vec3 k6 = dt*EVALUATE(t + dt*a6, x + k1*b61 + k2*b62 + k3*b63 + k4*b64 + k5*b65, PARAMETERS);

    SCOPE(RK45Result) result;
    result.fourthorder = x + k1*d1 + k2*d2 + k3*d3 + k4*d4 + k5*d5 + k6*d6;
    result.fifthorder = x + k1*c1 + k2*c2 + k3*c3 + k4*c4 + k5*c5 + k6*c6;
    return result;
}

// Runge-Kutta-Dormand-Prince method
// This is an embedded method that computes a fourth and fifth order answer
// Dormand and Prince wrote this to produce a fifth order answer with less error than RK45.
// It is a FSAL method (First Same As Last) so that you can rearrange this code to use 6 evaluations
// per step instead of 7 (although we don't take advantage of that yet).
struct SCOPE(DormandPrinceRungeKuttaResult)
{
    vec3 fourthorder;
    vec3 fifthorder;
};
INLINE SCOPE(DormandPrinceRungeKuttaResult)
SCOPE(dormandprincerungekutta)(float t, float dt, vec3 x, PARAMETERLIST)
{
    // https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method
    float a1 = 0;
    float a2 = 1 / 5.0f;
    float a3 = 3 / 10.0f;
    float a4 = 4 / 5.0f;
    float a5 = 8 / 9.0f;
    float a6 = 1.0f;
    float a7 = 1.0f;

    float b21 = 1 / 5.0f;
    float b31 = 3 / 40.0f;
    float b32 = 9 / 40.0f;
    float b41 = 44 / 45.0f;
    float b42 = -56 / 15.0f;
    float b43 = 32 / 9.0f;
    float b51 = 19372 / 6561.0f;
    float b52 = -25360 / 2187.0f;
    float b53 = 64448 / 6561.0f;
    float b54 = -212 / 729.0f;
    float b61 = 9017 / 3168.0f;
    float b62 = -355 / 33.0f;
    float b63 = 46732 / 5247.0f;
    float b64 = 49 / 176.0f;
    float b65 = -5103 / 18656.0f;
    float b71 = 35 / 384.0f;
    float b72 = 0.0f;
    float b73 = 500 / 1113.0f;
    float b74 = 125 / 192.0f;
    float b75 = -2187 / 6784.0f;
    float b76 = 11 / 84.0f;

    float c1 = 35 / 384.0f;    // fifth order answer
    float c2 = 0;
    float c3 = 500 / 1113.0f;
    float c4 = 125 / 192.0f;
    float c5 = -2187 / 6784.0f;
    float c6 = 11 / 84.0f;

    float d1 = 5179 / 57600.0f;    // fourth order answer
    float d2 = 0;
    float d3 = 7571 / 16695.0f;
    float d4 = 393 / 640.0f;
    float d5 = -92097 / 339200.0f;
    float d6 = 187 / 2100.0f;
    float d7 = 1 / 40.0f;

    vec3 k1 = dt*EVALUATE(t + dt*a1, x, PARAMETERS);
    vec3 k2 = dt*EVALUATE(t + dt*a2, x + k1*b21, PARAMETERS);
    vec3 k3 = dt*EVALUATE(t + dt*a3, x + k1*b31 + k2*b32, PARAMETERS);
    vec3 k4 = dt*EVALUATE(t + dt*a4, x + k1*b41 + k2*b42 + k3*b43, PARAMETERS);
    vec3 k5 = dt*EVALUATE(t + dt*a5, x + k1*b51 + k2*b52 + k3*b53 + k4*b54, PARAMETERS);
    vec3 k6 = dt*EVALUATE(t + dt*a6, x + k1*b61 + k2*b62 + k3*b63 + k4*b64 + k5*b65, PARAMETERS);
    vec3 k7 = dt*EVALUATE(t + dt*a7, x + k1*b71 + k2*b72 + k3*b73 + k4*b74 + k5*b75 + k5*b76, PARAMETERS);

    SCOPE(DormandPrinceRungeKuttaResult) result;
	result.fourthorder = x + k1*d1 + k2*d2 + k3*d3 + k4*d4 + k5*d5 + k6*d6 + k7*d7;
	result.fifthorder = x + k1*c1 + k2*c2 + k3*c3 + k4*c4 + k5*c5 + k6*c6;
    return result;
}

// Bogacki-Shampine 3(2) method
// This is an embedded method that computes a second and third order answer
// Given the same error tolerance, BS32 takes more steps than RKF45 or DP54, 
// but it is often cheaper to compute for each step since it takes fewer evaluations.
struct SCOPE(BogackiShampineRungeKuttaResult)
{
    vec3 secondorder;
    vec3 thirdorder;
};
INLINE SCOPE(BogackiShampineRungeKuttaResult)
SCOPE(bogackishampinerungekutta)(float t, float dt, vec3 x, PARAMETERLIST)
{
    // https://en.wikipedia.org/wiki/Bogacki%E2%80%93Shampine_method
    float a1 = 0;
    float a2 = 1 / 2.0f;
    float a3 = 3 / 4.0f;
    float a4 = 1.0f;

    float b21 = 1 / 2.0f;
    float b31 = 0;
    float b32 = 3 / 4.0f;
    float b41 = 2 / 9.0f;
    float b42 = 1 / 3.0f;
    float b43 = 4 / 9.0f;

    float c1 = 2 / 9.0f;    // third order answer
    float c2 = 1 / 3.0f;
    float c3 = 4 / 9.0f;

    float d1 = 7 / 24.0f;    // second order answer
    float d2 = 1 / 4.0f;
    float d3 = 1 / 3.0f;
    float d4 = 1 / 8.0f;

    vec3 k1 = dt*EVALUATE(t + dt*a1, x, PARAMETERS);
    vec3 k2 = dt*EVALUATE(t + dt*a2, x + k1*b21, PARAMETERS);
    vec3 k3 = dt*EVALUATE(t + dt*a3, x + k1*b31 + k2*b32, PARAMETERS);
    vec3 k4 = dt*EVALUATE(t + dt*a4, x + k1*b41 + k2*b42 + k3*b43, PARAMETERS);

    SCOPE(BogackiShampineRungeKuttaResult) result;
    result.secondorder = x + k1*d1 + k2*d2 + k3*d3 + k4*d4;
    result.thirdorder = x + k1*c1 + k2*c2 + k3*c3;
    return result;
}

///////////////////////////////////////////////////////
// Solvers
///////////////////////////////////////////////////////

// This takes 100 Euler steps. Included as a simple example.
INLINE vec3 SCOPE(_FixedEuler)(vec3 pos, float tstart, float tend, PARAMETERLIST)
{
    float t = tstart;
    float dt = (tend - tstart) * 0.01f;
    while (t < tend)
    {
        pos = SCOPE(euler)(t, dt, pos, PARAMETERS);
        t += dt;
    }

    return pos;
}

// This takes ten RK4 steps. Included as a simple example.
INLINE vec3 SCOPE(_FixedRungeKutta)(vec3 pos, float tstart, float tend, PARAMETERLIST)
{
    float t = tstart;
    float dt = (tend - tstart) * 0.1f;
    while (t < tend)
    {
        pos = SCOPE(rungekutta)(t, dt, pos, PARAMETERS);
        t += dt;
    }

    return pos;
}

// This takes a single RK4 step.
INLINE vec3 SCOPE(_RungeKutta)(vec3 pos, float tstart, float tend, PARAMETERLIST)
{
    float t = tstart;
    float dt = (tend - tstart);

    pos = SCOPE(rungekutta)(t, dt, pos, PARAMETERS);

    return pos;
}

INLINE vec3 SCOPE(_AdaptiveRK)(vec3 pos, float tstart, float tend, float maxerror, PARAMETERLIST)
{
    // Runge Kutta, adaptive by taking a step and then two half-steps

    float t = tstart;
    float dt = (tend - tstart) * ADAPTIVE_INTEGRATOR_INITIAL_DT;
    while (t < tend)
    {
        dt = min(dt, tend - t);

        vec3 fullstep = SCOPE(rungekutta)(t, dt, pos, PARAMETERS);
        vec3 halfstep = SCOPE(rungekutta)(t, dt / 2.0f, pos, PARAMETERS);
        vec3 twohalfsteps = SCOPE(rungekutta)(t + dt / 2.0f, dt / 2.0f, halfstep, PARAMETERS);

        // step doubling uses the formula (fullstep - twohalfsteps)/15
        // which is more accurate than the embedded methods, which use (fifthorderanswer - fourthorderanswer) as error
        // however, step doubling takes more calculations to compute
        float error = length(fullstep - twohalfsteps) / 15.0f / dt;

        float safety = 0.9f;
        float newdt = dt * safety * pow(maxerror / error, 0.2f);

        if (error <= maxerror || dt <= ADAPTIVE_INTEGRATOR_MINIMUM_DT)
        {
            pos = fullstep;
            t += dt;
            dt = newdt;
        }
        else
        {
            // we have a lot of error. if the new dt is (nearly) the same as our current dt,
            // then set a fixed step size to prevent an infinite loop
            dt = (abs(newdt - dt) < 0.00001f) ? dt / 2.f : newdt;
            dt = max(dt, ADAPTIVE_INTEGRATOR_MINIMUM_DT);
        }
    }

    return pos;
}

INLINE vec3 SCOPE(_AdaptiveRKF45)(vec3 pos, float tstart, float tend, float maxerror, PARAMETERLIST)
{
    float t = tstart;
    float dt = (tend - tstart) * ADAPTIVE_INTEGRATOR_INITIAL_DT;
    while (t < tend)
    {
        dt = min(dt, tend - t);

        SCOPE(RK45Result) rk45 = SCOPE(rungekuttafehlberg)(t, dt, pos, PARAMETERS);

        float error = length(rk45.fifthorder - rk45.fourthorder) / dt;

        float safety = 0.9f;
        float newdt = dt * safety * pow(maxerror / error, 0.20f);

        if (error <= maxerror || dt <= ADAPTIVE_INTEGRATOR_MINIMUM_DT)
        {
            pos = rk45.fifthorder;  // local extrapolation
            t += dt;
            dt = newdt;

        }
        else
        {
            // we have a lot of error. if the new dt is (nearly) the same as our current dt,
            // then use half our step size to prevent an infinite loop
            dt = (abs(newdt - dt) < 0.00001f) ? dt / 2.f : newdt;
            dt = max(dt, ADAPTIVE_INTEGRATOR_MINIMUM_DT);
        }
    }

    return pos;
}

INLINE vec3 SCOPE(_AdaptiveDP54)(vec3 pos, float tstart, float tend, float maxerror, PARAMETERLIST)
{
    float t = tstart;
    float dt = (tend - tstart) * ADAPTIVE_INTEGRATOR_INITIAL_DT;
    while (t < tend)
    {
        dt = min(dt, tend - t);

        SCOPE(DormandPrinceRungeKuttaResult) dprk = SCOPE(dormandprincerungekutta)(t, dt, pos, PARAMETERS);

        float error = length(dprk.fifthorder - dprk.fourthorder) / dt;

        float safety = 0.9f;
        float newdt = dt * safety * pow(maxerror / error, 0.2f);

        if (error <= maxerror || dt <= ADAPTIVE_INTEGRATOR_MINIMUM_DT)
        {
            pos = dprk.fifthorder;    // local extrapolation
            t += dt;
            dt = newdt;
        }
        else
        {
            // we have a lot of error. if the new dt is (nearly) the same as our current dt,
            // then use half our step size to prevent an infinite loop
            dt = (abs(newdt - dt) < 0.00001f) ? dt / 2.f : newdt;
            dt = max(dt, ADAPTIVE_INTEGRATOR_MINIMUM_DT);
        }
    }

    return pos;
}

INLINE vec3 SCOPE(_AdaptiveBS32)(vec3 pos, float tstart, float tend, float maxerror, PARAMETERLIST)
{
    float t = tstart;
    float dt = (tend - tstart) * ADAPTIVE_INTEGRATOR_INITIAL_DT;
    while (t < tend)
    {
        dt = min(dt, tend - t);

        SCOPE(BogackiShampineRungeKuttaResult) bsrk = SCOPE(bogackishampinerungekutta)(t, dt, pos, PARAMETERS);

        float error = length(bsrk.thirdorder - bsrk.secondorder) / dt;

        float safety = 0.9f;
        float newdt = dt * safety * pow(maxerror / error, 1/3.0f);

        if (error <= maxerror || dt <= ADAPTIVE_INTEGRATOR_MINIMUM_DT)
        {
            pos = bsrk.thirdorder;    // local extrapolation
            t += dt;
            dt = newdt;
        }
        else
        {
            // we have a lot of error. if the new dt is (nearly) the same as our current dt,
            // then use half our step size to prevent an infinite loop
            dt = (abs(newdt - dt) < 0.00001f) ? dt / 2.f : newdt;
            dt = max(dt, ADAPTIVE_INTEGRATOR_MINIMUM_DT);
        }
    }

    return pos;
}
