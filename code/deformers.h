// Copyright(c) Facebook, Inc. and its affiliates.
// All rights reserved.
// 
// This source code is licensed under the BSD - style license found in the
// LICENSE file in the root directory of this source tree.

INLINE float calcFalloff(vec3 position, vec3 origin, float innerRadius, float outerRadius)
{
    innerRadius = min(innerRadius, outerRadius - 0.000001f);

    float d = distance(position, origin);
    float t = (d - innerRadius) / (outerRadius - innerRadius);
    float falloff = 1 - saturate(t);
    falloff = smoothstep(0.0f, 1.0f, falloff);

    return falloff;
}

INLINE vec3 deformNonElastic(vec3 position, Deformation deformation, float innerRadius, float outerRadius, float maxerror)
{
    float falloff = calcFalloff(position, deformation.origin, innerRadius, outerRadius);

    if (falloff > 0.0f)
    {
        float endtime = lerp(deformation.time, deformation.time + deformation.dt, falloff);
        position = NonElastic_AdaptiveBS32(position, deformation.time, endtime, deformation, maxerror);
    }

    return position;
}

INLINE vec3 deformBlendedNonElastic(vec3 position, Deformation deformation0, Deformation deformation1, float innerRadius, float outerRadius, float maxerror)
{
    float falloff0 = calcFalloff(position, deformation0.origin, innerRadius, outerRadius);
    float falloff1 = calcFalloff(position, deformation1.origin, innerRadius, outerRadius);

    if (falloff0 > 0.0f || falloff1 > 0.0f)
    {
        NonElasticTwoDeformerParameters parameters;
        parameters.deformation[0] = deformation0;
        parameters.deformation[1] = deformation1;
        float endtime = lerp(deformation0.time, deformation0.time + deformation0.dt, min(falloff0, falloff1));
        position = NonElasticTwoDeformerParameters_AdaptiveBS32(position, deformation0.time, endtime, parameters, maxerror);

        // when the falloff of one deformer is less than the falloff of the second deformer,
        // then we need to integrate from the minfalloff to maxfalloff using the deformer
        // with the larger falloff
        if (falloff0 != falloff1)
        {
            Deformation deformation;
            float starttime = endtime;
            if (falloff0 > falloff1)
            {
                deformation = deformation0;
            } else {
                deformation = deformation1;
            }
            deformation.origin = deformation.origin + deformation.linearVelocity*(endtime - deformation.time);
            endtime = lerp(deformation.time, deformation.time + deformation.dt, max(falloff0, falloff1));
            position = NonElastic_AdaptiveBS32(position, starttime, endtime, deformation, maxerror);
        }
    }

    return position;
}

INLINE vec3 deformNonElastic(vec3 position, const vector<Deformation>& deformations, float innerRadius, float outerRadius, float maxerror)
{
    float falloff = calcFalloff(position, deformations[0].origin, innerRadius, outerRadius);

    if (falloff > 0.0f)
    {
        float starttime = deformations[0].time;
        float endtime = deformations[deformations.size()-1].time + deformations[deformations.size()-1].dt;
        
        float t = starttime;
        float maxt = lerp(starttime, endtime, falloff);
        int frame = 0;
        while (t < maxt)
        {
            float dt = deformations[frame].dt;
            dt = min(maxt - t, dt);

            position = NonElastic_RungeKutta(position, t, t + dt, deformations[frame]);

            frame++;
            t += dt;
        }
    }

    return position;
}

INLINE vec3 deformBlendedNonElastic(vec3 position, const vector<Deformation>& deformations0, const vector<Deformation>& deformations1, float innerRadius, float outerRadius, float maxerror)
{
    float falloff0 = calcFalloff(position, deformations0[0].origin, innerRadius, outerRadius);
    float falloff1 = calcFalloff(position, deformations1[0].origin, innerRadius, outerRadius);

    if (falloff0 > 0 || falloff1 > 0)
    {
        float starttime = deformations0[0].time;
        float endtime = deformations0[deformations0.size() - 1].time + deformations0[deformations0.size() - 1].dt;

        float minfalloff = min(falloff0, falloff1);

        float t = starttime;
        float maxt = lerp(starttime, endtime, minfalloff);
        int frame = 0;
        while (t < maxt)
        {
            float dt = deformations0[frame].dt;
            dt = min(maxt - t, dt);

            NonElasticTwoDeformerParameters parameters;
            parameters.deformation[0] = deformations0[frame];
            parameters.deformation[1] = deformations1[frame];

            position = NonElasticTwoDeformerParameters_RungeKutta(position, t, t + dt, parameters);

            frame++;
            t += dt;
        }

        // when the falloff of one deformer is less than the falloff of the second deformer,
        // then we need to integrate from the minfalloff to maxfalloff using the deformer
        // with the larger falloff
        frame = max(frame - 1, 0);
        float maxfalloff = max(falloff0, falloff1);
        maxt = lerp(starttime, endtime, maxfalloff);
        const vector<Deformation>& deformations = falloff0 > falloff1 ? deformations0 : deformations1;
        while (t < maxt)
        {
            float dt = deformations[frame].dt;
            dt = min(maxt - t, dt);

            Deformation deformation = deformations[frame];
            deformation.origin = deformation.origin + deformation.linearVelocity*(t-deformation.time);

            position = NonElastic_RungeKutta(position, t, t + dt, deformation);

            frame++;
            t += dt;
        }
    }

    return position;
}

INLINE vec3 deformKelvinlet(vec3 position, Kelvinlet kelvinlet, Material material, float radius, float maxerror)
{
    KelvinletParameters parameters;
    parameters.kelvinlet = kelvinlet;
    parameters.material = material;
    parameters.radius = radius;

    position = Kelvinlets_AdaptiveBS32(position, kelvinlet.time, kelvinlet.time + kelvinlet.dt, parameters, maxerror);

    return position;
}

INLINE vec3 deformBlendedKelvinlets(vec3 position, Kelvinlet kelvinlet0, Kelvinlet kelvinlet1, Material material, float radius, float maxerror)
{
    KelvinletTwoDeformerParameters parameters;
    parameters.kelvinlet[0] = kelvinlet0;
    parameters.kelvinlet[1] = kelvinlet1;
    parameters.material = material;
    parameters.radius = radius;

    position = KelvinletsTwoDeformers_AdaptiveBS32(position, kelvinlet0.time, kelvinlet0.time + kelvinlet0.dt, parameters, maxerror);

    return position;
}

INLINE vec3 deformKelvinlet(vec3 position, const vector<Kelvinlet>& kelvinlets, Material material, float radius, float maxerror)
{
    float starttime = kelvinlets[0].time;
    float endtime = kelvinlets[kelvinlets.size() - 1].time + kelvinlets[kelvinlets.size() - 1].dt;

    float t = starttime;
    float maxt = endtime;
    int frame = 0;
    while (t < maxt)
    {
        float dt = kelvinlets[frame].dt;
        dt = min(maxt - t, dt);

        KelvinletParameters parameters;
        parameters.kelvinlet = kelvinlets[frame];
        parameters.material = material;
        parameters.radius = radius;

        position = Kelvinlets_RungeKutta(position, t, t + dt, parameters);

        frame++;
        t += dt;
    }

    return position;
}

INLINE vec3 deformBlendedKelvinlets(vec3 position, const vector<Kelvinlet>& kelvinlets0, const vector<Kelvinlet>& kelvinlets1, Material material, float radius, float maxerror)
{
    float starttime = kelvinlets0[0].time;
    float endtime = kelvinlets0[kelvinlets0.size() - 1].time + kelvinlets0[kelvinlets0.size() - 1].dt;

    float t = starttime;
    float maxt = endtime;
    int frame = 0;
    while (t < maxt)
    {
        float dt = kelvinlets0[frame].dt;
        dt = min(maxt - t, dt);

        KelvinletTwoDeformerParameters parameters;
        parameters.kelvinlet[0] = kelvinlets0[frame];
        parameters.kelvinlet[1] = kelvinlets1[frame];
        parameters.material = material;
        parameters.radius = radius;

        position = KelvinletsTwoDeformers_RungeKutta(position, t, t + dt, parameters);

        frame++;
        t += dt;
    }

    return position;
}
