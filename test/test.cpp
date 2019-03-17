// Copyright(c) Facebook, Inc. and its affiliates.
// All rights reserved.
// 
// This source code is licensed under the BSD - style license found in the
// LICENSE file in the root directory of this source tree.

#include <cmath>
namespace deformation
{
    #include "../code/deformation.h"
};

#include <vector>
#include "assert.h"

using namespace std;
using namespace deformation;

using uint = unsigned int;

struct Vertex
{
    vec3 position;
};

struct Mesh
{
    vector<vec3> vertices;
    vector<uint> indices;
};

const uint datafileversion = 0;

Mesh readmesh(const char* filename)
{
    Mesh mesh;

    FILE* file;
    fopen_s(&file, filename, "rb");
    assert(file);

    uint version;
    fread(&version, 4, 1, file);
    assert(version == datafileversion);

    uint numvertices;
    fread(&numvertices, 4, 1, file);
    uint numindices;
    fread(&numindices, 4, 1, file);

    mesh.vertices.resize(numvertices);
    for (uint i = 0; i < numvertices; i++)
    {
        fread(&mesh.vertices[i], 12, 1, file);
    }

    mesh.indices.resize(numindices);
    for (uint i = 0; i < numindices; i++)
    {
        fread(&mesh.indices[i], 4, 1, file);
    }

    fclose(file);

    return mesh;
}

// A "stroke" is a set of poses that was recorded in Medium
// so that we have some test data.
struct Stroke
{
    vector<deformation::Pose> poses;
    float innerRadius;
    float outerRadius;
    bool  elastic;
    float stiffness;
    float compressibility;
};

Stroke readstroke(const char* filename)
{
    Stroke stroke;

    FILE* file;
    fopen_s(&file, filename, "rb");
    assert(file);

    uint version;
    fread(&version, 4, 1, file);
    assert(version == datafileversion);

    uint numposes;
    fread(&numposes, 4, 1, file);
    stroke.poses.resize(numposes);
    for (uint i = 0; i < numposes; i++)
    {
        fread(&stroke.poses[i].position, 12, 1, file);
        fread(&stroke.poses[i].orientation, 16, 1, file);
        fread(&stroke.poses[i].scale, 4, 1, file);
        fread(&stroke.poses[i].time, 4, 1, file);
    }
    fread(&stroke.innerRadius, 4, 1, file);
    fread(&stroke.outerRadius, 4, 1, file);
    uint elastic;
    fread(&elastic, 4, 1, file);
    stroke.elastic = elastic ? true : false;
    fread(&stroke.stiffness, 4, 1, file);
    fread(&stroke.compressibility, 4, 1, file);

    fclose(file);

    return stroke;
}

void writeobj(const char* filename, Mesh mesh)
{
    FILE* file;
    fopen_s(&file, filename, "wt");
    assert(file);

    // calculate vertex normals
    vector<vec3> vn;
    vn.resize(mesh.vertices.size());
    memset(vn.data(), 0, sizeof(vec3)*mesh.vertices.size());
    for (uint i = 0; i < mesh.indices.size(); i += 3)
    {
        vec3 v0 = mesh.vertices[mesh.indices[i + 0]];
        vec3 v1 = mesh.vertices[mesh.indices[i + 1]];
        vec3 v2 = mesh.vertices[mesh.indices[i + 2]];
        vec3 n = cross(v1 - v0, v2 - v0);
        vn[mesh.indices[i + 0]] = vn[mesh.indices[i + 0]] + n;
        vn[mesh.indices[i + 1]] = vn[mesh.indices[i + 0]] + n;
        vn[mesh.indices[i + 2]] = vn[mesh.indices[i + 0]] + n;
    }
    for (uint i = 0; i < mesh.vertices.size(); i++)
    {
        vn[i] = normalize(vn[i]);
    }

    // write out .obj
    for (uint i = 0; i < mesh.vertices.size(); i++)
    {
        fprintf(file, "v %f %f %f\n", mesh.vertices[i].x, mesh.vertices[i].y, mesh.vertices[i].z);
    }

    for (uint i = 0; i < mesh.vertices.size(); i++)
    {
        fprintf(file, "vn %f %f %f\n", vn[i].x, vn[i].y, vn[i].z);
    }

    for (uint i = 0; i < mesh.indices.size(); i+=3)
    {
        fprintf(file, "f %d//%d %d//%d %d//%d\n", 
            mesh.indices[i+0]+1, mesh.indices[i + 0] + 1, 
            mesh.indices[i+1]+1, mesh.indices[i + 1] + 1, 
            mesh.indices[i+2]+1, mesh.indices[i + 2] + 1);
    }

    fclose(file);
}

// fix any flips caused by quaternion double cover
vector<deformation::Pose> fixFlips(vector<deformation::Pose> poses)
{
    vector<deformation::Pose> fixedPoses = poses;

    for (uint i = 1; i < fixedPoses.size(); i++)
    {
        if (dot(fixedPoses[i - 1].orientation, fixedPoses[i].orientation) < 0)
        {
            fixedPoses[i].orientation *= -1;
        }
    }

    return fixedPoses;
}

// trim intermediate keyframes so only the start and end keyframe poses remain
vector<deformation::Pose> buildStartEndPoses(vector<deformation::Pose> poses)
{
    vector<deformation::Pose> startendposes;
    startendposes.push_back(poses[0]);
    startendposes.push_back(poses[poses.size() - 1]);
    return startendposes;
}

// build motion/deformation/kelvinlet data from poses
// each pair of poses has one motion/deformation/kelvinlet;
// for N poses, there are N-1 motion/deformation/kelvinlets
struct DataFromPoses
{
    vector<deformation::Motion> motions;
    vector<deformation::Deformation> deformations;
    vector<deformation::Kelvinlet> kelvinlets;
};

DataFromPoses buildDataFromPoses(Stroke stroke)
{
    DataFromPoses d;

    for (uint i = 0; i < stroke.poses.size() - 1; i++)
    {
        d.motions.push_back(deformation::buildMotion(stroke.poses[i], stroke.poses[i + 1]));
        d.deformations.push_back(deformation::buildDeformation(d.motions[i]));
        d.kelvinlets.push_back(deformation::buildKelvinlet(d.deformations[i], stroke.stiffness, stroke.compressibility, stroke.outerRadius));
    }

    return d;
}

// Calculate falloff. When the position is inside the inner radius, the falloff is 1;
// when the position is outside the outer radius, the falloff is 0. In between, the
// falloff smoothly falls off from 1 to 0.
float calcFalloff(vec3 position, vec3 origin, float innerRadius, float outerRadius)
{
	innerRadius = min(innerRadius, outerRadius - 0.000001f);

	float d = distance(position, origin);
	float t = (d - innerRadius) / (outerRadius - innerRadius);
	float falloff = 1 - saturate(t);
	falloff = smoothstep(0.0f, 1.0f, falloff);

	return falloff;
}

int main()
{
    // This is the same maxerror factor used in Medium
    // This is 0.1 of the voxel size in Medium, in world space
    // Medium uses meters for world space units
    // Different apps will have different max error requirements
    // and there is no one-size-fits-all number
    float maxerror = 0.00013f;

    // There are eight tests below this
    // They read in some data created from Medium, apply deformers,
    // and write out the results as testresult*.obj

    // --------------------
    // This demonstrates deforming with a start and end pose, like Medium's move tool.
    // This takes a stream of recorded poses, then uses the first and last one.
    // It uses adaptive integration to integrate between these two poses.
    // The deformation style is nonelastic (an affine transformation).
    // In Medium, when using the nonelastic Move Tool, this runs every frame
    // in a vertex shader. Here, we're doing it on the CPU as an example.
    if (true)
    {
        Mesh mesh = readmesh("data\\meshes\\test0_mesh.bin");
        Stroke stroke = readstroke("data\\strokes\\test0_righthandstroke.bin");

        stroke.poses = fixFlips(stroke.poses);

        stroke.poses = buildStartEndPoses(stroke.poses);

        deformation::Motion motion = buildMotion(stroke.poses[0], stroke.poses[1]);
        deformation::Deformation deformation = buildDeformation(motion);

        for (uint i = 0; i < mesh.vertices.size(); i++)
        {
			float falloff = calcFalloff(mesh.vertices[i], deformation.origin, stroke.innerRadius, stroke.outerRadius);

			if (falloff > 0.0f)
			{
				float endtime = lerp(deformation.time, deformation.time + deformation.dt, falloff);
				mesh.vertices[i] = IntegrateNonElastic_AdaptiveBS32(mesh.vertices[i], deformation.time, endtime, maxerror, deformation);
			}
        }

        writeobj("data\\testresult0.obj", mesh);
        printf("test0 success\n");
    }

    // --------------------
    // This demonstrates deforming with a start and end pose, like Medium's move tool.
    // This uses two streams of recorded poses, one for each hand.
    // It then uses the first and last pose for each hand.
    // It uses adaptive integration to integrate between these two poses.
    // This deforms with two blended nonelastic deformers.
    // In Medium, when using the nonelastic Move Tool, this runs every frame
    // in a vertex shader. Here, we're doing it on the CPU as an example.
    if (true)
    {
        Mesh mesh = readmesh("data\\meshes\\test1_mesh.bin");
        Stroke strokerighthand = readstroke("data\\strokes\\test1_righthandstroke.bin");
        Stroke strokelefthand = readstroke("data\\strokes\\test1_lefthandstroke.bin");

        strokerighthand.poses = fixFlips(strokerighthand.poses);
        strokelefthand.poses = fixFlips(strokelefthand.poses);

        strokerighthand.poses = buildStartEndPoses(strokerighthand.poses);
        strokelefthand.poses = buildStartEndPoses(strokelefthand.poses);

        deformation::Motion righthandmotion = buildMotion(strokerighthand.poses[0], strokerighthand.poses[1]);
        deformation::Deformation righthanddeformation = buildDeformation(righthandmotion);

        deformation::Motion lefthandmotion = buildMotion(strokelefthand.poses[0], strokelefthand.poses[1]);
        deformation::Deformation lefthanddeformation = buildDeformation(lefthandmotion);

		for (uint i = 0; i < mesh.vertices.size(); i++)
		{
			float falloff0 = calcFalloff(mesh.vertices[i], righthanddeformation.origin, strokerighthand.innerRadius, strokerighthand.outerRadius);
			float falloff1 = calcFalloff(mesh.vertices[i], lefthanddeformation.origin, strokelefthand.innerRadius, strokelefthand.outerRadius);

			if (falloff0 > 0.0f || falloff1 > 0.0f)
			{
				float endtime = lerp(righthanddeformation.time, righthanddeformation.time + righthanddeformation.dt, min(falloff0, falloff1));
				mesh.vertices[i] = IntegrateNonElasticTwoDeformers_AdaptiveBS32(mesh.vertices[i], righthanddeformation.time, endtime, maxerror, righthanddeformation, lefthanddeformation);

				// when the falloff of one deformer is less than the falloff of the second deformer,
				// then we need to integrate from the minfalloff to maxfalloff using the deformer
				// with the larger falloff
				if (falloff0 != falloff1)
				{
					deformation::Deformation deformation;
					float starttime = endtime;
					if (falloff0 > falloff1)
					{
						deformation = righthanddeformation;
					}
					else {
						deformation = lefthanddeformation;
					}
					deformation.origin = deformation.origin + deformation.linearVelocity*(endtime - deformation.time);
					endtime = lerp(deformation.time, deformation.time + deformation.dt, max(falloff0, falloff1));
					mesh.vertices[i] = IntegrateNonElastic_AdaptiveBS32(mesh.vertices[i], starttime, endtime, maxerror, deformation);
				}
			}
		}

        writeobj("data\\testresult1.obj", mesh);
        printf("test1 success\n");
    }

    // --------------------
    // This demonstrates deforming with a start and end pose, like Medium.
    // This takes a stream of recorded poses, then uses the first and last one.
    // It uses adaptive integration to integrate between these two poses.
    // This deforms with a single Kelvinlet.
    // In Medium, when using the elastic Move Tool, this runs every frame
    // in a vertex shader. Here, we're doing it on the CPU as an example.
    if (true)
    {
        Mesh mesh = readmesh("data\\meshes\\test0_mesh.bin");
        Stroke stroke = readstroke("data\\strokes\\test0_righthandstroke.bin");

        stroke.poses = fixFlips(stroke.poses);

        // Use the first and last pose as our two poses
        stroke.poses = buildStartEndPoses(stroke.poses);

        deformation::Motion motion = buildMotion(stroke.poses[0], stroke.poses[1]);
        deformation::Deformation deformation = buildDeformation(motion);
        deformation::Kelvinlet kelvinlet = buildKelvinlet(deformation, stroke.stiffness, stroke.compressibility, stroke.outerRadius);

        for (uint i = 0; i < mesh.vertices.size(); i++)
        {
			mesh.vertices[i] = IntegrateKelvinlets_AdaptiveBS32(mesh.vertices[i], kelvinlet.time, kelvinlet.time + kelvinlet.dt, maxerror, kelvinlet);
		}

        writeobj("data\\testresult2.obj", mesh);
        printf("test2 success\n");
    }

    // --------------------
    // This demonstrates deforming with a start and end pose, like Medium.
    // This uses two streams of recorded poses, one for each hand.
    // It then uses the first and last pose for each hand.
    // It uses adaptive integration to integrate between these two poses.
    // This deforms with two Kelvinlets.
    // In Medium, when using the elastic Move Tool, this runs every frame
    // in a vertex shader. Here, we're doing it on the CPU as an example.
    if (true)
    {
        Mesh mesh = readmesh("data\\meshes\\test1_mesh.bin");
        Stroke strokerighthand = readstroke("data\\strokes\\test1_righthandstroke.bin");
        Stroke strokelefthand = readstroke("data\\strokes\\test1_lefthandstroke.bin");

        strokerighthand.poses = fixFlips(strokerighthand.poses);
        strokelefthand.poses = fixFlips(strokelefthand.poses);

        strokerighthand.poses = buildStartEndPoses(strokerighthand.poses);
        strokelefthand.poses = buildStartEndPoses(strokelefthand.poses);

        deformation::Motion motionrighthand = buildMotion(strokerighthand.poses[0], strokerighthand.poses[1]);
        deformation::Deformation deformationrighthand = buildDeformation(motionrighthand);
        deformation::Kelvinlet kelvinletrighthand = buildKelvinlet(deformationrighthand, strokerighthand.stiffness, strokerighthand.compressibility, strokerighthand.outerRadius);

        deformation::Motion motionlefthand = buildMotion(strokelefthand.poses[0], strokelefthand.poses[1]);
        deformation::Deformation deformationlefthand = buildDeformation(motionlefthand);
        deformation::Kelvinlet kelvinletlefthand = buildKelvinlet(deformationlefthand, strokelefthand.stiffness, strokelefthand.compressibility, strokelefthand.outerRadius);

        for (uint i = 0; i < mesh.vertices.size(); i++)
        {
			mesh.vertices[i] = IntegrateKelvinletsTwoDeformers_AdaptiveBS32(mesh.vertices[i], kelvinletrighthand.time, kelvinletrighthand.time + kelvinletrighthand.dt, maxerror, kelvinletrighthand, kelvinletlefthand);
		}

        writeobj("data\\testresult3.obj", mesh);
        printf("test3 success\n");
    }

    // --------------------
    // This demonstrates deforming with motion capture data.
    // This was recorded in Medium by writing out the Pose structure every frame.
    // This deforms with a single nonelastic deformer.
	// This does a single RK4 integration per pair of keyframes.
	//
    // You probably don't want to do this on the CPU like this. Instead of accumulating
    // keyframe poses over time and applying them all at once, instead do a single Runge-Kutta
    // integration each frame. That will be much faster than this-- even cheaper on the GPU
    // than the adaptive integrators. However, this achieves the same results.
    if (true)
    {
        Mesh mesh = readmesh("data\\meshes\\test0_mesh.bin");
        Stroke stroke = readstroke("data\\strokes\\test0_righthandstroke.bin");

        stroke.poses = fixFlips(stroke.poses);
        DataFromPoses data = buildDataFromPoses(stroke);

        for (uint i = 0; i < mesh.vertices.size(); i++)
        {
			float falloff = calcFalloff(mesh.vertices[i], data.deformations[0].origin, stroke.innerRadius, stroke.outerRadius);

			if (falloff > 0.0f)
			{
				float starttime = data.deformations[0].time;
				float endtime = data.deformations[data.deformations.size() - 1].time + data.deformations[data.deformations.size() - 1].dt;

				float t = starttime;
				float maxt = lerp(starttime, endtime, falloff);
				int frame = 0;
				while (t < maxt)
				{
					float dt = data.deformations[frame].dt;
					dt = min(maxt - t, dt);

					mesh.vertices[i] = IntegrateNonElastic_RungeKutta(mesh.vertices[i], t, t + dt, data.deformations[frame]);

					frame++;
					t += dt;
				}
			}
		}

        writeobj("data\\testresult4.obj", mesh);
        printf("test4 success\n");
    }

    // --------------------
    // This demonstrates deforming with motion capture data.
    // Two streams of poses were recorded in Medium, one for each hand.
    // This deforms using two blended nonelastic.
	// This does a single RK4 integration per pair of keyframes.
	//
    // You probably don't want to do this on the CPU like this. Instead of accumulating
    // keyframe poses and applying them all at once, instead do a single Runge-Kutta
    // integration for that frame. That will be much faster than this-- even cheaper on the GPU
    // than the adaptive integrators. However, this achieves the same results.
    if (true)
    {
        Mesh mesh = readmesh("data\\meshes\\test1_mesh.bin");
        Stroke strokerighthand = readstroke("data\\strokes\\test1_righthandstroke.bin");
        Stroke strokelefthand = readstroke("data\\strokes\\test1_lefthandstroke.bin");

        strokerighthand.poses = fixFlips(strokerighthand.poses);
        strokelefthand.poses = fixFlips(strokelefthand.poses);
        DataFromPoses datarighthand = buildDataFromPoses(strokerighthand);
        DataFromPoses datalefthand = buildDataFromPoses(strokelefthand);

        for (uint i = 0; i < mesh.vertices.size(); i++)
        {
			float falloff0 = calcFalloff(mesh.vertices[i], datarighthand.deformations[0].origin, strokerighthand.innerRadius, strokerighthand.outerRadius);
			float falloff1 = calcFalloff(mesh.vertices[i], datalefthand.deformations[0].origin, strokelefthand.innerRadius, strokelefthand.outerRadius);

			if (falloff0 > 0 || falloff1 > 0)
			{
				float starttime = datarighthand.deformations[0].time;
				float endtime = datarighthand.deformations[datarighthand.deformations.size() - 1].time + datarighthand.deformations[datarighthand.deformations.size() - 1].dt;

				float minfalloff = min(falloff0, falloff1);

				float t = starttime;
				float maxt = lerp(starttime, endtime, minfalloff);
				int frame = 0;
				while (t < maxt)
				{
					float dt = datarighthand.deformations[frame].dt;
					dt = min(maxt - t, dt);

					mesh.vertices[i] = IntegrateNonElasticTwoDeformers_RungeKutta(mesh.vertices[i], t, t + dt, datarighthand.deformations[frame], datarighthand.deformations[frame]);

					frame++;
					t += dt;
				}

				// when the falloff of one deformer is less than the falloff of the second deformer,
				// then we need to integrate from the minfalloff to maxfalloff using the deformer
				// with the larger falloff
				frame = max(frame - 1, 0);
				float maxfalloff = max(falloff0, falloff1);
				maxt = lerp(starttime, endtime, maxfalloff);
				const vector<deformation::Deformation>& deformations = falloff0 > falloff1 ? datarighthand.deformations : datalefthand.deformations;
				while (t < maxt)
				{
					float dt = deformations[frame].dt;
					dt = min(maxt - t, dt);

					deformation::Deformation deformation = deformations[frame];
					deformation.origin = deformation.origin + deformation.linearVelocity*(t - deformation.time);

					mesh.vertices[i] = IntegrateNonElastic_RungeKutta(mesh.vertices[i], t, t + dt, deformation);

					frame++;
					t += dt;
				}
			}
		}

        writeobj("data\\testresult5.obj", mesh);
        printf("test5 success\n");
    }

    // --------------------
    // This demonstrates deforming with Kelvinlets.
    // It uses a single stream of poses from the right hand in Medium.
    // This deforms with a single Kelvinlet.
	// This does a single RK4 integration per pair of keyframes.
	//
    // You probably don't want to do this on the CPU like this. Instead of accumulating
    // keyframe poses and applying them all at once, instead do a single Runge-Kutta
    // integration for that frame. That will be much faster than this-- even cheaper on the GPU
    // than the adaptive integrators. However, this achieves the same results.
    if (true)
    {
        Mesh mesh = readmesh("data\\meshes\\test0_mesh.bin");
        Stroke stroke = readstroke("data\\strokes\\test0_righthandstroke.bin");

        stroke.poses = fixFlips(stroke.poses);
        DataFromPoses data = buildDataFromPoses(stroke);

        for (uint i = 0; i < mesh.vertices.size(); i++)
        {
			for (int frame = 0; frame < data.kelvinlets.size(); frame++)
			{
				mesh.vertices[i] = IntegrateKelvinlets_RungeKutta(mesh.vertices[i], data.kelvinlets[frame].time, data.kelvinlets[frame].time + data.kelvinlets[frame].dt, data.kelvinlets[frame]);
			}
		}

        writeobj("data\\testresult6.obj", mesh);
        printf("test6 success\n");
    }

    // --------------------
    // This demonstrates deforming with multiple Kelvinlets.
    // Two streams of poses were recorded in Medium, one for each hand.
    // This deforms with two Kelvinlets.
	// This does a single RK4 integration per pair of keyframes.
	//
    // You probably don't want to do this on the CPU like this. Instead of accumulating
    // keyframe poses and applying them all at once, instead do a single Runge-Kutta
    // integration for that frame. That will be much faster than this-- even cheaper on the GPU
    // than the adaptive integrators. However, this achieves the same results.
    if (true)
    {
        Mesh mesh = readmesh("data\\meshes\\test1_mesh.bin");
        Stroke strokerighthand = readstroke("data\\strokes\\test1_righthandstroke.bin");
        Stroke strokelefthand = readstroke("data\\strokes\\test1_lefthandstroke.bin");

        strokerighthand.poses = fixFlips(strokerighthand.poses);
        strokelefthand.poses = fixFlips(strokelefthand.poses);
        DataFromPoses datarighthand = buildDataFromPoses(strokerighthand);
        DataFromPoses datalefthand = buildDataFromPoses(strokelefthand);

        for (uint i = 0; i < mesh.vertices.size(); i++)
        {
			for (int frame = 0; frame < datarighthand.kelvinlets.size(); frame++)
			{
				mesh.vertices[i] = IntegrateKelvinletsTwoDeformers_RungeKutta(mesh.vertices[i], datarighthand.kelvinlets[frame].time, datarighthand.kelvinlets[frame].time + datarighthand.kelvinlets[frame].dt, datarighthand.kelvinlets[frame], datalefthand.kelvinlets[frame]);
			}
		}

        writeobj("data\\testresult7.obj", mesh);
        printf("test7 success\n");
    }

    printf("All tests successfully completed\n");

    return 0;
}
