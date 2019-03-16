// Copyright(c) Facebook, Inc. and its affiliates.
// All rights reserved.
// 
// This source code is licensed under the BSD - style license found in the
// LICENSE file in the root directory of this source tree.

struct vec3
{
    vec3() {};
    vec3(float a) : x(a), y(a), z(a) {};
    vec3(float a, float b, float c) : x(a), y(b), z(c) {};

    float x, y, z;
};

vec3 operator*(vec3 a, float b)
{
    return vec3(a.x * b, a.y * b, a.z * b);
}

vec3 operator*(float a, vec3 b)
{
    return vec3(b.x * a, b.y * a, b.z * a);
}

vec3 operator/(vec3 a, float b)
{
    return vec3(a.x / b, a.y / b, a.z / b);
}

vec3 operator+(vec3 a, vec3 b)
{
    return vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}

vec3 operator-(vec3 a, vec3 b)
{
    return vec3(a.x - b.x, a.y - b.y, a.z - b.z);
}

float dot(vec3 a, vec3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

vec3 cross(vec3 a, vec3 b)
{
    return vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

float length(vec3 v)
{
    return sqrt(dot(v, v));
}

vec3 normalize(vec3 a)
{
    return a / length(a);
}

float distance(vec3 a, vec3 b)
{
    return length(a - b);
}

float saturate(float t)
{
    if (t < 0)
        t = 0;
    if (t > 1)
        t = 1;
    return t;
}

float smoothstep(float x)
{
    return x * x * (3.f - 2.f * x);
}

float smoothstep(float min, float max, float t)
{
    return smoothstep(saturate((t - min) / (max - min)));
}

float lerp(float a, float b, float t)
{
    return a + (b - a) * t;
}

struct quat
{
    quat() {};

    quat(float a, float b, float c, float d) : r(a), i(b), j(c), k(d) {};

    float r, i, j, k;

    // Returns the imaginary rotation axis of the quaternion as a vector
    vec3 v() const { return vec3(i, j, k); }

    // Converts the quaternion to an axis-angle representation
    void toAxisAngle(vec3& out_axis, float& out_angle) const
    {
        vec3 a = v();
        float len = length(a);

        out_axis = a / len;
        out_angle = atan2(len, r) * 2.0f;
    }

    quat& operator*=(float rhs)
    {
        r *= rhs;
        i *= rhs;
        j *= rhs;
        k *= rhs;
        return *this;
    };
};

quat operator*(float a, quat q)
{
    return quat(a * q.r, a * q.i, a * q.j, a * q.k);
}

quat operator*(quat q, float b)
{
    return quat(q.r * b, q.i * b, q.j * b, q.k * b);
}

quat operator*(quat a, quat b)
{
    return quat(
        a.r * b.r - a.i * b.i - a.j * b.j - a.k * b.k, 
        a.r * b.i + a.i * b.r + a.j * b.k - a.k * b.j,
        a.r * b.j - a.i * b.k + a.j * b.r + a.k * b.i, 
        a.r * b.k + a.i * b.j - a.j * b.i + a.k * b.r);
}

float dot(quat a, quat b)
{
    return (a.r * b.r + a.i * b.i + a.j * b.j + a.k * b.k);
}

quat conj(quat q)
{
    return quat(+q.r, -q.i, -q.j, -q.k);
}

float rcp(float f)
{
    return 1 / f;
}

quat rcp(quat q)
{
    return conj(q) * rcp(q.r * q.r + q.i * q.i + q.j * q.j + q.k * q.k);
}

quat inverse(quat q)
{
    return rcp(q);
}

// uses column vectors and column major ordering
// to make a combined transform of matrix W, then V,
// then finally P:
// WVP = P * V * W;
// to transform point with transform T:
// v1 = T * v0;
struct mat3x3
{
    mat3x3() {};

    // from column vectors
    mat3x3(vec3 vx, vec3 vy, vec3 vz) : cx(vx), cy(vy), cz(vz) {};

    vec3 cx, cy, cz;
};

mat3x3 operator*(float a, mat3x3 b)
{
    return mat3x3(a * b.cx, a * b.cy, a * b.cz);
}

mat3x3 operator*(mat3x3 a, float b)
{
    return mat3x3(b * a.cx, b * a.cy, b * a.cz);
}

mat3x3 operator+(mat3x3 a, mat3x3 b)
{
    return mat3x3(a.cx + b.cx, a.cy + b.cy, a.cz + b.cz);
}

vec3 operator*(mat3x3 a, vec3 b)
{
    return b.x * a.cx + b.y * a.cy + b.z * a.cz;
}
