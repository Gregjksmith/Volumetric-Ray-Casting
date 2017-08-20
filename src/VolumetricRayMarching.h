#pragma once

/*
Volumetric Ray Marching

Greg Smith

2017
*/


#include <stdio.h>
#include <vector>
#include <opencv2/opencv.hpp>

#define PARALLEL_THRESHOLD 0.0000000001

namespace gs
{
	class Point;
	class VolumePoint;
	class Bounds;
	class RayHitInfo;
	class Ray;


	class Point
	{
	public:
		Point();
		Point(float x, float y, float z);
		virtual ~Point();
		inline Point(Point &rhs)
		{
			this->pos[0] = rhs.pos[0];
			this->pos[1] = rhs.pos[1];
			this->pos[2] = rhs.pos[2];
		}

		inline Point operator+(const Point& rhs)
		{
			Point r;
			r.pos[0] = pos[0] + rhs.pos[0];
			r.pos[1] = pos[1] + rhs.pos[1];
			r.pos[2] = pos[2] + rhs.pos[2];
			return r;
		}
		inline Point operator-(const Point& rhs)
		{
			Point r;
			r.pos[0] = pos[0] - rhs.pos[0];
			r.pos[1] = pos[1] - rhs.pos[1];
			r.pos[2] = pos[2] - rhs.pos[2];
			return r;
		}
		inline Point& operator=(const Point& rhs)
		{
			this->pos[0] = rhs.pos[0];
			this->pos[1] = rhs.pos[1];
			this->pos[2] = rhs.pos[2];

			return *this;
		}

		void clear();

		inline float x()
		{
			return pos[0];
		}
		inline float y()
		{
			return pos[1];
		}
		inline float z()
		{
			return pos[2];
		}
		float pos[3];
	};

	class VolumePoint : public Point
	{
	public:
		VolumePoint();
		VolumePoint(float* p, float* c, float* a);
		virtual ~VolumePoint();

		float color[3];
		float absorption[3];
	};

	class Light : public Point
	{
	public:
		Light();
		Light(float* p, float intensity);
		Light(float pX, float pY, float pZ, float intensity);
		virtual ~Light();

		float intensity;
	};

	class Bounds
	{
	public:
		Bounds();
		virtual ~Bounds();
		Point min;
		Point max;
	};

	class Ray
	{
	public:
		Ray();
		Ray(float* pos, float* dir);
		Ray(Ray& rhs);
		virtual ~Ray();
		void set(float* pos, float* dir);
		void clear();

		Point pos;
		Point dir;
	};

	class RayHitInfo
	{
	public:
		RayHitInfo();
		RayHitInfo(RayHitInfo &rhi);
		virtual ~RayHitInfo();
		
		void clear();
		Ray ray;
		Point max;
		Point min;
		float tmin;
		float tmax;

		float txmin;
		float txmax;
		float tymin;
		float tymax;
		float tzmin;
		float tzmax;

		void computeMinMax(Ray* r);
	};

	class VolumetricImageParams
	{
	public:
		VolumetricImageParams();
		virtual ~VolumetricImageParams();

		int screenWidth;
		int screenHeight;
		int rayMarchIterations;
	};

	inline void swap(float* x, float* y)
	{
		float temp = *x;
		*x = *y;
		*y = temp;
	}

	float length(float* x);
	void normalize(float* x);
	bool intersect(Ray* ray, Bounds* bounds, RayHitInfo* hitInfo);
	void computeBounds(std::vector<gs::VolumePoint*> &volumePoints, Bounds* bounds);
	void rayMarchVolume(std::vector<gs::VolumePoint*> &volumePoints, std::vector<gs::Light*> lights, VolumetricImageParams& imageParams, cv::Mat& volumetricImage);
}