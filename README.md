# Volumetric Ray Casting

Renders a translucent point cloud using a variant of volumetric ray casting. Each point (particle) in the point clound is given a 
spectral absorption coefficient which determines the media's opacity.

Radiance entering each pixel towards the eye is calculated by ray marching through the point cloud volume and gathering radiance along the ray. Direct radiance at each point is calculated using Beer-Lambert's law (exponential decay). Media radiance is calculated using local inscattering.

## API

``
void rayMarchVolume(std::vector<gs::VolumePoint*> &volumePoints, std::vector<gs::Light*> lights, VolumetricImageParams& imageParams, cv::Mat& volumetricImage);

volumePoints : media point cloud
light : vector of point lights
imageParams : struct of image paramters. Image width, image height and ray march iterations.
volumetricImage : image rendered here.
``

## Example

``
#include <stdio.h>
#include <vector>
#include "VolumetricRayMarching.h"

void createVolumePointCloud(std::vector<gs::VolumePoint*> &pointCloud)
{
	const int numSamplePoint = 50;
	const float samplingDist = 1.0 / (float)numSamplePoint;
	
	float c[3];
	float pos[3];
	float absorption[3];

	for (float x = 0.0; x <= 1.0; x = x + samplingDist)
	{
		for (float y = 0.0; y <= 1.0; y = y + samplingDist)
		{
			for (float z = 0.0; z <= 1.0; z = z + samplingDist)
			{
				pos[0] = x - 0.5;
				pos[1] = y - 0.8;
				pos[2] = z - 2.2;

				if (x > 0.25 && y > 0.25 && z > 0.25 && x < 0.75 && y < 0.75 && z < 0.75)
				{
					absorption[0] = 9.0f;
					absorption[1] = 9.0f;
					absorption[2] = 27.0f;
				}
				else
				{
					absorption[0] = 2.0f;
					absorption[1] = 2.0f;
					absorption[2] = 2.0f;
				}

				gs::VolumePoint* vp = new gs::VolumePoint(pos, c, absorption);
				pointCloud.push_back(vp);
			}
		}
	}
}

void main()
{
	std::vector<gs::VolumePoint*> pointCloud;
	createVolumePointCloud(pointCloud);
	
	gs::Light* light = new gs::Light(0.0f, 10.0f, 0.0f, 18.0f / 255.0f);
	std::vector<gs::Light*> lightVec;
	lightVec.push_back(light);
	
	gs::VolumetricImageParams vip;
	vip.rayMarchIterations = 100;
	cv::Mat volumetricImage;

	rayMarchVolume(pointCloud, lightVec, vip, volumetricImage);
	cv::imwrite(std::string("testImage.png"), volumetricImage);
}
``