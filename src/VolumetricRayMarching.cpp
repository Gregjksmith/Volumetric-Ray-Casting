#include "VolumetricRayMarching.h"

gs::Point::Point()
{
	pos[0] = 0.0;
	pos[1] = 0.0;
	pos[2] = 0.0;
}
gs::Point::Point(float x, float y, float z)
{
	pos[0] = x;
	pos[1] = y;
	pos[2] = z;
}
void gs::Point::clear()
{
	pos[0] = 0.0;
	pos[1] = 0.0;
	pos[2] = 0.0;
}

gs::VolumePoint::VolumePoint()
{
	pos[0] = 0.0;
	pos[1] = 0.0;
	pos[2] = 0.0;

	color[0] = 0.0;
	color[1] = 0.0;
	color[2] = 0.0;

	absorption[0] = 0.0;
	absorption[1] = 0.0;
	absorption[2] = 0.0;
}

gs::VolumePoint::VolumePoint(float* p, float* c, float* a) : Point(p[0],p[1],p[2])
{
	absorption[0] = a[0];
	absorption[1] = a[1];
	absorption[2] = a[2];

	color[0] = c[0];
	color[1] = c[1];
	color[2] = c[2];
}

gs::VolumePoint::~VolumePoint()
{

}

gs::Light::Light()
{

}
gs::Light::Light(float* p, float intensity) : Point(p[0], p[1], p[2])
{
	this->intensity = intensity;
}
gs::Light::Light(float pX, float pY, float pZ, float intensity) : Point(pX,pY,pZ)
{
	this->intensity = intensity;
}

gs::Light::~Light()
{

}

gs::Point::~Point()
{

}

gs::Bounds::Bounds()
{

}
gs::Bounds::~Bounds()
{

}


gs::RayHitInfo::RayHitInfo()
{
	clear();
}
gs::RayHitInfo::~RayHitInfo()
{

}
gs::RayHitInfo::RayHitInfo(RayHitInfo &rhi)
{
	this->ray = rhi.ray;
	this->max = rhi.max;
	this->min = rhi.min;

	this->tmin = rhi.tmin;
	this->tmax = rhi.tmax;
	this->txmin = rhi.txmin;
	this->txmax = rhi.txmax;
	this->tymin = rhi.tymin;
	this->tymax = rhi.tymax;
	this->tzmin = rhi.tzmin;
	this->tzmax = rhi.tzmax;
}
void gs::RayHitInfo::clear()
{
	this->ray.clear();
	this->max.clear();
	this->min.clear();

	this->tmin = 0.0;
	this->tmax = std::numeric_limits<float>::max();
	this->txmin = 0.0;
	this->txmax = std::numeric_limits<float>::max();
	this->tymin = 0.0;
	this->tymax = std::numeric_limits<float>::max();
	this->tzmin = 0.0;
	this->tzmax = std::numeric_limits<float>::max();
}
void gs::RayHitInfo::computeMinMax(Ray* r)
{
	ray.dir.pos[0] = r->dir.x();
	ray.dir.pos[1] = r->dir.y();
	ray.dir.pos[2] = r->dir.z();

	ray.pos.pos[0] = r->pos.x();
	ray.pos.pos[1] = r->pos.y();
	ray.pos.pos[2] = r->pos.z();

	min.pos[0] = ray.pos.x() + tmin*ray.dir.x();
	min.pos[1] = ray.pos.y() + tmin*ray.dir.y();
	min.pos[2] = ray.pos.z() + tmin*ray.dir.z();

	max.pos[0] = ray.pos.x() + tmax*ray.dir.x();
	max.pos[1] = ray.pos.y() + tmax*ray.dir.y();
	max.pos[2] = ray.pos.z() + tmax*ray.dir.z();
}


gs::Ray::Ray()
{
	clear();
}
gs::Ray::Ray(float* pos, float* dir)
{
	set(pos, dir);
}
gs::Ray::Ray(Ray& rhs)
{
	this->pos = rhs.pos;
	this->dir = rhs.dir;
}
gs::Ray::~Ray()
{

}
void gs::Ray::set(float* pos, float* dir)
{
	this->pos.pos[0] = pos[0];
	this->pos.pos[1] = pos[1];
	this->pos.pos[2] = pos[2];

	this->dir.pos[0] = dir[0];
	this->dir.pos[1] = dir[1];
	this->dir.pos[2] = dir[2];
}
void gs::Ray::clear()
{
	this->pos.pos[0] = 0.0;
	this->pos.pos[1] = 0.0;
	this->pos.pos[2] = 0.0;
}

gs::VolumetricImageParams::VolumetricImageParams()
{
	screenWidth = 800; // default values
	screenHeight = 600;
	rayMarchIterations = 50;
}
gs::VolumetricImageParams::~VolumetricImageParams()
{

}

void gs::computeBounds(std::vector<gs::VolumePoint*> &volumePoints, gs::Bounds* bounds)
{
	float minX, minY, minZ;
	float maxX, maxY, maxZ;

	minX = volumePoints[0]->pos[0];
	minY = volumePoints[0]->pos[1];
	minZ = volumePoints[0]->pos[2];
	maxX = minX;
	maxY = minY;
	maxZ = minZ;

	for (int i = 1; i < volumePoints.size(); i++)
	{
		minX = MIN(minX, volumePoints[i]->pos[0]);
		minY = MIN(minY, volumePoints[i]->pos[1]);
		minZ = MIN(minZ, volumePoints[i]->pos[2]);
		
		maxX = MAX(maxX, volumePoints[i]->pos[0]);
		maxY = MAX(maxY, volumePoints[i]->pos[1]);
		maxZ = MAX(maxZ, volumePoints[i]->pos[2]);
	}

	bounds->min.pos[0] = minX;
	bounds->min.pos[1] = minY;
	bounds->min.pos[2] = minZ;

	bounds->max.pos[0] = maxX;
	bounds->max.pos[1] = maxY;
	bounds->max.pos[2] = maxZ;
}

bool gs::intersect(Ray* ray, Bounds* bounds, RayHitInfo* hitInfo)
{
	float tmin = 0.0f;
	float tmax = std::numeric_limits<float>::max();
	float txmin = 0, tymin = 0, tzmin = 0;
	float txmax = std::numeric_limits<float>::max();
	float tymax = std::numeric_limits<float>::max();
	float tzmax = std::numeric_limits<float>::max();

	if (abs(ray->dir.x()) < PARALLEL_THRESHOLD)
	{
		// ray is parallel to slab.
		if (ray->pos.x() < bounds->min.x() || ray->pos.x() > bounds->max.x())
			return false;
	}
	if (abs(ray->dir.y()) < PARALLEL_THRESHOLD)
	{
		// ray is parallel to slab.
		if (ray->pos.y() < bounds->min.y() || ray->pos.y() > bounds->max.y())
			return false;
	}
	if (abs(ray->dir.z()) < PARALLEL_THRESHOLD)
	{
		// ray is parallel to slab.
		if (ray->pos.z() < bounds->min.x() || ray->pos.z() > bounds->max.z())
			return false;
	}

	txmin = (bounds->min.x() - ray->pos.x()) / ray->dir.x();
	txmax = (bounds->max.x() - ray->pos.x()) / ray->dir.x();
	if (txmin > txmax)
		swap(&txmin, &txmax);

	tymin = (bounds->min.y() - ray->pos.y()) / ray->dir.y();
	tymax = (bounds->max.y() - ray->pos.y()) / ray->dir.y();
	if (tymin > tymax)
		swap(&tymin, &tymax);


	tzmin = (bounds->min.z() - ray->pos.z()) / ray->dir.z();
	tzmax = (bounds->max.z() - ray->pos.z()) / ray->dir.z();
	if (tzmin > tzmax)
		swap(&tzmin, &tzmax);

	tmax = MIN(txmax, tymax);
	tmax = MIN(tzmax, tmax);
	tmin = MAX(txmin, tymin);
	tmin = MAX(tzmin, tmin);

	if (tmin > tmax)
		return false;

	hitInfo->tmax = tmax;
	hitInfo->tmin = tmin;
	hitInfo->txmax = txmax;
	hitInfo->txmin = txmin;
	hitInfo->tymax = tymax;
	hitInfo->tymin = tymin;
	hitInfo->tzmax = tzmax;
	hitInfo->tzmin = tzmin;
	hitInfo->computeMinMax(ray);
	return true;
}

float gs::length(float* x)
{
	float l = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
	return sqrt(l);
}

void gs::normalize(float* dir)
{
	float l = length(dir);
	if (l == 0.0)
		return;

	dir[0] = dir[0] / l;
	dir[1] = dir[1] / l;
	dir[2] = dir[2] / l;
}

void gs::rayMarchVolume(std::vector<gs::VolumePoint*> &volumePoints, std::vector<gs::Light*> lights, VolumetricImageParams& imageParams, cv::Mat& volumetricImage)
{
	if (volumePoints.size() <= 0)
		return;

	cv::Mat points(volumePoints.size(), 3, CV_32F);
	for (int i = 0; i < volumePoints.size(); i++)
	{
		points.at<float>(i, 0) = volumePoints[i]->pos[0];
		points.at<float>(i, 1) = volumePoints[i]->pos[1];
		points.at<float>(i, 2) = volumePoints[i]->pos[2];
	}
	cv::flann::Index* tree = new cv::flann::Index(points, cv::flann::KDTreeIndexParams(), cvflann::FLANN_DIST_EUCLIDEAN);

	Bounds bounds;
	computeBounds(volumePoints, &bounds);
	
	std::vector<float> reducedRadiance;
	Ray lightRay;
	lightRay.pos.pos[0] = 0.0;
	lightRay.pos.pos[1] = 10.0;
	lightRay.pos.pos[2] = 0.0;

	RayHitInfo rhi;

	for (int i = 0; i < volumePoints.size(); i++)
	{
		float radianceSample = 0.0f;
		for (int lightIndex = 0; lightIndex < lights.size(); lightIndex++)
		{
			lightRay.pos.pos[0] = lights[lightIndex]->pos[0];
			lightRay.pos.pos[1] = lights[lightIndex]->pos[1];
			lightRay.pos.pos[2] = lights[lightIndex]->pos[2];

			lightRay.dir.pos[0] = volumePoints[i]->pos[0] - lights[lightIndex]->pos[0];
			lightRay.dir.pos[1] = volumePoints[i]->pos[1] - lights[lightIndex]->pos[1];
			lightRay.dir.pos[2] = volumePoints[i]->pos[2] - lights[lightIndex]->pos[2];
			normalize(lightRay.dir.pos);

			if (intersect(&lightRay, &bounds, &rhi))
			{
				float lightRay[] = { volumePoints[i]->pos[0] - rhi.min.pos[0], volumePoints[i]->pos[1] - rhi.min.pos[1] , volumePoints[i]->pos[2] - rhi.min.pos[2] };
				float l = length(lightRay);
				float a = pow(volumePoints[i]->absorption[0], 2.0) + pow(volumePoints[i]->absorption[1], 2.0) + pow(volumePoints[i]->absorption[2], 2.0);
				a = sqrt(a);
				radianceSample += lights[lightIndex]->intensity*exp(-l*a);
			}
		}
		reducedRadiance.push_back(radianceSample);
	}

	int screenWidth = imageParams.screenWidth;
	int screenHeight = imageParams.screenHeight;
	int maxDepthMarchSteps = imageParams.rayMarchIterations;
	const int K = 10;

	float depthMarchStep = pow(bounds.max.x() - bounds.min.x(), 2.0) + pow(bounds.max.y() - bounds.min.y(), 2.0) + pow(bounds.max.z() - bounds.min.z(), 2.0);
	depthMarchStep = sqrt(depthMarchStep);
	depthMarchStep = depthMarchStep / (float)maxDepthMarchSteps;

	int imageSize = screenWidth*screenHeight * 4;

	cv::Mat image = cv::Mat::zeros(screenHeight, screenWidth, CV_8UC4);
	for (int i = 0; i < image.rows; i++)
	{
		for (int j = 0; j < image.cols; j++)
		{
			image.at<cv::Vec4b>(i, j) = cv::Vec4b(0, 0, 0, 255);
		}
	}

	cv::Mat ind;
	cv::Mat dist;
	cv::Mat rayMarchPosition(1,3,CV_32F);

	Ray ray;

	int screenWidthHalf = screenWidth / 2;
	int screenHeightHalf = screenHeight / 2;
	int pixelIndex = 0;
	for (int pixelY = 0; pixelY < screenWidth; pixelY++)
	{
		float x = (float)(pixelY - screenWidthHalf);
		x = x / screenWidthHalf;

		for (int pixelX = 0; pixelX < screenHeight; pixelX++)
		{
			float y = (float)(pixelX - screenHeightHalf);
			y = y / screenHeightHalf;

			float dir[] = { x, -y, -1.0 };
			normalize(dir);
			float pos[] = { 0.0,0.0,0.0 };
			ray.set(pos, dir);

			float color[] = { 0.0, 0.0, 0.0 };

			if (intersect(&ray, &bounds, &rhi))
			{
				float transport[] = { 1.0, 1.0, 1.0 };

				float depth = rhi.tmax - rhi.tmin;
				
				int totalSteps = ceil(depth / depthMarchStep);
				float step = depth / (float)totalSteps;

				for (int rayMarchIndex = 0; rayMarchIndex < totalSteps; rayMarchIndex++)
				{
					rayMarchPosition.at<float>(0, 0) = rhi.min.x() + dir[0] * rayMarchIndex*step;
					rayMarchPosition.at<float>(0, 1) = rhi.min.y() + dir[1] * rayMarchIndex*step;
					rayMarchPosition.at<float>(0, 2) = rhi.min.z() + dir[2] * rayMarchIndex*step;

					tree->knnSearch(rayMarchPosition, ind, dist, K);
					float absorp[3];
					int index = ind.at<int>(0,0);
					
					absorp[0] = volumePoints[index]->absorption[0];
					absorp[1] = volumePoints[index]->absorption[1];
					absorp[2] = volumePoints[index]->absorption[2];

					transport[0] = transport[0] * exp(-absorp[0] * step);
					transport[1] = transport[1] * exp(-absorp[1] * step);
					transport[2] = transport[2] * exp(-absorp[2] * step);

					float mediaRad[] = { 0.0,0.0,0.0 };

					for (int kIndex = 0; kIndex < K; kIndex++)
					{
						// 1 row, K cols.
						int index = ind.at<int>(0, kIndex);
						float d = dist.at<float>(0, kIndex);
						d = sqrt(d);

						mediaRad[0] += reducedRadiance[index] * exp(-absorp[0] * d);
						mediaRad[1] += reducedRadiance[index] * exp(-absorp[1] * d);
						mediaRad[2] += reducedRadiance[index] * exp(-absorp[2] * d);
					}
					color[0] += transport[0] * mediaRad[0] / (float)K;
					color[1] += transport[1] * mediaRad[1] / (float)K;
					color[2] += transport[2] * mediaRad[2] / (float)K;

				}

				cv::Vec4b sample;
				
				color[0] = MIN(color[0], 1.0);
				color[1] = MIN(color[1], 1.0);
				color[2] = MIN(color[2], 1.0);

				sample[0] = round(255 * color[0]);
				sample[1] = round(255 * color[1]);
				sample[2] = round(255 * color[2]);
				sample[3] = 255;
				image.at<cv::Vec4b>(pixelX, pixelY) = sample;
			}
			pixelIndex++;
		}
	}
	volumetricImage = image;
}