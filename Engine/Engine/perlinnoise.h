#pragma once
#include <cmath>

class perlinnoise
{
public:
	perlinnoise();
	~perlinnoise();
	float GetOffset(float x, float y, float z);
	double octavePerlin(double x, double y, double z, int octaves, double persistence, perlinnoise pn);
};

