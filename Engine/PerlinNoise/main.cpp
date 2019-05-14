#include "main.h"
#include "perlinnoise.h"
#include "ppm.h"

double octavePerlin(double x, double y, double z, int octaves, double persistence, perlinnoise pn)
{
	double total = 0;
	double frequency = 1;
	double amplitude = 1;
	double maxValue = 0;  // Used for normalizing result to 0.0 - 1.0
	for (int i = 0; i<octaves; i++) {
		total += pn.GetOffset(x * frequency, y * frequency, z * frequency) * amplitude;

		maxValue += amplitude;

		amplitude *= persistence;
		frequency *= 2;
	}

	return total / maxValue;
}


int main() {
	// Define the size of the image
	unsigned int width = 256, height = 256;

	// Create an empty PPM image
	ppm image(width, height);

	// Create a PerlinNoise object with a random permutation vector generated with seed
	unsigned int seed = 237;
	perlinnoise pn;

	unsigned int kk = 0;
	// Visit every pixel of the image and assign a color generated with Perlin noise
	for (unsigned int i = 0; i < height; ++i) {     // y
		for (unsigned int j = 0; j < width; ++j) {  // x
			double x = (double)j / ((double)width);
			double y = (double)i / ((double)height);

			// Typical Perlin noise
			double n = octavePerlin(2 * x, 2 * y, 0.8, 8, 0.7, pn);
			n = (n + 1) / 2;

			// Wood like structure
			//double n = 20 * pn.GetOffset(x, y, 0.8);
			//n = n - floor(n);
			//n = (n + 1) / 2;

			// Map the values to the [0, 255] interval, for simplicity we use 
			// tones of grey
			image.r[kk] = floor(255 * n);
			image.g[kk] = floor(255 * n);
			image.b[kk] = floor(255 * n);
			kk++;
		}
	}

	// Save the image in a binary PPM file
	image.write("result.ppm");

	return 0;
}