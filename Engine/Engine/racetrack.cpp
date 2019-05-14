////////////////////////////////////////////////////////////////////////////////
// Filename: racetrack.cpp
////////////////////////////////////////////////////////////////////////////////
#include "racetrack.h"
#include "perlinnoise.h"
#include <cmath>
#include <vector>
#include <math.h>
#include <ctime>

using namespace std;
#define PI 3.14159265

void Sort()
{
	
}

float sgn(float x)
{
	if (x > 0) return 1;
	if (x < 0) return -1;
	return 0;
}

racetrack::racetrack()
{
	m_vertexBuffer = 0;
	m_indexBuffer = 0;
	m_heightMap = 0;
	m_terrainGeneratedToggle = false;
}


racetrack::~racetrack()
{
}


bool racetrack::InitializeTrack(ID3D11Device* device, int terrainWidth, int terrainHeight)
{
	int index;
	float height = 0.0;
	bool result;

	// Save the dimensions of the terrain.
	m_terrainWidth = terrainWidth;
	m_terrainHeight = terrainHeight;

	// Create the structure to hold the terrain data.
	m_heightMap = new HeightMapType[m_terrainWidth * m_terrainHeight];
	if (!m_heightMap)
	{
		return false;
	}

	// Initialise the data in the height map (flat).
	for (int j = 0; j<m_terrainHeight; j++)
	{
		for (int i = 0; i<m_terrainWidth; i++)
		{
			
			index = (m_terrainHeight * j) + i;

			m_heightMap[index].x = (float)i;
			m_heightMap[index].y = 0.0f;
			m_heightMap[index].z = (float)j;

		}
	}

	DoWork(device, true);

	for (int i = 0; i < points.size(); i++)
	{
		index = points[i].x + (points[i].y * m_terrainWidth);

		m_heightMap[index].y = 10.0f * i;
	}

	

	//even though we are generating a flat terrain, we still need to normalise it. 
	// Calculate the normals for the terrain data.
	result = CalculateNormals();
	if (!result)
	{
		return false;
	}

	// Initialize the vertex and index buffer that hold the geometry for the terrain.
	result = InitializeBuffers(device);
	if (!result)
	{
		return false;
	}

	return true;
}


bool racetrack::Initialize(ID3D11Device* device, char* heightMapFilename)
{
	bool result;


	// Load in the height map for the terrain.
	result = LoadHeightMap(heightMapFilename);
	if (!result)
	{
		return false;
	}

	// Normalize the height of the height map.
	NormalizeHeightMap();

	// Calculate the normals for the terrain data.
	result = CalculateNormals();
	if (!result)
	{
		return false;
	}

	// Initialize the vertex and index buffer that hold the geometry for the terrain.
	result = InitializeBuffers(device);
	if (!result)
	{
		return false;
	}

	return true;
}


void racetrack::GenerateRandomPoints(ID3D11Device* device)
{
	int boundaryFraction = 5;

	for (int i = 0; i <= pointCount; i++)
	{
		float x = (rand() % (((m_terrainWidth / boundaryFraction)*(boundaryFraction-1))-(m_terrainWidth / boundaryFraction))) + (m_terrainWidth / boundaryFraction);
		float y = (rand() % (((m_terrainHeight / boundaryFraction)*(boundaryFraction - 1)) - (m_terrainHeight / boundaryFraction))) + (m_terrainHeight / boundaryFraction);
		points.push_back(D3DXVECTOR2(x, y));
	}
}


void racetrack::SortPoints(ID3D11Device* device)
{
	//Sort my points into counter-clockwise order
	D3DXVECTOR2 temp;
	for (int i = 0; i < points.size(); i++)
	{
		for (int j = 1; j < points.size(); j++)
		{

			if (points[j].x < points[j - 1].x)
			{
				temp = points[j];
				points[j] = points[j - 1];
				points[j - 1] = temp;
			}

			else if (points[j].x == points[j - 1].x)
			{
				if (points[j].y < points[j - 1].y)
				{
					temp = points[j];
					points[j] = points[j - 1];
					points[j - 1] = temp;
				}
			}
		}
	}
}


int CrossProd(const D3DXVECTOR2 &Orig, const D3DXVECTOR2 &A, const D3DXVECTOR2 &B)
{
	return (A.x - Orig.x) * (B.y - Orig.y) - (A.y - Orig.y) * (B.x - Orig.x);
}


// A utility function to return square of distance
// between p1 and p2
int distSq(D3DXVECTOR2 p1, D3DXVECTOR2 p2)
{
	return (p1.x - p2.x)*(p1.x - p2.x) +
		(p1.y - p2.y)*(p1.y - p2.y);
}


// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise
int orientation(D3DXVECTOR2 p, D3DXVECTOR2 q, D3DXVECTOR2 r)
{
	int val = (q.y - p.y) * (r.x - q.x) -
		(q.x - p.x) * (r.y - q.y);

	if (val == 0) return 0;  // colinear
	return (val > 0) ? 1 : 2; // clock or counterclock wise
}


void racetrack::ComplexHull(ID3D11Device* device)
{
	// Init Vars  
	float W = points.size(), Z = 0;
	std::vector<D3DXVECTOR2> H(pointCount * 2);

	// Build lower hull  
	for (int i = 0; i < W; i++) 
	{
		while (Z >= 2 && CrossProd(H[Z - 2], H[Z - 1], points[i]) <= 0)
			Z--;
		H[Z++] = points[i];
	}

	// Build upper hull  
	for (int i = W - 2, t = Z + 1; i >= 0; i--) 
	{
		while (Z >= t && CrossProd(H[Z - 2], H[Z - 1], points[i]) <= 0)
			Z--;
		H[Z++] = points[i];
	}

	H.resize(Z);
	H.erase(H.begin());
	points.erase(points.begin(), points.end());
	points.resize(Z);
	points = H;
}


void racetrack::ComplexHull2(ID3D11Device* device)
{
	// Initialize Result
	vector<D3DXVECTOR2> hull;

	// Find the leftmost point
	int l = 0;
	for (int i = 1; i < points.size(); i++)
		if (points[i].x < points[l].x)
			l = i;

	// Start from leftmost point, keep moving counterclockwise
	// until reach the start point again.  This loop runs O(h)
	// times where h is number of points in result or output.
	int p = l, q;
	do
	{
		// Add current point to result
		hull.push_back(points[p]);

		// Search for a point 'q' such that orientation(p, x,
		// q) is counterclockwise for all points 'x'. The idea
		// is to keep track of last visited most counterclock-
		// wise point in q. If any point 'i' is more counterclock-
		// wise than q, then update q.
		q = (p + 1) % points.size();
		for (int i = 0; i < points.size(); i++)
		{
			// If i is more counterclockwise than current q, then
			// update q
			if (orientation(points[p], points[i], points[q]) == 2)
				q = i;
		}

		// Now q is the most counterclockwise with respect to p
		// Set p as q for next iteration, so that q is added to
		// result 'hull'
		p = q;

	} while (p != l);  // While we don't come to first point

	points.erase(points.begin(), points.end());
	points.resize(hull.size());
	points = hull;
}


void racetrack::SeperatePoints(ID3D11Device* device)
{
	float dst = 3; //I found that 5 is a good value, though maybe, depending on your scale you'll need other value.  
	float dst2 = dst*dst;

	for (int i = 0; i < points.size(); i++)
	{
		int j = i + 1;
		if (j == points.size())
		{
			j = 0;
		}
		if ( ((points[j].x - points[i].x) * (points[j].x - points[i].x)) + ((points[j].y - points[i].y) * (points[j].y - points[i].y)) < dst2)
		{
			float hx = points[j].x - points[i].x;
			float hy = points[j].y - points[i].y;
			float hl = (float)std::sqrt((hx*hx) + (hy*hy));
			hx /= hl;
			hy /= hl;
			float dif = dst - hl;
			hx *= dif;
			hy *= dif;
			points[j].x += hx;
			points[j].y += hy;
			points[i].x -= hx;
			points[i].y -= hy;
		}
	}
}


//void racetrack::AddMidpoints(ID3D11Device* device)
//{
//	D3DXVECTOR2 disp;
//	std::vector<D3DXVECTOR2> rSet(points.size() *2 + 1);
//	float difficulty = 1.0f; //the closer the value is to 0, the harder the track should be. Grows exponentially.  
//	float maxDisp = 0.0f; // Again, this may change to fit your units.  
//	for (int i = 0; i < points.size(); i++)
//	{
//		float dispLen = (float)pow((rand()%10 / 10), difficulty) * maxDisp;
//		disp.x = 0.0f;
//		disp.y = 1.0f;
//		float randomAngle = (float)((rand() % 10 / 10) * 360);
//		disp.x *= sin(randomAngle*PI / 180);
//		disp.y *= cos(randomAngle*PI / 180);
//		disp *= dispLen;
//		rSet[i * 2] = points[i];
//		if (i == points.size() - 1)
//		{
//			rSet[i * 2 + 1] = ((points[i] + points[0]) / 2);
//		}
//		else
//		{
//			rSet[i * 2 + 1] = ((points[i] + points[i + 1]) / 2);
//		}
//		//Explaining: a mid point can be found with (dataSet[i]+dataSet[i+1])/2.  
//		//Then we just add the displacement.  
//	}
//	rSet.pop_back();
//	points.erase(points.begin(), points.end());
//	points.resize(rSet.size());
//	points = rSet;
//	//push apart again, so we can stabilize the points distances.  
//}

void racetrack::AddMidpoints(ID3D11Device* device)
{
	std::vector<D3DXVECTOR2> rSet(points.size() * 2 + 1);
	for (int i = 0; i < points.size(); i++)
	{
		rSet[i * 2] = points[i];
		if (i == points.size() - 1)
		{
			rSet[i * 2 + 1] = ((points[i] + points[0]) / 2);
		}
		else
		{
			rSet[i * 2 + 1] = ((points[i] + points[i + 1]) / 2);
		}
		//Explaining: a mid point can be found with (dataSet[i]+dataSet[i+1])/2.  
		//Then we just add the displacement.  
	}
	rSet.pop_back();
	points.erase(points.begin(), points.end());
	points.resize(rSet.size());
	points = rSet;
	//push apart again, so we can stabilize the points distances.  
}

void racetrack::FixAngles(ID3D11Device* device)
{
	for (int i = 0; i < points.size(); i++)
	{
		int previous = (i - 1);
		if (i == 0)
		{
			previous = (points.size() - 1);
		}
		int next = (i + 1);
		if (i == points.size() - 1)
		{
			next = 0;
		}
		float px = points[i].x - points[previous].x;
		float py = points[i].y - points[previous].y;
		float pl = (float)sqrt(px*px + py*py);
		px /= pl;
		py /= pl;

		float nx = points[i].x - points[next].x;
		float ny = points[i].y - points[next].y;
		nx = -nx;
		ny = -ny;
		float nl = (float)sqrt(nx*nx + ny*ny);
		nx /= nl;
		ny /= nl;
		//I got a vector going to the next and to the previous points, normalised.  

		float a = (float)atan2(px * ny - py * nx, px * nx + py * ny); // perp dot product between the previous and next point.

		if (abs((a * PI) / 180) <= 100)
		{
			float nA = ((100 * sgn(a)) * PI) / 180;
			float diff = nA - a;
			float cosA = (float)cos(diff);
			float sinA = (float)sin(diff);
			float newX = nx * cosA - ny * sinA;
			float newY = nx * sinA + ny * cosA;
			newX *= nl;
			newY *= nl;
			points[next].x = points[i].x + newX;
			points[next].y = points[i].y + newY;
		}
	}
}


void racetrack::Shutdown()
{
	// Release the vertex and index buffer.
	ShutdownBuffers();

	// Release the height map data.
	ShutdownHeightMap();

	return;
}


void racetrack::Render(ID3D11DeviceContext* deviceContext)
{
	// Put the vertex and index buffers on the graphics pipeline to prepare them for drawing.
	RenderBuffers(deviceContext);

	return;
}


int racetrack::GetIndexCount()
{
	return m_indexCount;
}


bool racetrack::GenerateHeightMap(ID3D11Device* device, bool keydown)
{

	bool result;
	//the toggle is just a bool that I use to make sure this is only called ONCE when you press a key
	//until you release the key and start again. We dont want to be generating the terrain 500
	//times per second. 
	if (keydown && (!m_terrainGeneratedToggle))
	{
		int index;
		float height = 0.0;


		//loop through the terrain and set the hieghts how we want. This is where we generate the terrain
		//in this case I will run a sin-wave through the terrain in one axis.

		for (int j = 0; j<m_terrainHeight; j++)
		{
			for (int i = 0; i<m_terrainWidth; i++)
			{
				index = (m_terrainHeight * j) + i;

				//Sin Wave
				//m_heightMap[index].x = (float)i;
				//m_heightMap[index].y = (float)(sin((float)j/(m_terrainWidth/12))*3.0) + (cos((float)i / (m_terrainWidth / 12))*3.0); //magic numbers ahoy, just to ramp up the height of the sin function so its visible.
				//m_heightMap[index].z = (float)j;

				//Random Height
				m_heightMap[index].x = (float)i;
				if (index % ((int)rand() % 200 + 1) == 0)
				{
					m_heightMap[index].y = (int)rand() % 20;
				}
				m_heightMap[index].z = (float)j;

				//Perlin Noise
				//m_heightMap[index].x = (float)i;
				//m_heightMap[index].y = (GetNoise(i, 0.0, j))*2;
				//m_heightMap[index].z = (float)j;
			}
		}

		result = CalculateNormals();
		if (!result)
		{
			return false;
		}

		// Initialize the vertex and index buffer that hold the geometry for the terrain.
		result = InitializeBuffers(device);
		if (!result)
		{
			return false;
		}

		m_terrainGeneratedToggle = true;
	}
	else
	{
		m_terrainGeneratedToggle = false;
	}



	return true;
}


bool racetrack::AddNoise(ID3D11Device* device, bool keydown)
{
	bool result;
	int index;

	//loop through the terrain and set offset the terrain based on perlin noise
	if (keydown)
	{
		for (int j = 0; j < m_terrainHeight; j++)
		{
			for (int i = 0; i < m_terrainWidth; i++)
			{
				index = (m_terrainHeight * j) + i;
				m_heightMap[index].x = (float)i;
				m_heightMap[index].y += (GetNoise(i, 0.0, j)) * 2;
				m_heightMap[index].z = (float)j;
			}
		}
	}

	result = CalculateNormals();
	if (!result)
	{
		return false;
	}

	// Initialize the vertex and index buffer that hold the geometry for the terrain.
	result = InitializeBuffers(device);
	if (!result)
	{
		return false;
	}

	return true;
}


bool racetrack::Smooth(ID3D11Device* device, bool keydown)
{
	int counter = 1;
	int index;
	bool result;
	if (keydown)
	{
		for (int j = 0; j < m_terrainHeight; j++)
		{
			for (int i = 0; i < m_terrainWidth; i++)
			{
				index = (m_terrainHeight * j) + i;

				if (index - 129 >= 0)
				{
					m_heightMap[index].y += m_heightMap[index - 129].y;
					counter++;
				}

				if (index - 128 >= 0)
				{
					m_heightMap[index].y += m_heightMap[index - 128].y;
					counter++;
				}

				if (index - 127 >= 0)
				{
					m_heightMap[index].y += m_heightMap[index - 127].y;
					counter++;
				}

				if (index - 1 >= 0)
				{
					m_heightMap[index].y += m_heightMap[index - 1].y;
					counter++;
				}

				if (index + 1 <= ((m_terrainHeight*m_terrainWidth) - 1))
				{
					m_heightMap[index].y += m_heightMap[index + 1].y;
					counter++;
				}

				if (index + 127 <= ((m_terrainHeight*m_terrainWidth) - 1))
				{
					m_heightMap[index].y += m_heightMap[index + 127].y;
					counter++;
				}

				if (index + 128 <= ((m_terrainHeight*m_terrainWidth) - 1))
				{
					m_heightMap[index].y += m_heightMap[index + 128].y;
					counter++;
				}

				if (index + 129 <= ((m_terrainHeight*m_terrainWidth) - 1))
				{
					m_heightMap[index].y += m_heightMap[index + 129].y;
					counter++;
				}

				m_heightMap[index].y /= counter;
				counter = 1;
			}
		}

		result = CalculateNormals();
		if (!result)
		{
			return false;
		}

		// Initialize the vertex and index buffer that hold the geometry for the terrain.
		result = InitializeBuffers(device);
		if (!result)
		{
			return false;
		}
		m_terrainGeneratedToggle = true;
	}
	else
	{
		m_terrainGeneratedToggle = false;
	}
	return true;
}


float racetrack::GetNoise(float x, float y, float z)
{
	perlinnoise offset;
	return (offset.GetOffset(x / 8, y, z / 8));
}


bool racetrack::LoadHeightMap(char* filename)
{
	FILE* filePtr;
	int error;
	unsigned int count;
	BITMAPFILEHEADER bitmapFileHeader;
	BITMAPINFOHEADER bitmapInfoHeader;
	int imageSize, i, j, k, index;
	unsigned char* bitmapImage;
	unsigned char height;


	// Open the height map file in binary.
	error = fopen_s(&filePtr, filename, "rb");
	if (error != 0)
	{
		return false;
	}

	// Read in the file header.
	count = fread(&bitmapFileHeader, sizeof(BITMAPFILEHEADER), 1, filePtr);
	if (count != 1)
	{
		return false;
	}

	// Read in the bitmap info header.
	count = fread(&bitmapInfoHeader, sizeof(BITMAPINFOHEADER), 1, filePtr);
	if (count != 1)
	{
		return false;
	}

	// Save the dimensions of the terrain.
	m_terrainWidth = bitmapInfoHeader.biWidth;
	m_terrainHeight = bitmapInfoHeader.biHeight;

	// Calculate the size of the bitmap image data.
	imageSize = m_terrainWidth * m_terrainHeight * 3;

	// Allocate memory for the bitmap image data.
	bitmapImage = new unsigned char[imageSize];
	if (!bitmapImage)
	{
		return false;
	}

	// Move to the beginning of the bitmap data.
	fseek(filePtr, bitmapFileHeader.bfOffBits, SEEK_SET);

	// Read in the bitmap image data.
	count = fread(bitmapImage, 1, imageSize, filePtr);
	if (count != imageSize)
	{
		return false;
	}

	// Close the file.
	error = fclose(filePtr);
	if (error != 0)
	{
		return false;
	}

	// Create the structure to hold the height map data.
	m_heightMap = new HeightMapType[m_terrainWidth * m_terrainHeight];
	if (!m_heightMap)
	{
		return false;
	}

	// Initialize the position in the image data buffer.
	k = 0;

	// Read the image data into the height map.
	for (j = 0; j<m_terrainHeight; j++)
	{
		for (i = 0; i<m_terrainWidth; i++)
		{
			height = bitmapImage[k];

			index = (m_terrainHeight * j) + i;

			m_heightMap[index].x = (float)i;
			m_heightMap[index].y = (float)height;
			m_heightMap[index].z = (float)j;

			k += 3;
		}
	}

	// Release the bitmap image data.
	delete[] bitmapImage;
	bitmapImage = 0;

	return true;
}


void racetrack::NormalizeHeightMap()
{
	int i, j;


	for (j = 0; j<m_terrainHeight; j++)
	{
		for (i = 0; i<m_terrainWidth; i++)
		{
			m_heightMap[(m_terrainHeight * j) + i].y /= 15.0f;
		}
	}

	return;
}


bool racetrack::CalculateNormals()
{
	int i, j, index1, index2, index3, index, count;
	float vertex1[3], vertex2[3], vertex3[3], vector1[3], vector2[3], sum[3], length;
	VectorType* normals;


	// Create a temporary array to hold the un-normalized normal vectors.
	normals = new VectorType[(m_terrainHeight - 1) * (m_terrainWidth - 1)];
	if (!normals)
	{
		return false;
	}

	// Go through all the faces in the mesh and calculate their normals.
	for (j = 0; j<(m_terrainHeight - 1); j++)
	{
		for (i = 0; i<(m_terrainWidth - 1); i++)
		{
			index1 = (j * m_terrainHeight) + i;
			index2 = (j * m_terrainHeight) + (i + 1);
			index3 = ((j + 1) * m_terrainHeight) + i;

			// Get three vertices from the face.
			vertex1[0] = m_heightMap[index1].x;
			vertex1[1] = m_heightMap[index1].y;
			vertex1[2] = m_heightMap[index1].z;

			vertex2[0] = m_heightMap[index2].x;
			vertex2[1] = m_heightMap[index2].y;
			vertex2[2] = m_heightMap[index2].z;

			vertex3[0] = m_heightMap[index3].x;
			vertex3[1] = m_heightMap[index3].y;
			vertex3[2] = m_heightMap[index3].z;

			// Calculate the two vectors for this face.
			vector1[0] = vertex1[0] - vertex3[0];
			vector1[1] = vertex1[1] - vertex3[1];
			vector1[2] = vertex1[2] - vertex3[2];
			vector2[0] = vertex3[0] - vertex2[0];
			vector2[1] = vertex3[1] - vertex2[1];
			vector2[2] = vertex3[2] - vertex2[2];

			index = (j * (m_terrainHeight - 1)) + i;

			// Calculate the cross product of those two vectors to get the un-normalized value for this face normal.
			normals[index].x = (vector1[1] * vector2[2]) - (vector1[2] * vector2[1]);
			normals[index].y = (vector1[2] * vector2[0]) - (vector1[0] * vector2[2]);
			normals[index].z = (vector1[0] * vector2[1]) - (vector1[1] * vector2[0]);
		}
	}

	// Now go through all the vertices and take an average of each face normal 	
	// that the vertex touches to get the averaged normal for that vertex.
	for (j = 0; j<m_terrainHeight; j++)
	{
		for (i = 0; i<m_terrainWidth; i++)
		{
			// Initialize the sum.
			sum[0] = 0.0f;
			sum[1] = 0.0f;
			sum[2] = 0.0f;

			// Initialize the count.
			count = 0;

			// Bottom left face.
			if (((i - 1) >= 0) && ((j - 1) >= 0))
			{
				index = ((j - 1) * (m_terrainHeight - 1)) + (i - 1);

				sum[0] += normals[index].x;
				sum[1] += normals[index].y;
				sum[2] += normals[index].z;
				count++;
			}

			// Bottom right face.
			if ((i < (m_terrainWidth - 1)) && ((j - 1) >= 0))
			{
				index = ((j - 1) * (m_terrainHeight - 1)) + i;

				sum[0] += normals[index].x;
				sum[1] += normals[index].y;
				sum[2] += normals[index].z;
				count++;
			}

			// Upper left face.
			if (((i - 1) >= 0) && (j < (m_terrainHeight - 1)))
			{
				index = (j * (m_terrainHeight - 1)) + (i - 1);

				sum[0] += normals[index].x;
				sum[1] += normals[index].y;
				sum[2] += normals[index].z;
				count++;
			}

			// Upper right face.
			if ((i < (m_terrainWidth - 1)) && (j < (m_terrainHeight - 1)))
			{
				index = (j * (m_terrainHeight - 1)) + i;

				sum[0] += normals[index].x;
				sum[1] += normals[index].y;
				sum[2] += normals[index].z;
				count++;
			}

			// Take the average of the faces touching this vertex.
			sum[0] = (sum[0] / (float)count);
			sum[1] = (sum[1] / (float)count);
			sum[2] = (sum[2] / (float)count);

			// Calculate the length of this normal.
			length = sqrt((sum[0] * sum[0]) + (sum[1] * sum[1]) + (sum[2] * sum[2]));

			// Get an index to the vertex location in the height map array.
			index = (j * m_terrainHeight) + i;

			// Normalize the final shared normal for this vertex and store it in the height map array.
			m_heightMap[index].nx = (sum[0] / length);
			m_heightMap[index].ny = (sum[1] / length);
			m_heightMap[index].nz = (sum[2] / length);
		}
	}

	// Release the temporary normals.
	delete[] normals;
	normals = 0;

	return true;
}


void racetrack::ShutdownHeightMap()
{
	if (m_heightMap)
	{
		delete[] m_heightMap;
		m_heightMap = 0;
	}

	return;
}


bool racetrack::InitializeBuffers(ID3D11Device* device)
{
	VertexType* vertices;
	unsigned long* indices;
	int index, i, j;
	D3D11_BUFFER_DESC vertexBufferDesc, indexBufferDesc;
	D3D11_SUBRESOURCE_DATA vertexData, indexData;
	HRESULT result;
	int index1, index2, index3, index4;


	// Calculate the number of vertices in the terrain mesh.
	m_vertexCount = (m_terrainWidth - 1) * (m_terrainHeight - 1) * 6;

	// Set the index count to the same as the vertex count.
	m_indexCount = m_vertexCount;

	// Create the vertex array.
	vertices = new VertexType[m_vertexCount];
	if (!vertices)
	{
		return false;
	}

	// Create the index array.
	indices = new unsigned long[m_indexCount];
	if (!indices)
	{
		return false;
	}

	// Initialize the index to the vertex buffer.
	index = 0;

	// Load the vertex and index array with the terrain data.
	for (j = 0; j<(m_terrainHeight - 1); j++)
	{
		for (i = 0; i<(m_terrainWidth - 1); i++)
		{
			if (i % 2 == 0)
			{
				if (j % 2 == 0)
				{
					index4 = (m_terrainHeight * j) + i;          // Bottom left.
					index3 = (m_terrainHeight * j) + (i + 1);      // Bottom right.
					index2 = (m_terrainHeight * (j + 1)) + i;      // Upper left.
					index1 = (m_terrainHeight * (j + 1)) + (i + 1);  // Upper right.
				}
				else
				{
					index3 = (m_terrainHeight * j) + i;          // Bottom left.
					index1 = (m_terrainHeight * j) + (i + 1);      // Bottom right.
					index4 = (m_terrainHeight * (j + 1)) + i;      // Upper left.
					index2 = (m_terrainHeight * (j + 1)) + (i + 1);  // Upper right.
				}
			}
			else
			{
				if (j % 2 == 0)
				{
					index2 = (m_terrainHeight * j) + i;          // Bottom left.
					index4 = (m_terrainHeight * j) + (i + 1);      // Bottom right.
					index1 = (m_terrainHeight * (j + 1)) + i;      // Upper left.
					index3 = (m_terrainHeight * (j + 1)) + (i + 1);  // Upper right.
				}
				else
				{
					index4 = (m_terrainHeight * j) + i;          // Bottom left.
					index3 = (m_terrainHeight * j) + (i + 1);      // Bottom right.
					index2 = (m_terrainHeight * (j + 1)) + i;      // Upper left.
					index1 = (m_terrainHeight * (j + 1)) + (i + 1);  // Upper right.
				}
			}

			// Upper left.
			vertices[index].position = D3DXVECTOR3(m_heightMap[index3].x, m_heightMap[index3].y, m_heightMap[index3].z);
			vertices[index].normal = D3DXVECTOR3(m_heightMap[index3].nx, m_heightMap[index3].ny, m_heightMap[index3].nz);
			indices[index] = index;
			index++;

			// Upper right.
			vertices[index].position = D3DXVECTOR3(m_heightMap[index4].x, m_heightMap[index4].y, m_heightMap[index4].z);
			vertices[index].normal = D3DXVECTOR3(m_heightMap[index4].nx, m_heightMap[index4].ny, m_heightMap[index4].nz);
			indices[index] = index;
			index++;

			// Bottom left.
			vertices[index].position = D3DXVECTOR3(m_heightMap[index1].x, m_heightMap[index1].y, m_heightMap[index1].z);
			vertices[index].normal = D3DXVECTOR3(m_heightMap[index1].nx, m_heightMap[index1].ny, m_heightMap[index1].nz);
			indices[index] = index;
			index++;

			// Bottom left.
			vertices[index].position = D3DXVECTOR3(m_heightMap[index1].x, m_heightMap[index1].y, m_heightMap[index1].z);
			vertices[index].normal = D3DXVECTOR3(m_heightMap[index1].nx, m_heightMap[index1].ny, m_heightMap[index1].nz);
			indices[index] = index;
			index++;

			// Upper right.
			vertices[index].position = D3DXVECTOR3(m_heightMap[index4].x, m_heightMap[index4].y, m_heightMap[index4].z);
			vertices[index].normal = D3DXVECTOR3(m_heightMap[index4].nx, m_heightMap[index4].ny, m_heightMap[index4].nz);
			indices[index] = index;
			index++;

			// Bottom right.
			vertices[index].position = D3DXVECTOR3(m_heightMap[index2].x, m_heightMap[index2].y, m_heightMap[index2].z);
			vertices[index].normal = D3DXVECTOR3(m_heightMap[index2].nx, m_heightMap[index2].ny, m_heightMap[index2].nz);
			indices[index] = index;
			index++;
		}
	}

	// Set up the description of the static vertex buffer.
	vertexBufferDesc.Usage = D3D11_USAGE_DEFAULT;
	vertexBufferDesc.ByteWidth = sizeof(VertexType) * m_vertexCount;
	vertexBufferDesc.BindFlags = D3D11_BIND_VERTEX_BUFFER;
	vertexBufferDesc.CPUAccessFlags = 0;
	vertexBufferDesc.MiscFlags = 0;
	vertexBufferDesc.StructureByteStride = 0;

	// Give the subresource structure a pointer to the vertex data.
	vertexData.pSysMem = vertices;
	vertexData.SysMemPitch = 0;
	vertexData.SysMemSlicePitch = 0;

	// Now create the vertex buffer.
	result = device->CreateBuffer(&vertexBufferDesc, &vertexData, &m_vertexBuffer);
	if (FAILED(result))
	{
		return false;
	}

	// Set up the description of the static index buffer.
	indexBufferDesc.Usage = D3D11_USAGE_DEFAULT;
	indexBufferDesc.ByteWidth = sizeof(unsigned long) * m_indexCount;
	indexBufferDesc.BindFlags = D3D11_BIND_INDEX_BUFFER;
	indexBufferDesc.CPUAccessFlags = 0;
	indexBufferDesc.MiscFlags = 0;
	indexBufferDesc.StructureByteStride = 0;

	// Give the subresource structure a pointer to the index data.
	indexData.pSysMem = indices;
	indexData.SysMemPitch = 0;
	indexData.SysMemSlicePitch = 0;

	// Create the index buffer.
	result = device->CreateBuffer(&indexBufferDesc, &indexData, &m_indexBuffer);
	if (FAILED(result))
	{
		return false;
	}

	// Release the arrays now that the buffers have been created and loaded.
	delete[] vertices;
	vertices = 0;

	delete[] indices;
	indices = 0;

	return true;
}


void racetrack::ShutdownBuffers()
{
	// Release the index buffer.
	if (m_indexBuffer)
	{
		m_indexBuffer->Release();
		m_indexBuffer = 0;
	}

	// Release the vertex buffer.
	if (m_vertexBuffer)
	{
		m_vertexBuffer->Release();
		m_vertexBuffer = 0;
	}

	return;
}


void racetrack::RenderBuffers(ID3D11DeviceContext* deviceContext)
{
	unsigned int stride;
	unsigned int offset;


	// Set vertex buffer stride and offset.
	stride = sizeof(VertexType);
	offset = 0;

	// Set the vertex buffer to active in the input assembler so it can be rendered.
	deviceContext->IASetVertexBuffers(0, 1, &m_vertexBuffer, &stride, &offset);

	// Set the index buffer to active in the input assembler so it can be rendered.
	deviceContext->IASetIndexBuffer(m_indexBuffer, DXGI_FORMAT_R32_UINT, 0);

	// Set the type of primitive that should be rendered from this vertex buffer, in this case triangles.
	deviceContext->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);

	return;
}


void racetrack::DoWork(ID3D11Device* device, bool keydown)
{
	if (keydown)
	{
		srand(time(NULL));
		GenerateRandomPoints(device);
		SortPoints(device);
		ComplexHull2(device);
		for (int i = 0; i < 3; i++)
		{
			//SeperatePoints(device);
		}
		//AddMidpoints(device);
		for (int i = 0; i < 3; i++)
		{
			//SeperatePoints(device);
		}
		//SortPoints(device);
		//ComplexHull2(device);
		//FixAngles(device);
		//D3DXVec2CatmullRom();
	}
}