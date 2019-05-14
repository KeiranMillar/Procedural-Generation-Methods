////////////////////////////////////////////////////////////////////////////////
// Filename: racetrack.h
////////////////////////////////////////////////////////////////////////////////
#ifndef _racetrack_H_
#define _racetrack_H_


//////////////
// INCLUDES //
//////////////
#include <d3d11.h>
#include <d3dx10math.h>
#include <stdio.h>
#include <vector>


////////////////////////////////////////////////////////////////////////////////
// Class name: racetrack
////////////////////////////////////////////////////////////////////////////////
class racetrack
{
private:
	struct VertexType
	{
		D3DXVECTOR3 position;
		D3DXVECTOR3 normal;
	};

	struct HeightMapType
	{
		float x, y, z;
		float nx, ny, nz;
	};

	struct VectorType
	{
		float x, y, z;
	};

public:
	racetrack();
	~racetrack();

	bool Initialize(ID3D11Device*, char*);
	bool InitializeTrack(ID3D11Device*, int terrainWidth, int terrainHeight);
	void GenerateRandomPoints(ID3D11Device*);
	void SortPoints(ID3D11Device* device);
	void ComplexHull(ID3D11Device* device);
	void ComplexHull2(ID3D11Device* device);
	void SeperatePoints(ID3D11Device* device);
	void AddMidpoints(ID3D11Device* device);
	void FixAngles(ID3D11Device* device);
	void Shutdown();
	void DoWork(ID3D11Device* device, bool keydown);
	void Render(ID3D11DeviceContext*);
	bool GenerateHeightMap(ID3D11Device* device, bool keydown);
	bool Smooth(ID3D11Device* device, bool keydown);
	bool AddNoise(ID3D11Device* device, bool keydown);
	float GetNoise(float x, float y, float z);
	int  GetIndexCount();

	//int pointCount = rand() % 10 + 11;
	int pointCount = 20;
	std::vector<D3DXVECTOR2> points;

private:
	bool LoadHeightMap(char*);
	void NormalizeHeightMap();
	bool CalculateNormals();
	void ShutdownHeightMap();

	bool InitializeBuffers(ID3D11Device*);
	void ShutdownBuffers();
	void RenderBuffers(ID3D11DeviceContext*);

private:
	bool m_terrainGeneratedToggle;
	int m_terrainWidth, m_terrainHeight;
	int m_vertexCount, m_indexCount;
	ID3D11Buffer *m_vertexBuffer, *m_indexBuffer;
	HeightMapType* m_heightMap;
};

#endif