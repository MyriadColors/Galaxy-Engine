#pragma once

#include "Particles/particle.h"
#include "Physics/physics.h"

struct UpdateVariables;

struct CollisionGrid {
private:
	// Memory-efficient grid storage: flat vector instead of nested vectors
	std::vector<size_t> cellData;
	std::vector<size_t> cellOffsets;
	
	// Cache particle radii to avoid repeated calculations
	std::vector<float> cachedRadii;
	
	// Grid state caching
	int lastCellAmountX = 0;
	int lastCellAmountY = 0;
	float lastCellSize = 0.0f;
	bool needsCellSizeRecalculation = true;
	bool needsRadiiRecalculation = true;
	
	// Pre-calculated neighbor offsets for better cache performance
	std::vector<int> neighborOffsets;
	
	void updateCachedRadii(const std::vector<ParticleRendering>& rParticles, 
	                      const UpdateVariables& myVar);
	void initializeNeighborOffsets(int cellAmountX);

public:
	float cellSize = 0.0f;

	void buildGrid(std::vector<ParticlePhysics>& pParticles, std::vector<ParticleRendering>& rParticles,
		Physics& physics, UpdateVariables& myVar, glm::vec2& gridSize, float& dt);
	
	void invalidateCellSize() {
		needsCellSizeRecalculation = true;
	}
	
	void invalidateRadii() {
		needsRadiiRecalculation = true;
	}
};