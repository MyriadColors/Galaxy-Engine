#include "pch.h"
#include "Physics/collisionGrid.h"
#include "parameters.h"

void CollisionGrid::updateCachedRadii(const std::vector<ParticleRendering>& rParticles, 
                                     const UpdateVariables& myVar) {
	cachedRadii.resize(rParticles.size());
	
	const float radiusMultiplier = myVar.particleSizeMultiplier * myVar.particleTextureHalfSize;
	
	for (size_t i = 0; i < rParticles.size(); ++i) {
		cachedRadii[i] = rParticles[i].size * radiusMultiplier;
	}
	needsRadiiRecalculation = false;
}

void CollisionGrid::initializeNeighborOffsets(int cellAmountX) {
	neighborOffsets.clear();
	neighborOffsets.reserve(9); // Maximum 9 neighbors (3x3 grid)
	
	for (int dx = -1; dx <= 1; ++dx) {
		for (int dy = -1; dy <= 1; ++dy) {
			neighborOffsets.push_back(dx + dy * cellAmountX);
		}
	}
}

void CollisionGrid::buildGrid(std::vector<ParticlePhysics>& pParticles, std::vector<ParticleRendering>& rParticles,
	Physics& physics, UpdateVariables& myVar, glm::vec2& gridSize, float& dt) {

	// OPTIMIZATION: Only recalculate cell size when needed
	if (needsCellSizeRecalculation) {
		float maxSize = 0.0f;
		for (size_t i = 0; i < rParticles.size(); i++) {
			float thisSize = rParticles[i].totalRadius * 4.0f;
			maxSize = std::max(maxSize, thisSize);
		}
		cellSize = maxSize;
		lastCellSize = cellSize;
		needsCellSizeRecalculation = false;
	} else {
		cellSize = lastCellSize;
	}

	int cellAmountX = static_cast<int>(gridSize.x / cellSize);
	int cellAmountY = static_cast<int>(gridSize.y / cellSize);
	int totalCells = cellAmountX * cellAmountY;

	// OPTIMIZATION: Update cached radii only when needed
	if (needsRadiiRecalculation || cachedRadii.size() != rParticles.size()) {
		updateCachedRadii(rParticles, myVar);
	}

	// OPTIMIZATION: Initialize neighbor offsets if grid dimensions changed
	if (lastCellAmountX != cellAmountX || neighborOffsets.empty()) {
		initializeNeighborOffsets(cellAmountX);
		lastCellAmountX = cellAmountX;
		lastCellAmountY = cellAmountY;
	}

	// OPTIMIZATION: Reuse memory with flat storage structure
	cellOffsets.assign(totalCells + 1, 0);
	
	// Count particles per cell
	for (size_t i = 0; i < pParticles.size(); ++i) {
		int xIdx = static_cast<int>(pParticles[i].pos.x / cellSize);
		int yIdx = static_cast<int>(pParticles[i].pos.y / cellSize);

		if (xIdx >= 0 && xIdx < cellAmountX && yIdx >= 0 && yIdx < cellAmountY) {
			int cellId = xIdx + yIdx * cellAmountX;
			cellOffsets[cellId + 1]++;
		}
	}

	// Convert counts to offsets (prefix sum)
	for (int i = 1; i <= totalCells; ++i) {
		cellOffsets[i] += cellOffsets[i - 1];
	}

	// Resize cellData and fill with particle indices
	cellData.resize(cellOffsets[totalCells]);
	std::vector<size_t> currentOffsets = cellOffsets; // Working copy

	for (size_t i = 0; i < pParticles.size(); ++i) {
		int xIdx = static_cast<int>(pParticles[i].pos.x / cellSize);
		int yIdx = static_cast<int>(pParticles[i].pos.y / cellSize);

		if (xIdx >= 0 && xIdx < cellAmountX && yIdx >= 0 && yIdx < cellAmountY) {
			int cellId = xIdx + yIdx * cellAmountX;
			cellData[currentOffsets[cellId]++] = i;
		}
	}

	// OPTIMIZATION: Simplified collision checking lambda with cached radii
	auto checkCollision = [&](size_t a, size_t b) {
		if (a == b) return;

		const float sumR = cachedRadii[a] + cachedRadii[b];
		const glm::vec2 delta = pParticles[a].pos - pParticles[b].pos;
		const float distSq = delta.x * delta.x + delta.y * delta.y;

		if (distSq < sumR * sumR) {
			physics.collisions(pParticles[a], pParticles[b],
			                  rParticles[a], rParticles[b], myVar, dt);
		}
	};

	// OPTIMIZATION: Parallel collision detection with improved load balancing
#pragma omp parallel for schedule(dynamic, 64)
	for (int cellId = 0; cellId < totalCells; ++cellId) {
		const size_t cellStart = cellOffsets[cellId];
		const size_t cellEnd = cellOffsets[cellId + 1];
		
		if (cellStart == cellEnd) continue; // Skip empty cells

		const int cellX = cellId % cellAmountX;
		const int cellY = cellId / cellAmountX;

		// Check collisions within the same cell
		for (size_t i = cellStart; i < cellEnd; ++i) {
			for (size_t j = i + 1; j < cellEnd; ++j) {
				checkCollision(cellData[i], cellData[j]);
			}
		}

		// Check collisions with neighboring cells (only forward neighbors to avoid duplicates)
		for (int dx = 0; dx <= 1; ++dx) {
			for (int dy = (dx == 0) ? 1 : -1; dy <= 1; ++dy) {
				const int nx = cellX + dx;
				const int ny = cellY + dy;
				
				if (nx >= cellAmountX || ny < 0 || ny >= cellAmountY) continue;

				const int neighborId = nx + ny * cellAmountX;
				const size_t neighborStart = cellOffsets[neighborId];
				const size_t neighborEnd = cellOffsets[neighborId + 1];

				for (size_t i = cellStart; i < cellEnd; ++i) {
					for (size_t j = neighborStart; j < neighborEnd; ++j) {
						checkCollision(cellData[i], cellData[j]);
					}
				}
			}
		}
	}
}
