import numpy as np

class GridAgent:
    def __init__(self, resolution=1.0):
        self.resolution = resolution

    def generate_grid(self, coordinates, features, padding=2.0):
        """
        Generates a 3D grid representation of features.
        coordinates: list of (x, y, z) tuples
        features: list of feature values (e.g., charge, hydrophobicity) corresponding to coordinates
        """
        coords = np.array(coordinates)
        if coords.size == 0:
            return None

        min_coords = coords.min(axis=0) - padding
        max_coords = coords.max(axis=0) + padding

        dims = np.ceil((max_coords - min_coords) / self.resolution).astype(int)
        grid = np.zeros(dims)

        for coord, feat in zip(coords, features):
            grid_idx = ((coord - min_coords) / self.resolution).astype(int)
            # Simple occupancy/feature assignment.
            # In a real system, we might use Gaussian splatting.
            if 0 <= grid_idx[0] < dims[0] and \
               0 <= grid_idx[1] < dims[1] and \
               0 <= grid_idx[2] < dims[2]:
                grid[tuple(grid_idx)] += feat

        return {
            'grid': grid.tolist(), # Convert to list for JSON serialization
            'origin': min_coords.tolist(),
            'dimensions': dims.tolist(),
            'resolution': self.resolution
        }
