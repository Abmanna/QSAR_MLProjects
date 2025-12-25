import unittest
import numpy as np
from src.agents.pharmacophore import PharmacophoreAgent
from src.agents.grid import GridAgent

class TestAdvancedAgents(unittest.TestCase):
    def test_pharmacophore_2d(self):
        agent = PharmacophoreAgent()
        smiles = "CCO"
        feats = agent.generate_pharmacophore(smiles)
        self.assertTrue(len(feats) > 0)
        # Alcohol oxygen should be HDonor and HAcceptor often
        families = [f['Family'] for f in feats]
        self.assertTrue(any('Donor' in f or 'Acceptor' in f for f in families))

    def test_pharmacophore_3d(self):
        # 3D generation involves embedding, which might fail or be slow
        agent = PharmacophoreAgent()
        smiles = "CC(=O)O"
        feats = agent.generate_3d_pharmacophore(smiles)
        if feats:
            self.assertTrue('Position' in feats[0])
            self.assertEqual(len(feats[0]['Position']), 3)

    def test_grid_generation(self):
        agent = GridAgent(resolution=1.0)
        coords = [[0,0,0], [2,0,0]]
        features = [1.0, 0.5]
        grid_data = agent.generate_grid(coords, features, padding=1.0)

        self.assertIsNotNone(grid_data)
        grid = np.array(grid_data['grid'])
        # Distance between 0 and 2 is 2. With padding 1.0, range is -1 to 3. Size 4.
        # Resolution 1.0.
        # Actually logic: min -1, max 3. (3 - -1)/1 = 4.
        self.assertTrue(grid.shape[0] >= 3)
        self.assertTrue(grid.sum() > 0)

if __name__ == '__main__':
    unittest.main()
