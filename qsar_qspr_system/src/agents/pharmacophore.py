from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig
import os

class PharmacophoreAgent:
    def __init__(self):
        # Load default RDKit feature definitions
        fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
        self.featFactory = ChemicalFeatures.BuildFeatureFactory(fdefName)

    def generate_pharmacophore(self, smiles):
        """Generates pharmacophore features for a given SMILES string. (2D only)"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            print(f"Invalid SMILES: {smiles}")
            return []

        features = self.featFactory.GetFeaturesForMol(mol)

        pharma_features = []
        for f in features:
            # For 2D, we cannot get positions safely if no conformer exists.
            # GetPos throws exception if no conformer.
            pharma_features.append({
                'Family': f.GetFamily(),
                'Type': f.GetType(),
                'AtomIds': f.GetAtomIds(),
                'Position': None
            })

        return pharma_features

    def generate_3d_pharmacophore(self, smiles):
        """Generates 3D pharmacophore features."""
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return []

        mol = Chem.AddHs(mol)
        try:
            # Generate conformer
            params = AllChem.ETKDG()
            AllChem.EmbedMolecule(mol, params)
        except ValueError:
             print("Could not embed molecule")
             return []

        if mol.GetNumConformers() == 0:
            return []

        features = self.featFactory.GetFeaturesForMol(mol)
        pharma_features = []
        for f in features:
            pharma_features.append({
                'Family': f.GetFamily(),
                'Type': f.GetType(),
                'Position': list(f.GetPos())
            })
        return pharma_features
