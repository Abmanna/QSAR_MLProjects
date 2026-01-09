from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, DataStructs, rdFingerprintGenerator
import pandas as pd
import numpy as np

class PreprocessingAgent:
    def __init__(self):
        pass

    def calculate_descriptors(self, df, smiles_col='SMILES', descriptor_type='basic'):
        """Calculates molecular descriptors for a given dataframe with SMILES."""
        if descriptor_type == 'fingerprints':
            return self.calculate_morgan_fingerprints(df, smiles_col)

        descriptors = []
        valid_indices = []

        for index, row in df.iterrows():
            smiles = row[smiles_col]
            mol = Chem.MolFromSmiles(smiles)

            if mol:
                desc = {
                    'MolWt': Descriptors.MolWt(mol),
                    'LogP': Descriptors.MolLogP(mol),
                    'NumHDonors': Descriptors.NumHDonors(mol),
                    'NumHAcceptors': Descriptors.NumHAcceptors(mol),
                    'TPSA': Descriptors.TPSA(mol)
                }
                descriptors.append(desc)
                valid_indices.append(index)
            else:
                print(f"Invalid SMILES at index {index}: {smiles}")

        desc_df = pd.DataFrame(descriptors, index=valid_indices)
        result_df = df.loc[valid_indices].join(desc_df)

        return result_df

    def calculate_morgan_fingerprints(self, df, smiles_col='SMILES', radius=2, nBits=2048):
        """Calculates Morgan Fingerprints (ECFP)."""
        fingerprints = []
        valid_indices = []

        mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=nBits)

        for index, row in df.iterrows():
            smiles = row[smiles_col]
            mol = Chem.MolFromSmiles(smiles)

            if mol:
                fp = mfpgen.GetFingerprint(mol)
                # Correctly allocate array of size nBits
                arr = np.zeros((nBits,), dtype=np.int8)
                DataStructs.ConvertToNumpyArray(fp, arr)
                fingerprints.append(arr)
                valid_indices.append(index)
            else:
                print(f"Invalid SMILES at index {index}: {smiles}")

        # Create columns for each bit
        col_names = [f'FP_{i}' for i in range(nBits)]
        fp_df = pd.DataFrame(fingerprints, index=valid_indices, columns=col_names)

        result_df = df.loc[valid_indices].join(fp_df)
        return result_df

    def clean_data(self, df):
        """Handle missing values, etc."""
        # Simple dropna for now
        return df.dropna()
