from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd
import numpy as np

class PreprocessingAgent:
    def __init__(self):
        pass

    def calculate_descriptors(self, df, smiles_col='SMILES'):
        """Calculates molecular descriptors for a given dataframe with SMILES."""
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

        # Drop the original SMILES column if needed, but keeping it for reference is often good
        # result_df = result_df.drop(columns=[smiles_col])

        return result_df

    def clean_data(self, df):
        """Handle missing values, etc."""
        # Simple dropna for now
        return df.dropna()
