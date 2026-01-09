from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np

class PreprocessingAgent:
    def __init__(self):
        pass

    def calculate_descriptors(self, df, smiles_col='SMILES', method='basic'):
        """
        Calculates molecular descriptors for a given dataframe with SMILES.

        Args:
            df (pd.DataFrame): Input dataframe containing SMILES strings.
            smiles_col (str): Column name for SMILES strings.
            method (str): Descriptor calculation method ('basic', 'extended', 'fingerprints').

        Returns:
            pd.DataFrame: Dataframe with calculated descriptors.
        """
        descriptors = []
        valid_indices = []

        for index, row in df.iterrows():
            smiles = row[smiles_col]
            mol = Chem.MolFromSmiles(smiles)

            if mol:
                desc = {}
                if method == 'basic':
                    desc = {
                        'MolWt': Descriptors.MolWt(mol),
                        'LogP': Descriptors.MolLogP(mol),
                        'NumHDonors': Descriptors.NumHDonors(mol),
                        'NumHAcceptors': Descriptors.NumHAcceptors(mol),
                        'TPSA': Descriptors.TPSA(mol)
                    }
                elif method == 'extended':
                    desc = {
                        'MolWt': Descriptors.MolWt(mol),
                        'LogP': Descriptors.MolLogP(mol),
                        'NumHDonors': Descriptors.NumHDonors(mol),
                        'NumHAcceptors': Descriptors.NumHAcceptors(mol),
                        'TPSA': Descriptors.TPSA(mol),
                        'BertzCT': Descriptors.BertzCT(mol),
                        'MolMR': Descriptors.MolMR(mol),
                        'HeavyAtomCount': Descriptors.HeavyAtomCount(mol),
                        'RotatableBonds': Descriptors.NumRotatableBonds(mol),
                        'FractionCSP3': Descriptors.FractionCSP3(mol),
                        'HallKierAlpha': Descriptors.HallKierAlpha(mol),
                        'Kappa1': Descriptors.Kappa1(mol)
                    }
                elif method == 'fingerprints':
                    # Radius 2 (ECFP4), 2048 bits
                    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                    desc = {f'fp_{i}': int(bit) for i, bit in enumerate(fp)}

                descriptors.append(desc)
                valid_indices.append(index)
            else:
                print(f"Invalid SMILES at index {index}: {smiles}")

        desc_df = pd.DataFrame(descriptors, index=valid_indices)
        result_df = df.loc[valid_indices].join(desc_df)

        return result_df

    def clean_data(self, df):
        """Handle missing values, etc."""
        # Simple dropna for now
        return df.dropna()
