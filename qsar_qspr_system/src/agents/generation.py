from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd

class GenerationAgent:
    def __init__(self):
        self.rxn_smarts = [
            # Aromatic C-H to Methyl
            ('[cH&$(c(c)c):1]>>[c:1](C)', 'Methylation'),
            # Aromatic C-H to F
            ('[cH&$(c(c)c):1]>>[c:1](F)', 'Fluorination'),
            # Aromatic C-H to Cl
            ('[cH&$(c(c)c):1]>>[c:1](Cl)', 'Chlorination'),
            # Aromatic C-H to OH
            ('[cH&$(c(c)c):1]>>[c:1](O)', 'Hydroxylation'),
            # Aliphatic C-H to Methyl
            ('[#6&!H0:1]>>[#6:1](C)', 'Aliphatic Methylation'),
            # Aliphatic C-H to F
            ('[#6&!H0:1]>>[#6:1](F)', 'Aliphatic Fluorination'),
            # Acid to Amide
            ('C(=O)[OH]>>C(=O)N', 'Amidation'),
            # Amide to Acid
            ('C(=O)N>>C(=O)[OH]', 'Hydrolysis')
        ]
        self.rxns = [(AllChem.ReactionFromSmarts(x[0]), x[1]) for x in self.rxn_smarts]

    def generate_analogs(self, smiles_list, num_analogs=50):
        """
        Generates analogs for a list of SMILES strings.

        Args:
            smiles_list (list): List of SMILES strings to modify.
            num_analogs (int): Approximate number of analogs to generate per input molecule.

        Returns:
            pd.DataFrame: DataFrame containing new molecules (SMILES) and their source.
        """
        generated_mols = []

        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                continue

            # Try to apply each transformation
            mol_analogs = set()

            # Since reaction application is deterministic for a specific set of reactants,
            # we run all defined reactions.

            for rxn, name in self.rxns:
                try:
                    # Run reaction
                    products = rxn.RunReactants((mol,))
                    for product_tuple in products:
                        for product in product_tuple:
                            try:
                                Chem.SanitizeMol(product)
                                new_smiles = Chem.MolToSmiles(product)
                                if new_smiles != smiles:
                                    mol_analogs.add((new_smiles, name, smiles))
                            except:
                                pass
                except Exception as e:
                    # Reaction failed, skip
                    pass

            # Add to main list
            generated_mols.extend(list(mol_analogs))

        # Limit if too many, but for now return all unique ones
        # If we wanted strictly `num_analogs` per molecule, we would sample.
        # But since we use deterministic reactions, we just return what we found.

        df = pd.DataFrame(generated_mols, columns=['SMILES', 'Transformation', 'Source_SMILES'])
        df = df.drop_duplicates(subset=['SMILES'])

        # If we have too many, we could sample, but usually comprehensive enumeration of 1-step analogs is better.
        if len(df) > num_analogs * len(smiles_list):
            df = df.sample(n=num_analogs * len(smiles_list), random_state=42)

        return df
