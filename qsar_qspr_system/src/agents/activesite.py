from Bio.PDB import PDBParser, Selection
import warnings

class ActiveSiteAgent:
    def __init__(self):
        self.parser = PDBParser(QUIET=True)

    def extract_pocket(self, pdb_file, ligand_resname=None, radius=5.0):
        """
        Extracts residues within a certain radius of a ligand.
        If ligand_resname is not provided, tries to find a HETATM that looks like a ligand.
        """
        try:
            structure = self.parser.get_structure('protein', pdb_file)
        except Exception as e:
            print(f"Error parsing PDB: {e}")
            return []

        model = structure[0]
        atoms = Selection.unfold_entities(model, 'A')

        ligand_atoms = []
        protein_atoms = []

        for atom in atoms:
            res = atom.get_parent()
            # Simple heuristic: HETATM is ligand unless water
            if res.id[0].startswith('H_'):
                if res.resname != 'HOH':
                    if ligand_resname:
                        if res.resname == ligand_resname:
                            ligand_atoms.append(atom)
                    else:
                        # Take the first non-water hetatm chain as ligand for simplicity
                        # In a real app, we'd be more selective
                        ligand_atoms.append(atom)
            else:
                protein_atoms.append(atom)

        if not ligand_atoms:
            print("No ligand found.")
            return []

        # Find neighbors
        # Using simple distance check for MVP. Bio.PDB has NeighborSearch but this is quick enough for small proteins.
        pocket_residues = set()
        for latom in ligand_atoms:
            for patom in protein_atoms:
                if (latom - patom) < radius:
                    pocket_residues.add(patom.get_parent())

        result = []
        for res in pocket_residues:
            result.append({
                'Residue': res.resname,
                'ID': res.id[1],
                'Chain': res.get_parent().id
            })

        return result
