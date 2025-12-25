import argparse
import os
import json
from agents.ingestion import DataIngestionAgent
from agents.preprocessing import PreprocessingAgent
from agents.modeling import ModelBuilderAgent
from agents.evaluation import EvaluationAgent
from agents.pharmacophore import PharmacophoreAgent
from agents.activesite import ActiveSiteAgent
from agents.grid import GridAgent

def main():
    parser = argparse.ArgumentParser(description="QSAR/QSPR Multi-Agent System")
    subparsers = parser.add_subparsers(dest='command', help='Sub-commands')

    # QSAR Command
    qsar_parser = subparsers.add_parser('qsar', help='Run QSAR pipeline')
    qsar_parser.add_argument('--data', type=str, required=True, help='Path to the input CSV file')
    qsar_parser.add_argument('--target', type=str, required=True, help='Name of the target column')

    # Pharmacophore Command
    pharm_parser = subparsers.add_parser('pharmacophore', help='Generate Pharmacophore')
    pharm_parser.add_argument('--smiles', type=str, required=True, help='SMILES string')
    pharm_parser.add_argument('--mode', type=str, choices=['2d', '3d'], default='2d')

    # Active Site Command
    site_parser = subparsers.add_parser('activesite', help='Analyze Active Site')
    site_parser.add_argument('--pdb', type=str, required=True, help='Path to PDB file')
    site_parser.add_argument('--ligand', type=str, help='Ligand Residue Name')

    # Grid Command
    grid_parser = subparsers.add_parser('grid', help='Generate Grid')
    grid_parser.add_argument('--coords', type=str, help='JSON string of coordinates list')

    args = parser.parse_args()

    if args.command == 'qsar':
        run_qsar(args)
    elif args.command == 'pharmacophore':
        run_pharmacophore(args)
    elif args.command == 'activesite':
        run_activesite(args)
    elif args.command == 'grid':
        run_grid(args)
    else:
        parser.print_help()

def run_qsar(args):
    # 1. Data Ingestion
    print("--- Starting Data Ingestion ---")
    ingestion_agent = DataIngestionAgent(args.data)
    df = ingestion_agent.load_data()
    if df is None:
        return

    # 2. Preprocessing
    print("\n--- Starting Preprocessing ---")
    preprocessing_agent = PreprocessingAgent()
    df_processed = preprocessing_agent.calculate_descriptors(df)
    df_clean = preprocessing_agent.clean_data(df_processed)
    print(f"Data processed. Shape: {df_clean.shape}")

    # 3. Model Building
    print("\n--- Starting Model Building ---")
    model_agent = ModelBuilderAgent(target_col=args.target)
    model, X_test, y_test = model_agent.train_model(df_clean)

    # 4. Evaluation
    print("\n--- Starting Evaluation ---")
    eval_agent = EvaluationAgent()
    metrics = eval_agent.evaluate_model(model, X_test, y_test)

    # Save model
    model_path = "model.pkl"
    model_agent.save_model(model_path)

def run_pharmacophore(args):
    print(f"--- Generating {args.mode.upper()} Pharmacophore ---")
    agent = PharmacophoreAgent()
    if args.mode == '2d':
        # Note: 2D generation doesn't give positions in space usually, RDKit features are 3D based.
        # But we can list families.
        feats = agent.generate_pharmacophore(args.smiles)
    else:
        feats = agent.generate_3d_pharmacophore(args.smiles)

    print(json.dumps(feats, indent=2))

def run_activesite(args):
    print(f"--- Analyzing Active Site in {args.pdb} ---")
    agent = ActiveSiteAgent()
    residues = agent.extract_pocket(args.pdb, ligand_resname=args.ligand)
    print(f"Found {len(residues)} pocket residues.")
    print(json.dumps(residues, indent=2))

def run_grid(args):
    print("--- Generating Grid ---")
    # Mock data if not provided
    if not args.coords:
        coords = [[0,0,0], [1,1,1], [2,0,0]]
        features = [1, 1, 0.5]
    else:
        coords = json.loads(args.coords)
        features = [1] * len(coords) # Dummy features

    agent = GridAgent(resolution=0.5)
    grid_data = agent.generate_grid(coords, features)
    print(f"Grid Dimensions: {grid_data['dimensions']}")

if __name__ == "__main__":
    main()
