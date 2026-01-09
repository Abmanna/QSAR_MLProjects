import argparse
import os
import pandas as pd
from agents.ingestion import DataIngestionAgent
from agents.preprocessing import PreprocessingAgent
from agents.modeling import ModelBuilderAgent
from agents.evaluation import EvaluationAgent
from agents.generation import GenerationAgent

def main():
    parser = argparse.ArgumentParser(description="QSAR/QSPR Multi-Agent System")
    parser.add_argument('--data', type=str, required=True, help='Path to the input CSV file')
    parser.add_argument('--target', type=str, required=True, help='Name of the target column')
    parser.add_argument('--descriptors', type=str, default='basic', choices=['basic', 'extended', 'fingerprints'], help='Type of descriptors to calculate')
    parser.add_argument('--design', action='store_true', help='Enable design/generation of new molecular entities')
    args = parser.parse_args()

    # 1. Data Ingestion
    print("--- Starting Data Ingestion ---")
    ingestion_agent = DataIngestionAgent(args.data)
    df = ingestion_agent.load_data()
    if df is None:
        return

    # 2. Preprocessing
    print(f"\n--- Starting Preprocessing (Method: {args.descriptors}) ---")
    preprocessing_agent = PreprocessingAgent()
    df_processed = preprocessing_agent.calculate_descriptors(df, method=args.descriptors)
    df_clean = preprocessing_agent.clean_data(df_processed)
    print(f"Data processed. Shape: {df_clean.shape}")
    print(df_clean.head())

    # 3. Model Building
    print("\n--- Starting Model Building ---")
    model_agent = ModelBuilderAgent(target_col=args.target)
    model, X_test, y_test = model_agent.train_model(df_clean)

    # 4. Evaluation
    print("\n--- Starting Evaluation ---")
    eval_agent = EvaluationAgent()
    metrics = eval_agent.evaluate_model(model, X_test, y_test)
    print("Metrics:", metrics)

    # Save model (optional)
    model_path = "model.pkl"
    model_agent.save_model(model_path)

    # 5. Design New Molecules (Optional)
    if args.design:
        print("\n--- Starting Molecule Design ---")
        print("Note: Selecting top active molecules (Maximized Target) as seeds.")

        # Identify top molecules to use as seeds
        # Assume target is 'activity', so we want higher or lower?
        # Let's assume user wants to MAXIMIZE the target value.

        # Sort by target
        if args.target in df_clean.columns:
            top_mols = df_clean.sort_values(by=args.target, ascending=False).head(10)
            seed_smiles = top_mols['SMILES'].tolist()
            print(f"Selected {len(seed_smiles)} top molecules as seeds.")

            # Generate Analogs
            gen_agent = GenerationAgent()
            generated_df = gen_agent.generate_analogs(seed_smiles)
            print(f"Generated {len(generated_df)} unique analogs.")

            if not generated_df.empty:
                # Calculate descriptors for new molecules
                # We need to use the SAME descriptor method
                print("Calculating descriptors for new molecules...")

                # generated_df has 'SMILES', 'Transformation', 'Source_SMILES'
                # preprocessing_agent expects 'SMILES' column

                gen_processed = preprocessing_agent.calculate_descriptors(generated_df, method=args.descriptors)
                # Cleaning might drop some if descriptors fail
                gen_clean = preprocessing_agent.clean_data(gen_processed)

                if not gen_clean.empty:
                    # Predict
                    print("Predicting activity for new molecules...")
                    predictions = model_agent.predict(gen_clean)
                    gen_clean['Predicted_' + args.target] = predictions

                    # Sort by prediction
                    best_new_mols = gen_clean.sort_values(by='Predicted_' + args.target, ascending=False).head(10)

                    print("\nTop 10 Designed Molecules:")
                    cols_to_show = ['SMILES', 'Source_SMILES', 'Transformation', 'Predicted_' + args.target]
                    print(best_new_mols[cols_to_show].to_string(index=False))

                    # Save results
                    gen_clean.to_csv("designed_molecules.csv", index=False)
                    print("\nAll designed molecules saved to designed_molecules.csv")
                else:
                    print("No valid descriptors calculated for generated molecules.")
            else:
                print("No analogs generated.")
        else:
            print("Target column not found in dataframe, skipping design.")

if __name__ == "__main__":
    main()
