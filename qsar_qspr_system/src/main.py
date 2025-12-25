import argparse
import os
from agents.ingestion import DataIngestionAgent
from agents.preprocessing import PreprocessingAgent
from agents.modeling import ModelBuilderAgent
from agents.evaluation import EvaluationAgent

def main():
    parser = argparse.ArgumentParser(description="QSAR/QSPR Multi-Agent System")
    parser.add_argument('--data', type=str, required=True, help='Path to the input CSV file')
    parser.add_argument('--target', type=str, required=True, help='Name of the target column')
    args = parser.parse_args()

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
    print(df_clean.head())

    # 3. Model Building
    print("\n--- Starting Model Building ---")
    model_agent = ModelBuilderAgent(target_col=args.target)
    model, X_test, y_test = model_agent.train_model(df_clean)

    # 4. Evaluation
    print("\n--- Starting Evaluation ---")
    eval_agent = EvaluationAgent()
    metrics = eval_agent.evaluate_model(model, X_test, y_test)

    # Save model (optional)
    model_path = "model.pkl"
    model_agent.save_model(model_path)

if __name__ == "__main__":
    main()
