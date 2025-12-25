import unittest
import pandas as pd
import os
from src.agents.ingestion import DataIngestionAgent
from src.agents.preprocessing import PreprocessingAgent
from src.agents.modeling import ModelBuilderAgent
from src.agents.evaluation import EvaluationAgent

class TestQSARSystem(unittest.TestCase):
    def setUp(self):
        # Create a dummy CSV file
        self.csv_path = "test_data.csv"
        self.model_path = "test_model.pkl"
        self.df = pd.DataFrame({
            'SMILES': ['CCO', 'CCN', 'CC(=O)O'],
            'Activity': [1.2, 1.5, 0.8]
        })
        self.df.to_csv(self.csv_path, index=False)

    def tearDown(self):
        if os.path.exists(self.csv_path):
            os.remove(self.csv_path)
        if os.path.exists(self.model_path):
            os.remove(self.model_path)

    def test_ingestion(self):
        agent = DataIngestionAgent(self.csv_path)
        df = agent.load_data()
        self.assertEqual(len(df), 3)

    def test_preprocessing(self):
        df = pd.read_csv(self.csv_path)
        agent = PreprocessingAgent()
        df_proc = agent.calculate_descriptors(df)
        self.assertIn('MolWt', df_proc.columns)
        self.assertIn('LogP', df_proc.columns)

    def test_modeling(self):
        df = pd.read_csv(self.csv_path)
        # Add descriptors manually for testing or run preprocessing
        agent = PreprocessingAgent()
        df_proc = agent.calculate_descriptors(df)

        model_agent = ModelBuilderAgent(target_col='Activity')
        model, X_test, y_test = model_agent.train_model(df_proc)
        self.assertIsNotNone(model)

        # Test saving
        model_agent.save_model(self.model_path)
        self.assertTrue(os.path.exists(self.model_path))

    def test_evaluation(self):
        # Mock model and data
        class MockModel:
            def predict(self, X):
                return [1.2] * len(X)

        model = MockModel()
        X_test = [[1, 2]]
        y_test = [1.2]

        agent = EvaluationAgent()
        metrics = agent.evaluate_model(model, X_test, y_test)
        self.assertEqual(metrics['RMSE'], 0.0)

if __name__ == '__main__':
    unittest.main()
