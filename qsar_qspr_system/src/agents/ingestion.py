import pandas as pd

class DataIngestionAgent:
    def __init__(self, filepath):
        self.filepath = filepath

    def load_data(self):
        """Loads data from a CSV file."""
        try:
            df = pd.read_csv(self.filepath)
            print(f"Data loaded successfully from {self.filepath}")
            return df
        except Exception as e:
            print(f"Error loading data: {e}")
            return None
