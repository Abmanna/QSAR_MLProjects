from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler
import joblib
import pandas as pd

class ModelBuilderAgent:
    def __init__(self, target_col):
        self.target_col = target_col
        self.model = RandomForestRegressor(n_estimators=100, random_state=42)
        self.scaler = StandardScaler()
        self.feature_columns = None # To keep track of column order

    def train_model(self, df):
        """Trains a Random Forest model."""
        # Assume all non-target columns are features (except SMILES if present)
        features = df.drop(columns=[self.target_col])
        if 'SMILES' in features.columns:
            features = features.drop(columns=['SMILES'])

        self.feature_columns = features.columns.tolist()

        X = features
        y = df[self.target_col]

        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

        # Scale features
        X_train_scaled = self.scaler.fit_transform(X_train)
        X_test_scaled = self.scaler.transform(X_test)

        self.model.fit(X_train_scaled, y_train)
        print("Model trained.")

        return self.model, X_test_scaled, y_test

    def predict(self, df):
        """
        Predicts target values for new data.
        Args:
            df (pd.DataFrame): Dataframe containing features.
        Returns:
            np.array: Predicted values.
        """
        # Ensure we have the same columns
        if 'SMILES' in df.columns:
            features = df.drop(columns=['SMILES'])
        else:
            features = df.copy()

        # If target col is present (e.g. during testing), drop it
        if self.target_col in features.columns:
            features = features.drop(columns=[self.target_col])

        # Align columns
        if self.feature_columns:
            # Add missing cols with 0
            for col in self.feature_columns:
                if col not in features.columns:
                    features[col] = 0
            # Reorder and select
            features = features[self.feature_columns]

        X_scaled = self.scaler.transform(features)
        return self.model.predict(X_scaled)

    def save_model(self, filepath):
        joblib.dump(self.model, filepath)
        joblib.dump(self.scaler, filepath + ".scaler")
        joblib.dump(self.feature_columns, filepath + ".cols")
        print(f"Model saved to {filepath}")
