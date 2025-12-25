from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler
import joblib

class ModelBuilderAgent:
    def __init__(self, target_col):
        self.target_col = target_col
        self.model = RandomForestRegressor(n_estimators=100, random_state=42)
        self.scaler = StandardScaler()

    def train_model(self, df):
        """Trains a Random Forest model."""
        # Assume all non-target columns are features (except SMILES if present)
        features = df.drop(columns=[self.target_col])
        if 'SMILES' in features.columns:
            features = features.drop(columns=['SMILES'])

        X = features
        y = df[self.target_col]

        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

        # Scale features
        X_train_scaled = self.scaler.fit_transform(X_train)
        X_test_scaled = self.scaler.transform(X_test)

        self.model.fit(X_train_scaled, y_train)
        print("Model trained.")

        return self.model, X_test_scaled, y_test

    def save_model(self, filepath):
        joblib.dump(self.model, filepath)
        print(f"Model saved to {filepath}")
