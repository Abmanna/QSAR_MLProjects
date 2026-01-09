from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.svm import SVR
from sklearn.preprocessing import StandardScaler
import joblib

class ModelBuilderAgent:
    def __init__(self, target_col, model_type='rf'):
        self.target_col = target_col
        self.model_type = model_type
        self.scaler = StandardScaler()
        self.model = self._get_model(model_type)

    def _get_model(self, model_type):
        if model_type == 'rf':
            return RandomForestRegressor(random_state=42)
        elif model_type == 'gb':
            return GradientBoostingRegressor(random_state=42)
        elif model_type == 'svr':
            return SVR()
        else:
            raise ValueError(f"Unknown model type: {model_type}")

    def tune_hyperparameters(self, X_train, y_train):
        """Performs hyperparameter tuning using GridSearchCV."""
        param_grid = {}
        if self.model_type == 'rf':
            param_grid = {
                'n_estimators': [50, 100, 200],
                'max_depth': [None, 10, 20],
                'min_samples_split': [2, 5]
            }
        elif self.model_type == 'gb':
            param_grid = {
                'n_estimators': [50, 100, 200],
                'learning_rate': [0.01, 0.1, 0.2],
                'max_depth': [3, 5, 7]
            }
        elif self.model_type == 'svr':
            param_grid = {
                'C': [0.1, 1, 10],
                'gamma': ['scale', 'auto'],
                'kernel': ['rbf', 'linear']
            }

        print(f"Tuning hyperparameters for {self.model_type}...")
        grid_search = GridSearchCV(self.model, param_grid, cv=3, scoring='neg_mean_squared_error', n_jobs=-1)
        grid_search.fit(X_train, y_train)

        print(f"Best params: {grid_search.best_params_}")
        return grid_search.best_estimator_

    def train_model(self, df, tune=False):
        """Trains the model."""
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

        if tune:
            self.model = self.tune_hyperparameters(X_train_scaled, y_train)
        else:
            self.model.fit(X_train_scaled, y_train)

        print("Model trained.")

        return self.model, X_test_scaled, y_test

    def save_model(self, filepath):
        joblib.dump(self.model, filepath)
        print(f"Model saved to {filepath}")
