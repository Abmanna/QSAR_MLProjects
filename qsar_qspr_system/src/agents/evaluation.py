from sklearn.metrics import mean_squared_error, r2_score
import numpy as np

class EvaluationAgent:
    def __init__(self):
        pass

    def evaluate_model(self, model, X_test, y_test):
        """Evaluates the model and prints metrics."""
        predictions = model.predict(X_test)

        mse = mean_squared_error(y_test, predictions)
        rmse = np.sqrt(mse)
        r2 = r2_score(y_test, predictions)

        print(f"RMSE: {rmse}")
        print(f"R2 Score: {r2}")

        return {'RMSE': rmse, 'R2': r2}
