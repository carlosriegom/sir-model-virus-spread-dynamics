import numpy as np
from scipy.optimize import fsolve

def find_a_general(k):
    """Calculate the parameter 'a' for a given recovery rate 'k'.
    
    Args:
    k (float): Recovery rate, defined as 1/recovery time in days.
    
    Returns:
    float: Value of 'a' that sets the 50% recovery probability at 1/k days.
    """
    # Convert '1/k' recovery rate into hours
    recovery_time_hours = (1 / k) * 24  # Convert days to hours
    
    # Define the function to solve for 'a' such that the recovery probability is 0.5 at 'recovery_time_hours'
    def equation(a):
        total_time = 720  # 30 days in hours for normalization
        return (1 - np.exp(-a * recovery_time_hours)) / (1 - np.exp(-a * total_time)) - 0.5
    
    # Initial guess for 'a'
    a_initial_guess = 0.01
    # Solve for 'a'
    a_solution = fsolve(equation, a_initial_guess)[0]
    return a_solution

# Example usage of the function
if __name__ == "__main__":
    # Define a recovery rate 'k', e.g., k = 1/5 for a recovery time of 5 days
    k_value = 1/7  # Change this value to test different recovery rates
    a_value = find_a_general(k_value)
    print(f"The value of 'a' for a recovery time of {1/k_value:.1f} days is approximately {a_value:.5f}")
