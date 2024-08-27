import numpy as np
import scipy.interpolate as interp
import scipy.optimize as opt
import argparse

def read_data(filename):
    """Reads the data from a file and returns q and I(q) as numpy arrays."""
    data = np.loadtxt(filename)
    q = data[:, 0]
    Iq = data[:, 1]
    return q, Iq

def interpolate_function(q, Iq, q_new):
    """Interpolates the function using a spline and evaluates it on a new grid."""
    spline = interp.InterpolatedUnivariateSpline(q, Iq, k=3)  # Cubic spline
    Iq_new = spline(q_new)
    return Iq_new

def fit_function(q, Iq, q_min=1.7, q_max=2.5):
    """Fits the function in the specified q range using a linear model."""
    mask = (q > q_min) & (q <= q_max)
    popt, _ = opt.curve_fit(lambda x, a, b: a * x + b, q[mask], Iq[mask])
    return popt

def main():
    parser = argparse.ArgumentParser(description="Process two files containing function data.")
    parser.add_argument("file1", help="Path to the first file (F_exp)")
    parser.add_argument("file2", help="Path to the second file (F_sol)")
    parser.add_argument("output", help="Path to the output file")
    parser.add_argument("--q_min", type=float, default=0.0, help="Minimum q value for the grid (default 0.0)")
    parser.add_argument("--q_max", type=float, default=2.5, help="Maximum q value for the grid (default 2.5)")
    parser.add_argument("--q_step", type=float, default=0.01, help="Step size for the grid (default 0.01)")

    args = parser.parse_args()

    # Read data from files
    q1, Iq1 = read_data(args.file1)
    q2, Iq2 = read_data(args.file2)

    # Create a regular grid
    q_new = np.arange(args.q_min, args.q_max + args.q_step, args.q_step)

    # Interpolate the functions onto the regular grid
    Iq1_new = interpolate_function(q1, Iq1, q_new)
    Iq2_new = interpolate_function(q2, Iq2, q_new)

    # Fit both functions in the region q > 1.7 and q ≤ 2.5
    fit1 = fit_function(q_new, Iq1_new)
    fit2 = fit_function(q_new, Iq2_new)

    # Apply the fits to the functions
    Iq1_fit = fit1[0] * q_new + fit1[1]
    Iq2_fit = fit2[0] * q_new + fit2[1]

    # Compute the difference
    difference = Iq1_fit - Iq2_fit

    # Write the result to the output file
    np.savetxt(args.output, np.column_stack((q_new, difference)), fmt="%.6f", header="q I_diff")

if __name__ == "__main__":
    main()
