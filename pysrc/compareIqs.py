import argparse
import numpy as np
from scipy import interpolate
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def load_data(filename):
    return np.loadtxt(filename, unpack=True)

def interpolate_data(x, y, new_x):
    spline = interpolate.splrep(x, y, s=0)
    return interpolate.splev(new_x, spline, der=0)

def fit_function(x, a, b, c):
    return a * np.exp(-b * x) + c

def process_functions(file1, file2, output_file, plot=False):
    # Load data
    x1, y1 = load_data(file1)
    x2, y2 = load_data(file2)

    # Create a common x grid
    x_min = max(x1.min(), x2.min())
    x_max = min(x1.max(), x2.max())
    x_common = np.linspace(x_min, x_max, 1000)

    # Interpolate both functions to the common grid
    y1_interp = interpolate_data(x1, y1, x_common)
    y2_interp = interpolate_data(x2, y2, x_common)

    # Fit both functions in the region where q > 1.7 and q <= 2.5
    fit_mask = (x_common > 1.7) & (x_common <= 2.5)
    x_fit = x_common[fit_mask]
    y1_fit = y1_interp[fit_mask]
    y2_fit = y2_interp[fit_mask]

    popt1, _ = curve_fit(fit_function, x_fit, y1_fit)
    popt2, _ = curve_fit(fit_function, x_fit, y2_fit)

    # Compute the difference between the fitted functions
    y1_fitted = fit_function(x_common, *popt1)
    y2_fitted = fit_function(x_common, *popt2)
    y_diff = y1_fitted - y2_fitted

    # Save the result to a file
    np.savetxt(output_file, np.column_stack((x_common, y_diff)), delimiter='\t',
               header='q\tDifference', comments='')

    if plot:
        plt.figure(figsize=(10, 6))
        plt.plot(x_common, y1_interp, label='F_exp interpolated')
        plt.plot(x_common, y2_interp, label='F_sol interpolated')
        plt.plot(x_common, y1_fitted, label='F_exp fitted')
        plt.plot(x_common, y2_fitted, label='F_sol fitted')
        plt.plot(x_common, y_diff, label='Difference')
        plt.xlabel('q')
        plt.ylabel('I(q)')
        plt.legend()
        plt.title('Function Comparison')
        plt.show()

def main():
    parser = argparse.ArgumentParser(description='Process and compare two functions.')
    parser.add_argument('file1', help='Path to the first input file (F_exp)')
    parser.add_argument('file2', help='Path to the second input file (F_sol)')
    parser.add_argument('output', help='Path to the output file')
    parser.add_argument('--plot', action='store_true', help='Plot the results')

    args = parser.parse_args()

    process_functions(args.file1, args.file2, args.output, args.plot)

if __name__ == '__main__':
    main()