""" Gaussian function with noise. """
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Create data
x = np.linspace(0, 10, 100)
y = 3 * np.exp(-x) + np.random.normal(0, 0.5, 100)

# Fit the data
def func(x, a, b, c):
    return a * np.exp(-b * x) + c

popt, pcov = curve_fit(func, x, y)

# Plot the data
plt.plot(x, y, 'ko', label="Original Noised Data")

# Plot the fitted curve
plt.plot(x, func(x, *popt), 'r-', label="Fitted Curve")

plt.legend()
plt.show()

# Print the parameters
print("a = %s, b = %s, c = %s" % (popt[0], popt[1], popt[2]))

# Print the covariance
print(pcov)

# Print the standard deviation
perr = np.sqrt(np.diag(pcov))
print(perr)

# Print the reduced chi-squared
residuals = y - func(x, *popt)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((y-np.mean(y))**2)
r_squared = 1 - (ss_res / ss_tot)
print(r_squared)

# Print the standard error
from scipy.stats import chi2
chi2_stat = np.sum(((y - func(x, *popt)) / perr)**2)
print(chi2_stat)
print(chi2.ppf(0.95, 97))

# Print the confidence interval
from scipy.stats import t
confidence = 0.95
n = len(y)    # number of data points
m = popt.size # number of parameters
df = n - m    # number of degrees of freedom
tval = t.ppf(confidence, df)

for i, p,var in zip(range(m), popt, np.diag(pcov)):
    sigma = var**0.5
    print('p{0}: {1} [{2}  {3}]'.format(i, p,
                                        p - sigma*tval,
                                        p + sigma*tval))

