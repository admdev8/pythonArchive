import os 
import numpy as np
import scipy.optimize as optimize
from scipy import linspace,polyfit, polyval, poly1d, divide, mean

##################
# HELPER METHODS #
##################
# Method to calculate the taylor derivitive
def diffTaylor(v,i): 
    coeffs = polyfit(v,i,3)
    vFit = linspace(min(v),max(v),100) 
    iFit = polyval(coeffs, vFit) 
    polynomial = poly1d(coeffs)
   
    #didv = lambda ix: 6*coeffs[0]*((ix)**5) + 5*coeffs[1]*((ix)**4) + 4*coeffs[2]*((ix)**3) + 3*coeffs[3]*((ix)**2) + 2*coeffs[4]*(ix) + coeffs[5]
    didv = lambda ix:  3*coeffs[0]*((ix)**2) + 2*coeffs[1]*(ix) + coeffs[2]
    g = [float(didv(val)) for val in vFit]  
    rFit = divide(1,g)
    return vFit, iFit, rFit


# Method to check specific paths against what one is intersted in
def checkString(modifier, path): 
    for ck in modifier: 
        if ck in path:
            pass
        else:
            return None
 
    return True

# Calculate Correlation coefficient
def correlation(dat1,dat2): 
    top=mean(np.multiply((dat1-mean(dat1)),(dat2-mean(dat2))))
    bot=np.std(dat1)*np.std(dat2)
    return top/bot


if __name__ == "__main__": 

    root = os.getcwd()
    modifier = ['postBCB', 'chipJ.png']

    # Get the Hall path
    for path, dirs, files in os.walk(root):
        for filename in files:
            if checkString(modifier,os.path.join(path, filename)):
                hallPath = os.path.join(path, filename)
                print hallPath



