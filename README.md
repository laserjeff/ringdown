# ringdown
Environment for ringdown analysis

Fits data to an Exponential modified Gaussian (EMG) and filters based on R^2.

EMG:
def exp_mod_gauss(x, b, m, s, l):
    y = b*(0.5*l*np.exp(0.5*l*(2*m+l*s*s-2*x))*erfc((m+l*s*s-x)/(np.sqrt(2)*s)))
    return y
    #l=Lambda, s=Sigma, m=Mu, #b=scaling
    
Calculates average ringdown time for remaining samples.
