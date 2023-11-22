from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt

from EfficiencyFactor import Functions
func = Functions()
#test
print(func.pi_n(3,30))

def PSD(D,j,K=1):
    """
    particle size distribution(PSD) N(D) = K*D^-j
    D: diameter of the particle
    j: parameter1
    K(optional): parameter2
    """
    return K * D ** (-j)


#derive Qbb
def Qbb(D,Lambda,m):
    x = 2 * np.pi * D / Lambda
    y = m * x
    Q_backscattering = func.Qbsca(x,y,m)
    return Q_backscattering


def Cbb(Lambda,m,j,dmin,dmax):
    def numerator(D):
        Q_bb = Qbb(D,Lambda,m)
        ND = PSD(D,j)
        return Q_bb * np.pi * D**2 * ND /4
    
    def denominator(D):
        ND = PSD(D,j)
        return ND
    
    integrated_numerator,err_numerator = integrate.quad(numerator,dmin,dmax)
    intergrated_denominator,err_denominator = integrate.quad(denominator,dmin,dmax)
    return integrated_numerator/intergrated_denominator


def G(j,dmin,dmax):
    def numerator(D):
        ND = PSD(D,j)
        return ND * (D **2)
    
    def denominator(D):
        ND = PSD(D,j)
        return ND
    
    integrated_numerator,err_numerator = integrate.quad(numerator,dmin,dmax)
    intergrated_denominator,err_denominator = integrate.quad(denominator,dmin,dmax)
    return integrated_numerator/intergrated_denominator


def bbp_SIOP(Lambda,j,m,pho,dmin,dmax):
    Qbb_ave = Cbb(Lambda,m,j,dmin,dmax)/G(j,dmin,dmax)

    def numerator(D):
        ND = PSD(D,j)
        return ND * (D **2)
    
    def denominator(D):
        ND = PSD(D,j)
        return ND * (D **3)
    
    integrated_numerator,err_numerator = integrate.quad(numerator,dmin,dmax)
    intergrated_denominator,err_denominator = integrate.quad(denominator,dmin,dmax) 

    return 3 * Qbb_ave * integrated_numerator / (2 * pho * intergrated_denominator)  



dmin = 0.27 * (1e-6)
dmax = 240 * (1e-6)
j = -4 # -4 is a case for open ocean
pho = 2.5
lambda_min = 4 #* (1e-9)
lambda_max = 8 #* (1e-9)
m = [complex(1.17, 0.015*((l*1e-7) ** -0.004))  for l in range(lambda_min,lambda_max)]
Lambda = [l*(1e-7) for l  in range(lambda_min,lambda_max)]

SIOP = [bbp_SIOP(l,j,m,pho,dmin,dmax) for l,m in zip(Lambda,m)]


fig,ax = plt.subplots()
ax.plot(Lambda,SIOP)
fig.savefig("../png/SIOP_backscattering_particle.png",dpi=600)