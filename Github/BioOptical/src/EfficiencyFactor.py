#this is the python script version of deriveSIOP_Mie.ipynb

import numpy as np
import sympy as sym
from scipy.special import legendre,jv,hankel2
#import matplotlib.pyplot as plt
from scipy import integrate
class Functions:
    def __init__(self):
        return
    
    def pi_n(self,n,theta):
        """
        thetaは度数で与えること
        p124中頃の式だとsinθで割っているが、
        これはルジャンドル多項式の変数がcosのときに微分することでsinの項が出てくるため。（おそらく-を忘れている誤植。）
        """
        cos = np.cos(np.radians(theta))
        legendre_func = legendre(n)
        legendre_func_derivative = legendre_func.deriv()
        return legendre_func_derivative(cos)
    
    def tau_n(self,n,theta):
        cos = np.cos(np.radians(theta))
        sin = np.sin(np.radians(theta))
        legendre_func = legendre(n)
        legendre_func_derivative = legendre_func.deriv()
        legendre_func_derivative_2 = legendre_func_derivative.deriv()

        return cos * self.pi_n(n,theta)- sin * legendre_func_derivative_2(cos)
    
    def psi_n(self,n,x):
        """
        psi_n is a key function of an and bn
        n(int): the order of the bessel function
        x(float): x= 2pi・a/λ
        a(float):The diameter of a single particle 
        """
        bessel = jv(n+1/2, x)
        return np.sqrt(np.pi*x/2)*bessel

    def xi_n(self,n,x):
        """
        xi_n is a key function of an and bn
        n(int): the order of the bessel function
        x(float): x= 2pi・a/λ
        a(float):The diameter of a single particle 
        """
        bessel = hankel2(n+1/2, x)
        return np.sqrt(np.pi*x/2)*bessel
        
    def diff_f(self,f,n,x,h):
        """
        psi_n,xi_nの数値微分用の関数。hを決めると微分係数の近似値が出る。
        bessel関数なので数値微分するしかない
        f(func): function you want to deriverate
        n(int): the order of the bessel function
        x(float): x= 2pi・a/λ
        h(float): h for limit h→0 (h should be near 0)
        """
        return (f(n,x+h)-f(n,x-h))/(2.0*h)
    
    def a_n(self,n,x,y,m,h=1e-6):
        """
        S1(θ),S2(θ)を求めるための関数。
        n(int): the order of the bessel function
        x(float): x= 2pi・a/λ
        y(float):y = mka
        m(float or complex?):refractive index (complex?) or real part?
        11/22 mはcomplex と思われる
        """
        psi_derivative_x = self.diff_f(self.psi_n,n,x,h) #psiのx=xでの微分係数
        xi_derivative_x = self.diff_f(self.xi_n,n,x,h)#xiのx=xでの微分係数
        psi_derivative_y = self.diff_f(self.psi_n,n,y,h)#psiのy=yでの微分係数
        #xi_derivative_y = diff_f(xi_n,n,y,h)

        numerator = psi_derivative_y * self.psi_n(n,x) - m * self.psi_n(n,y) *psi_derivative_x # anの分子
        denominator = psi_derivative_y * self.xi_n(n,x) - m * self.psi_n(n,y) *xi_derivative_x # anの分母
        
        return numerator/denominator

    
    def b_n(self,n,x,y,m,h=1e-6):
        """
        S1(θ),S2(θ)を求めるための関数。
        n(int): the order of the bessel function
        x(float): x= 2pi・a/λ
        y(float):y = mka
        m(float or complex?):refractive index (complex?) or real part?   
        """
        psi_derivative_x = self.diff_f(self.psi_n,n,x,h)
        xi_derivative_x = self.diff_f(self.xi_n,n,x,h)
        psi_derivative_y = self.diff_f(self.psi_n,n,y,h)
        #xi_derivative_y = diff_f(xi_n,n,y,h)

        numerator = m * psi_derivative_y * self.psi_n(n,x) -  self.psi_n(n,y) *psi_derivative_x
        denominator = m * psi_derivative_y * self.xi_n(n,x) -  self.psi_n(n,y) *xi_derivative_x
        
        return numerator/denominator
    
    def S_1(self,theta,x,y,m,nmax=30):
        """
        S1(θ)を求めるための関数。
        theta(int): 散乱角。入射角の方向が0度である。
        x(float): x= 2pi・a/λ
        y(float):y = mka
        m(float or complex?): refractive index (complex?) or real part?  
        nmax(int,optional): since this function is infinite series, the upper limit for iteration should be determined.
        """
        sum = 0
        for n in range(1,nmax):
            numerator_1 = 2*n+1
            numerator_2 = self.a_n(n,x,y,m) * self.pi_n(n,theta) + self.b_n(n,x,y,m) * self.tau_n(n,theta)
            denominator = n * (n+1)
            #print( numerator_1,numerator_2,denominator)
            sum += (numerator_1*numerator_2)/denominator
        return sum

    def S_2(self,theta,x,y,m,nmax=30):
        """
        S2(θ)を求めるための関数。
        theta(int): 散乱角。入射角の方向が0度である。
        x(float): x= 2pi・a/λ
        y(float):y = mka
        m(float or complex?): refractive index (complex?) or real part?  
        nmax(int,optional): since this function is infinite series, the upper limit for iteration should be determined.
        """
        sum = 0
        for n in range(1,nmax):
            numerator_1 = 2*n+1
            numerator_2 = self.b_n(n,x,y,m) * self.pi_n(n,theta) + self.a_n(n,x,y,m) * self.tau_n(n,theta)
            denominator = n * (n+1)
            sum += (numerator_1*numerator_2)/denominator
        return sum    

    def cal_abs(self,theta,x,y,m):
        """
        i_1,i_2を求める。S1,S2は虚数なので絶対値の二乗を求める。これにより実数となる。
        theta(int): 散乱角。入射角の方向が0度である。
        x(float): x= 2pi・a/λ
        y(float):y = mka
        m(float or complex?): refractive index (complex?) or real part?  
        
        return:
            i_1,i_2 (tuple)
        """
        i_1 = abs(self.S_1(theta,x,y,m)) ** 2
        i_2 = abs(self.S_1(theta,x,y,m)) ** 2
        return i_1,i_2


    def Qsca(self,x,y,m):
        def integrate_scattering(theta):
            return (self.cal_abs(theta,x,y,m)[0] + self.cal_abs(theta,x,y,m)[0]) * np.sin(np.radians(theta)) #thetaのいち変数とするための関数。

        integrated,err = integrate.quad(integrate_scattering,0,np.pi) #sympyの数値積分メソッドを使用。同時に推定誤差も返される。
        return integrated/(x**2)
    
    def Qbsca(self,x,y,m):
        def integrate_scattering(theta):
            return (self.cal_abs(theta,x,y,m)[0] + self.cal_abs(theta,x,y,m)[0]) * np.sin(np.radians(theta))

        integrated,err = integrate.quad(integrate_scattering,np.pi/2,np.pi)
        return integrated/(x**2)
    
    def Intensity(self,I0,theta,x,y,a,m,r):
        i_1,i_2 = self.cal_abs(theta,x,y,m)
        k = x/a
        numerator = I0 * (i_1+i_2)
        denominator = 2 * (k**2) * (r**2)

        return numerator/denominator
        
