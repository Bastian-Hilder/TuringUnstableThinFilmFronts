from sympy import *
from sympy.parsing.mathematica import parse_mathematica
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from Difference_M import Mdiff

b, g = var('b g')

######### read mathematica expression from external files and convert to python expressions ##########
print("loading mathematica expressions\n")

# load quadratic coefficient hexagon
with open('quadratic-coeff-hex.txt','r') as file:
    NCoeff = file.read()

NCoeff = NCoeff.replace('{','')
NCoeff = NCoeff.replace('}','')
N_python = parse_mathematica(NCoeff)

N_fast_eval = lambdify([b,g],N_python)  

# load self coefficient on N=0
with open('K0-coeffs-on-Nzero.txt','r') as file:
    K0N0 = file.read()

K0N0 = K0N0.replace('{','')
K0N0 = K0N0.replace('}','')
K0N0_python = parse_mathematica(K0N0)

K0N0_fast_eval = lambdify([g],K0N0_python)  

# load cross coefficient on N=0
with open('K2-coeffs-on-Nzero.txt','r') as file:
    K2N0 = file.read()

K2N0 = K2N0.replace('{','')
K2N0 = K2N0.replace('}','')
K2N0_python = parse_mathematica(K2N0)

K2N0_fast_eval = lambdify([g],K2N0_python) 

############## generate grid ##############
nbeta = 250
ng = 250
betaGrid = np.linspace(0.01,5,nbeta)
gGrid = np.linspace(0.01,20,ng)
X , Y = np.meshgrid(betaGrid,gGrid,indexing='xy')
n = np.zeros((len(gGrid),len(betaGrid)))

# define function for Mm*-Mo* to determine if points are in the relevant parameter area
MmMinusMo = Mdiff()

############## evaluate function on grid ##############
# write z-values by evaluating kappajj
for ii in range(len(betaGrid)):
    print(str(ii)+"/"+str(len(betaGrid)))
    for jj in range(len(gGrid)):
        if MmMinusMo(X[jj,ii],Y[jj,ii]) >= 0:
            n[jj,ii] = N_fast_eval(X[jj,ii],Y[jj,ii])
        else:
            n[jj,ii] = np.nan


# print(str(2*np.heaviside(-K0_fast_eval(9,0.5),0.5) + np.heaviside(-K0_fast_eval(9,0.5)-K1_fast_eval(9,0.5),0.5)))
# print(str(2*np.heaviside(-K0_fast_eval(9,3),0.5) + np.heaviside(-K0_fast_eval(9,3)-K1_fast_eval(9,3),0.5)))
# print(str(2*np.heaviside(-K0_fast_eval(9,5),0.5) + np.heaviside(-K0_fast_eval(9,5)-K1_fast_eval(9,5),0.5)))

# plot and save results
# coefficients: N
fig, ax = plt.subplots()
densityN = ax.pcolormesh(X,Y,n,cmap='jet')
ax.contour(X,Y,n,[0],colors='w',linestyles='dashed')
fig.colorbar(densityN)

plt.title(r'$N$')
plt.xlabel('β')
plt.ylabel('g')
plt.savefig("DensityMapN.png",dpi=500)


################ plot coefficients on the line {N=0} ##################
x = np.linspace(0.01,17.99,500)
k0n0 = np.zeros(len(x))
k2n0 = np.zeros(len(x))


for ii in range(len(x)):
    k0n0[ii] = K0N0_fast_eval(x[ii])
    k2n0[ii] = K2N0_fast_eval(x[ii])

fig, ax = plt.subplots()
ax.plot(x,k0n0,label=r'$K_0|_{\{N=0\}}$')
ax.plot(x,k2n0,label=r'$K_2|_{\{N=0\}}$')
plt.xlim(0,18)
plt.xlabel('g')
plt.title(r'Coefficients on $\{N=0\}$')
plt.legend()

plt.savefig("Coeffs-on-Nzero.png",dpi=500)