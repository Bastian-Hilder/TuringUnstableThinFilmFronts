from sympy import *
from sympy.parsing.mathematica import parse_mathematica
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from Difference_M import Mdiff

b, g = var('b g')

# specify colors
tumblue = '#64A0C8'
tumorange = '#E37222'
tumgreen = '#A2AD00'
tumivory = '#DAD7CB'

######### read mathematica expression from external files and convert to python expressions ##########
print("loading mathematica expressions\n")
# load self coefficient #
with open('Self-cubic-coeffs-hex.txt','r') as file:
    K0 = file.read()

K0 = K0.replace('{','')
K0 = K0.replace('}','')
K0_python = parse_mathematica(K0)

K0_fast_eval = lambdify([b,g],K0_python)

# load cross coefficient square
with open('cross-cubic-coeffs-square.txt','r') as file:
    K1 = file.read()

K1 = K1.replace('{','')
K1 = K1.replace('}','')
K1_python = parse_mathematica(K1)

K1_fast_eval = lambdify([b,g],K1_python)    

############## generate grid ##############
nbeta = 2000
ng = 1000
betaGrid = np.linspace(0.01,50,nbeta)
gGrid = np.linspace(0.1,20,ng)
X , Y = np.meshgrid(betaGrid,gGrid,indexing='xy')
existence = np.zeros((len(gGrid),len(betaGrid)))
Mdifference = np.zeros((len(gGrid),len(betaGrid)))
k0 = np.zeros((len(gGrid),len(betaGrid)))
k1 = np.zeros((len(gGrid),len(betaGrid)))

# define function for Mm*-Mo* to determine if points are in the relevant parameter area
MmMinusMo = Mdiff()

############## evaluate function on grid ##############
# write z-values by evaluating kappajj
for ii in range(len(betaGrid)):
    print(str(ii)+"/"+str(len(betaGrid)))
    for jj in range(len(gGrid)):
        k0[jj,ii] = K0_fast_eval(X[jj,ii],Y[jj,ii])
        k1[jj,ii] = K1_fast_eval(X[jj,ii],Y[jj,ii])
        Mdifference[jj,ii] = MmMinusMo(X[jj,ii],Y[jj,ii])
        if Mdifference[jj,ii] >= 0:
            existence[jj,ii] = np.heaviside(-k0[jj,ii],0) + 2*np.heaviside(-k0[jj,ii]-k1[jj,ii],0) + 1
        else:
            k0[jj,ii] = np.nan
            k1[jj,ii] = np.nan


# print(str(2*np.heaviside(-K0_fast_eval(9,0.5),0.5) + np.heaviside(-K0_fast_eval(9,0.5)-K1_fast_eval(9,0.5),0.5)))
# print(str(2*np.heaviside(-K0_fast_eval(9,3),0.5) + np.heaviside(-K0_fast_eval(9,3)-K1_fast_eval(9,3),0.5)))
# print(str(2*np.heaviside(-K0_fast_eval(9,5),0.5) + np.heaviside(-K0_fast_eval(9,5)-K1_fast_eval(9,5),0.5)))

# plot and save results
fig, ax = plt.subplots()
cs_existence = ax.pcolormesh(X,Y,existence,cmap=mpl.colors.ListedColormap(['white',tumblue,tumgreen,tumivory,tumorange]))
proxy = [plt.Rectangle((0, 0), 1, 1, fc=fc) for fc in ['white',tumblue,tumgreen,tumivory,tumorange]]
plt.legend(proxy, [r'$\Omega_o$',"none", "rolls", "squares", "squares & rolls"])
# plt.show()

plt.title('Existence for square lattice')
plt.xlabel('β')
plt.ylabel('g')
plt.savefig("ExistenceSquareLatticeLargeDomain.png",dpi=500)

# coefficients: K0
fig, ax = plt.subplots()
densityK0 = ax.pcolormesh(X,Y,k0,cmap='jet')
ax.contour(X,Y,k0,[0],colors='w',linestyles='dashed')
fig.colorbar(densityK0)

plt.title(r'$K_0$')
plt.xlabel('β')
plt.ylabel('g')
plt.savefig("DensityMapK0LargeDomain.png",dpi=500)

# coefficients: K1
fig, ax = plt.subplots()
densityK1 = ax.pcolormesh(X,Y,k1,cmap='jet')
ax.contour(X,Y,k1,[0],colors='w',linestyles='dashed')
fig.colorbar(densityK1)

plt.title(r'$K_1$')
plt.xlabel('β')
plt.ylabel('g')
plt.savefig("DensityMapK1LargeDomain.png",dpi=500)