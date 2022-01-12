import numpy as np
from scipy import linalg as LA
import math
def lebedev():
    file_name = "F:/calculate/lebedev/lebedev1.txt"
    file_object = open(file_name)
    data1 = file_object.readlines()
    data = []
    for i in data1:
        j = i.split()
        l = [float(j[0])/180*math.pi,float(j[1])/180*math.pi,float(j[2])]
        data.append(l)
    return data

def Judge1(x,y,z,i,position):
    position = [[0.6986449763, -5.515579545, 0.07220265577], [-0.4350206231, -5.515579545, 0.07220265577]]
    tot = 1
    for j in range(len(position)):
        if j != i:
            mu = (math.sqrt(x**2+y**2+z**2)-math.sqrt((x+position[i][0]-position[j][0])**2+(y+position[i][1]-position[j][1])**2+(z+position[i][2]-position[j][2])**2))/\
                 math.sqrt((position[j][0]-position[i][0])**2+(position[j][1]-position[i][1])**2+(position[j][2]-position[i][2])**2)
            tot = tot*(0.5-0.5*RE(mu,3))
    total = []
    total.append(tot)
    for j in range(len(position)):
        tot1 = 1
        if j!= i :
            for k in range(len(position)):
                if k != j:
                    mu =(math.sqrt((x+position[i][0]-position[j][0])**2+(y+position[i][1]-position[j][1])**2+(z+position[i][2]-position[j][2])**2)-math.sqrt((x+position[i][0]-position[k][0])**2+(y+position[i][1]-position[k][1])**2+(z+position[i][2]-position[k][2])**2))/ \
                         math.sqrt((position[j][0]-position[k][0])**2+(position[j][1]-position[k][1])**2+(position[j][2]-position[k][2])**2)
                    tot1 = tot1 * (0.5 - 0.5 * RE(mu, 3))
            total.append(tot1)
    print(total)
    return tot/sum(total)

def RE(mu,i):
    if i == 0 :
        return mu
    return RE(1.5*mu-0.5*mu**3,i-1)

def density():
    p = 1
    n = 30
    D = 0
    L = lebedev()
    charge = [1,1]
    for i in range(len(charge)):
        for k in range(1,n+1):
            x = math.cos((2*k-1)/2/n*math.pi)
            for s in L:
                fai = s[1]
                theta = s[0]
                D = D + s[2]*4*math.pi*math.pi/n*2*p**3*(1+x)**2.5/(1-x)**3.5*f(p*(1+x)/(1-x)*math.sin(fai)*math.cos(theta),p*(1+x)/(1-x)*math.sin(fai)*math.sin(theta),p*(1+x)/(1-x)*math.cos(fai),i)*\
                Judge1(p*(1+x)/(1-x)*math.sin(fai)*math.cos(theta),p*(1+x)/(1-x)*math.sin(fai)*math.sin(theta),p*(1+x)/(1-x)*math.cos(fai),i)
    return D
def test1(x,y,z,pos):
    g = 22.5377*math.exp(-99.9997*(x**2+y**2+z**2))
    return g**2     # this code is designed just to check the correctness of Density Calculation.
def f(x,y,z,pos):
    charge = [1,1]
    position = [[0.6986449763, -5.515579545, 0.07220265577], [-0.4350206231, -5.515579545, 0.07220265577]]
    orbital_coe = [[0.2769343551, 0.2678388516, 0.08347367113], [0.2769343551, 0.2678388516, 0.08347367113]]
    orbital_exp = [[3.425250914, 0.6239137298, 0.1688554040], [3.425250914, 0.6239137298, 0.1688554040]]
    orbital_pos = [[0.6986449763, -5.515579545, 0.07220265577], [-0.4350206231, -5.515579545, 0.07220265577]]
    orbital_typ = [[0,0,0],[0,0,0]]
    X = [[-0.53429946, -1.41836484],
         [-0.53429946, 1.41836484]]
    rho_sqr = 0
    rho = 0
    for i in range(int(sum(charge)/2)):
        for j in range(len(X)):
            for k in range(len(orbital_coe[j])):
                rho_sqr = rho_sqr+X[j][i]*orbital_coe[j][k]*(x+position[pos][0]-orbital_pos[j][0])**orbital_typ[j][0]*(y+position[pos][1]-orbital_pos[j][1])**orbital_typ[j][1]*(z+position[pos][2]-orbital_pos[j][2])**orbital_typ[j][2]*math.exp(-orbital_exp[j][k]*((x+position[pos][0]-orbital_pos[j][0])**2+(y+position[pos][1]-orbital_pos[j][1])**2+(z+position[pos][2]-orbital_pos[j][2])**2))
        rho = rho+rho_sqr**2
        rho_sqr = 0
    return rho*2
print(density())