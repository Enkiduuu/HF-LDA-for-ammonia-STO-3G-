import numpy as np
from numpy import linalg as LA
import math

def overlap(x1,x2,p1,p2,alpha,beta,Judge,y1=0,y2=0,p3=0,p4=0,z1=0,z2=0,p5=0,p6=0):
    if p1<0 or p2 <0 :
        return 0
    if p1>0 and p2 == 0:
        return -beta/(alpha+beta)*(x1-x2)*overlap(x1,x2,p1-1,p2,alpha,beta,Judge,y1=y1,y2=y2,p3=p3,p4=p4,z1=z1,z2=z2,p5=p5,p6=p6)\
                   +(p1-1)/(2*(alpha+beta))*overlap(x1,x2,p1-2,p2,alpha,beta,Judge,y1=y1,y2=y2,p3=p3,p4=p4,z1=z1,z2=z2,p5=p5,p6=p6)\
                   +p2/(2*(alpha+beta))*overlap(x1,x2,p1-1,p2-1,alpha,beta,Judge,y1=y1,y2=y2,p3=p3,p4=p4,z1=z1,z2=z2,p5=p5,p6=p6)
    if p1 == 0 and p2 >0:
        return alpha/(alpha+beta)*(x1-x2)*overlap(x1,x2,p1,p2-1,alpha,beta,Judge,y1=y1,y2=y2,p3=p3,p4=p4,z1=z1,z2=z2,p5=p5,p6=p6)\
                   +p1/(2*(alpha+beta))*overlap(x1,x2,p1-1,p2-1,alpha,beta,Judge,y1=y1,y2=y2,p3=p3,p4=p4,z1=z1,z2=z2,p5=p5,p6=p6)\
                   +(p2-1)/(2*(alpha+beta))*overlap(x1,x2,p1,p2-2,alpha,beta,Judge,y1=y1,y2=y2,p3=p3,p4=p4,z1=z1,z2=z2,p5=p5,p6=p6)
    if p1==0 and p2==0 and Judge==True:
        return Gauss(x1,x2,p1,p2,alpha,beta)
    if p1==0 and p2==0 and p3==0 and p4==0 and p5==0 and p6==0 and Judge==False:
        return F(x1,x2,p1,p2,alpha,beta,y1,y2,p3,p4,z1,z2,p5,p6)
    if Judge== False and p1==0 and p2==0 :
        return overlap(y1,y2,p3,p4,alpha,beta,Judge,y1=z1,y2=z2,p3=p5,p4=p6,z1=x1,z2=x2,p5=p1,p6=p2)
    return -beta/(alpha+beta)*(x1-x2)*overlap(x1,x2,p1-1,p2,alpha,beta,Judge,y1=y1,y2=y2,p3=p3,p4=p4,z1=z1,z2=z2,p5=p5,p6=p6)\
            +(p1-1)/(2*(alpha+beta))*overlap(x1,x2,p1-2,p2,alpha,beta,Judge,y1=y1,y2=y2,p3=p3,p4=p4,z1=z1,z2=z2,p5=p5,p6=p6)\
            +p2/(2*(alpha+beta))*overlap(x1,x2,p1-1,p2-1,alpha,beta,Judge,y1=y1,y2=y2,p3=p3,p4=p4,z1=z1,z2=z2,p5=p5,p6=p6)
def Gauss(x1,x2,p1,p2,alpha,beta):
    return math.sqrt(math.pi/(alpha+beta))*math.exp(-(x1-x2)**2*alpha*beta/(alpha+beta))
def F(x1,x2,p1,p2,alpha,beta,y1,y2,p3,p4,z1,z2,p5,p6):
    F0=0
    l = alpha*beta/(alpha+beta)*((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
    if l <20 :
        for i in range(0,int(2*l+8)):
            F0=F0+2*(math.pi)**2.5/alpha/beta/math.sqrt(alpha+beta)*(-l)**i/math.factorial(i)/(2*i+1)
    else:
        F0 = math.sqrt(1/(alpha+beta)/l)*math.pi**3/alpha/beta
    return F0
def reexpress(x1,x2,p1,p2,alpha,beta):
    y = (alpha*x1+beta*x2)/(alpha+beta)
    k = []
    for i in range(0,p1+p2+1):
        q = 0
        for j in range(max(0,i-p2),min(i,p1)+1):
            q = q+math.factorial(p1)/math.factorial(j)/math.factorial(p1-j)*(y-x1)**(p1-j)*math.factorial(p2)\
                /math.factorial(i-j)/math.factorial(p2-i+j)*(y-x2)**(p2-i+j)  #i power in total, j power of the first formula
        k.append(q)
    return k
def test1(): # It's just a function designed to test the correctness of program in the process of coding.
    a = input()
    b = input()
    c = input()
    d = input()
    a = a.split()
    b = b.split()
    c = c.split()
    d = d.split()
    a = [float(i) for i in a]
    b = [float(i) for i in b]
    c = [float(i) for i in c]
    d = [float(i) for i in d]
    sum1 = 0
    for j in range(len(a)):
        for l in range(len(a)):
            for m in range(len(a)):
                for n in range(len(a)):
                    sum1 = sum1+b[j]*b[l]*b[m]*b[n]*test2(0.30926,0.30926,0,0,a[j],a[l],0.30926,0.30926,0,0,a[m],a[n],0.30926,0.30926,0.30926,0.30926,0,0,0,0,0.30926,0.30926,0.30926,0.30926,0,0,0,0)
    print("complete")
    print(sum1)
    return 0
def test2(x1,x2,p1,p2,alpha,beta,y1,y2,q1,q2,alpha1,beta1,x3,x4,x5,x6,p3,p4,p5,p6,y3,y4,y5,y6,q3,q4,q5,q6): # This function is used to calculate the double electron integrals!
    a1 = reexpress(x1,x2,p1,p2,alpha,beta)
    b1 = reexpress(y1,y2,q1,q2,alpha1,beta1)
    a2 = reexpress(x3,x4,p3,p4,alpha,beta)
    b2 = reexpress(y3,y4,q3,q4,alpha1,beta1)
    a3 = reexpress(x5,x6,p5,p6,alpha,beta)
    b3 = reexpress(y5,y6,q5,q6,alpha1,beta1)
    tot = 0
    for j in range(len(a1)):
        for k in range(len(b1)):
            for l in range(len(a2)):
                for m in range(len(b2)):
                    for n in range(len(a3)):
                        for o in range(len(b3)):
                            tot = tot+a1[j]*b1[k]*a2[l]*b2[m]*a3[n]*b3[o]*math.exp(-alpha*beta/(alpha+beta)*((x1-x2)**2+(x3-x4)**2+(x5-x6)**2))\
                                  *math.exp(-alpha1*beta1/(alpha1+beta1)*((y1-y2)**2+(y3-y4)**2+(y5-y6)**2))*\
                                  overlap((alpha*x1+beta*x2)/(alpha+beta),(alpha1*y1+beta1*y2)/(alpha1+beta1),j,k,alpha+beta,alpha1+beta1,
                                          Judge=False,y1=(alpha*x3+beta*x4)/(alpha+beta),y2=(alpha1*y3+beta1*y4)/(alpha1+beta1),p3=l,p4=m,
                                          z1=(alpha*x5+beta*x6)/(alpha+beta),z2=(alpha1*y5+beta1*y6)/(alpha1+beta1),p5=n,p6=o)
           # you cannot separate these three integrals, attention please!
    return tot
def H(x1,y1,z1,p1,p2,p3,x2,y2,z2,q1,q2,q3,alpha,beta,position,charge):
    Tx = (4*alpha*beta*overlap(x1,x2,p1+1,q1+1,alpha,beta,Judge=True)+p1*(-2*beta)*overlap(x1,x2,p1-1,q1+1,alpha,beta,Judge=True)+\
        p1*q1*overlap(x1,x2,p1-1,q1-1,alpha,beta,Judge=True)+(-2*alpha)*q1*overlap(x1,x2,p1+1,q1-1,alpha,beta,Judge=True))\
         *overlap(y1,y2,p2,q2,alpha,beta,Judge=True)*overlap(z1,z2,p3,q3,alpha,beta,Judge=True)
    Ty = (4*alpha*beta*overlap(y1,y2,p2+1,q2+1,alpha,beta,Judge=True)+p2*(-2*beta)*overlap(y1,y2,p2-1,q2+1,alpha,beta,Judge=True)+\
        p2*q2*overlap(y1,y2,p2-1,q2-1,alpha,beta,Judge=True)+(-2*alpha)*q2*overlap(y1,y2,p2+1,q2-1,alpha,beta,Judge=True))\
         *overlap(x1,x2,p1,q1,alpha,beta,Judge=True)*overlap(z1,z2,p3,q3,alpha,beta,Judge=True)
    Tz = (4*alpha*beta*overlap(z1,z2,p3+1,q3+1,alpha,beta,Judge=True)+p3*(-2*beta)*overlap(z1,z2,p3-1,q3+1,alpha,beta,Judge=True)+\
        p3*q3*overlap(z1,z2,p3-1,q3-1,alpha,beta,Judge=True)+(-2*alpha)*q3*overlap(z1,z2,p3+1,q3-1,alpha,beta,Judge=True))\
         *overlap(x1,x2,p1,q1,alpha,beta,Judge=True)*overlap(y1,y2,p2,q2,alpha,beta,Judge=True)
    V = 0
    a1 = reexpress(x1, x2, p1, q1, alpha, beta)
    b1 = reexpress(y1, y2, p2, q2, alpha, beta)
    c1 = reexpress(z1, z2, p3, q3, alpha, beta)
    for a in range(len(a1)):
        for b in range(len(b1)):
            for c in range(len(c1)):
                for i in range(len(position)):
                    V=V -a1[a]*b1[b]*c1[c]*charge[i]*math.exp(-alpha*beta/(alpha+beta)*((x1-x2)**2+(y1-y2)**2+(z1-z2)**2))*Pot((alpha*x1+beta*x2)/(alpha+beta),a,alpha+beta,position[i][0],0,0,(alpha*y1+beta*y2)/(alpha+beta),b,position[i][1],0,(alpha*z1+beta*z2)/(alpha+beta),c,position[i][2],0)
    return 0.5*(Tx+Ty+Tz)+V
def H1(x1,y1,z1,p1,p2,p3,x2,y2,z2,q1,q2,q3,alpha,beta,position,charge):
    Tx = (4*alpha*beta*overlap(x1,x2,p1+1,q1+1,alpha,beta,Judge=True)+p1*(-2*beta)*overlap(x1,x2,p1-1,q1+1,alpha,beta,Judge=True)+\
        p1*q1*overlap(x1,x2,p1-1,q1-1,alpha,beta,Judge=True)+(-2*alpha)*q1*overlap(x1,x2,p1+1,q1-1,alpha,beta,Judge=True))\
         *overlap(y1,y2,p2,q2,alpha,beta,Judge=True)*overlap(z1,z2,p3,q3,alpha,beta,Judge=True)
    Ty = (4*alpha*beta*overlap(y1,y2,p2+1,q2+1,alpha,beta,Judge=True)+p2*(-2*beta)*overlap(y1,y2,p2-1,q2+1,alpha,beta,Judge=True)+\
        p2*q2*overlap(y1,y2,p2-1,q2-1,alpha,beta,Judge=True)+(-2*alpha)*q2*overlap(y1,y2,p2+1,q2-1,alpha,beta,Judge=True))\
         *overlap(x1,x2,p1,q1,alpha,beta,Judge=True)*overlap(z1,z2,p3,q3,alpha,beta,Judge=True)
    Tz = (4*alpha*beta*overlap(z1,z2,p3+1,q3+1,alpha,beta,Judge=True)+p3*(-2*beta)*overlap(z1,z2,p3-1,q3+1,alpha,beta,Judge=True)+\
        p3*q3*overlap(z1,z2,p3-1,q3-1,alpha,beta,Judge=True)+(-2*alpha)*q3*overlap(z1,z2,p3+1,q3-1,alpha,beta,Judge=True))\
         *overlap(x1,x2,p1,q1,alpha,beta,Judge=True)*overlap(y1,y2,p2,q2,alpha,beta,Judge=True)
    V = 0
    for i in range(len(position)):
        V=V -charge[i]*Pot1(x1,p1,alpha,x2,q1,beta,0,y1,p2,y2,q2,z1,p3,z2,q3,position[i][0],position[i][1],position[i][2])
    return 0.5*(Tx+Ty+Tz)+V
def Pot1(x1,p1,alpha,x2,p2,beta,pt,y1,p3,y2,p4,z1,p5,z2,p6,xc,yc,zc):
    if p1 == 0 and p2 == 0 and p3 == 0 and p4 ==0 and p5 ==0 and p6==0 :
        return 2*math.pi/(alpha+beta)*F2((alpha*x1+beta*x2)/(alpha+beta),alpha+beta,xc,pt,(alpha*y1+beta*y2)/(alpha+beta),yc,(alpha*z1+beta*z2)/(alpha+beta),zc)*math.exp(-alpha*beta/(alpha+beta)*((x1-x2)**2+(y1-y2)**2+(z1-z2)**2))
    if p1 <0 or p2<0:
        return 0
    if p1>0 and p2 == 0:
        return -(x1-x2)*beta/(alpha+beta)*Pot1(x1,p1-1,alpha,x2,p2,beta,pt,y1,p3,y2,p4,z1,p5,z2,p6,xc,yc,zc)+(p1-1)/2/(alpha+beta)*Pot1(x1,p1-2,alpha,x2,p2,beta,pt,y1,p3,y2,p4,z1,p5,z2,p6,xc,yc,zc)+p2/2/(alpha+beta)*Pot1(x1,p1-1,alpha,x2,p2-1,beta,pt,y1,p3,y2,p4,z1,p5,z2,p6,xc,yc,zc)\
               +(xc-(alpha*x1+beta*x2)/(alpha+beta))*Pot1(x1,p1-1,alpha,x2,p2,beta,pt+1,y1,p3,y2,p4,z1,p5,z2,p6,xc,yc,zc)-(p1-1)/2/(alpha+beta)*Pot1(x1,p1-2,alpha,x2,p2,beta,pt+1,y1,p3,y2,p4,z1,p5,z2,p6,xc,yc,zc)-p2/2/(alpha+beta)*Pot1(x1,p1-1,alpha,x2,p2-1,beta,pt+1,y1,p3,y2,p4,z1,p5,z2,p6,xc,yc,zc)
    if p1 == 0 and p2 >0:
        return -(x2-x1)*alpha/(alpha+beta)*Pot1(x1,p1,alpha,x2,p2-1,beta,pt,y1,p3,y2,p4,z1,p5,z2,p6,xc,yc,zc)+(p2-1)/2/(alpha+beta)*Pot1(x1,p1,alpha,x2,p2-2,beta,pt,y1,p3,y2,p4,z1,p5,z2,p6,xc,yc,zc)+p1/2/(alpha+beta)*Pot1(x1,p1-1,alpha,x2,p2-1,beta,pt,y1,p3,y2,p4,z1,p5,z2,p6,xc,yc,zc)\
                +(xc-(alpha*x1+beta*x2)/(alpha+beta))*Pot1(x1,p1,alpha,x2,p2-1,beta,pt+1,y1,p3,y2,p4,z1,p5,z2,p6,xc,yc,zc)-(p2-1)/2/(alpha+beta)*Pot1(x1,p1,alpha,x2,p2-2,beta,pt+1,y1,p3,y2,p4,z1,p5,z2,p6,xc,yc,zc)-p1/2/(alpha+beta)*Pot1(x1,p1-1,alpha,x2,p2-1,beta,pt+1,y1,p3,y2,p4,z1,p5,z2,p6,xc,yc,zc)
    if p1 == 0 and p2 == 0 :
        return Pot1(y1,p3,alpha,y2,p4,beta,pt,z1,p5,z2,p6,x1,p1,x2,p2,yc,zc,xc)
    return -(x1-x2)*beta/(alpha+beta)*Pot1(x1,p1-1,alpha,x2,p2,beta,pt,y1,p3,y2,p4,z1,p5,z2,p6,xc,yc,zc)+(p1-1)/2/(alpha+beta)*Pot1(x1,p1-2,alpha,x2,p2,beta,pt,y1,p3,y2,p4,z1,p5,z2,p6,xc,yc,zc)+p2/2/(alpha+beta)*Pot1(x1,p1-1,alpha,x2,p2-1,beta,pt,y1,p3,y2,p4,z1,p5,z2,p6,xc,yc,zc)\
                +(xc-(alpha*x1+beta*x2)/(alpha+beta))*Pot1(x1,p1-1,alpha,x2,p2,beta,pt+1,y1,p3,y2,p4,z1,p5,z2,p6,xc,yc,zc)-(p1-1)/2/(alpha+beta)*Pot1(x1,p1-2,alpha,x2,p2,beta,pt+1,y1,p3,y2,p4,z1,p5,z2,p6,xc,yc,zc)-p2/2/(alpha+beta)*Pot1(x1,p1-1,alpha,x2,p2-1,beta,pt+1,y1,p3,y2,p4,z1,p5,z2,p6,xc,yc,zc)
def F2(x1,alpha,x2,pt,y1,y2,z1,z2):
    l = alpha*((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
    F0 = 0
    if l < 20:
        for i in range(0,int(2*l+8+pt/2)):
            F0 = F0 + (-l)**i/math.factorial(i)/(2*i+2*pt+1)
    else:
        if pt== 0 :
            F0 = (math.pi)**0.5/math.sqrt(l)/2
        else:
            F0 = math.factorial(2*pt-1)/math.factorial(pt-1)/4**pt*(math.pi)**0.5/l**(pt+0.5)
    return F0
def Pot(x1,p1,alpha,x2,p2,pt,y1,p3,y2,p4,z1,p5,z2,p6):
    if p1<0 :
        return 0
    if p2<0 :
        return 0
    if p1==0 and p2==0 and p3==0 and p4==0 and p5==0 and p6==0 :
        return F1(x1,p1,alpha,x2,p2,pt,y1,p3,y2,p4,z1,p5,z2,p6)
    if p1==0 and p2==0:
        return Pot(y1,p3,alpha,y2,p4,pt,z1,p5,z2,p6,x1,p1,x2,p2)
    return (p2*Pot(x1,p1-1,alpha,x2,p2-1,pt,y1,p3,y2,p4,z1,p5,z2,p6)+(p1-1)*Pot(x1,p1-2,alpha,x2,p2,pt,y1,p3,y2,p4,z1,p5,z2,p6)-2*(x1-x2)*Pot(x1,p1-1,alpha,x2,p2,pt+2,y1,p3,y2,p4,z1,p5,z2,p6))/(2+2*alpha)
def F1(x1,p1,alpha,x2,p2,pt,y1,p3,y2,p4,z1,p5,z2,p6):
    l = alpha*((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
    F0 = 0
    if l < 20:
        for i in range(0,int(2*l+8+pt/2)):
            F0 = F0 + 2*math.pi/alpha*(-l)**i/math.factorial(i)/(2*i+2*pt+1)
    else:
        if pt== 0 :
            F0 = (math.pi)**1.5/math.sqrt(l)/alpha
        else:
            F0 = math.factorial(2*pt-1)/math.factorial(pt-1)/4**pt*(math.pi)**1.5/alpha/l**(pt+0.5)*2
    return F0
def test3():
    charge = [7,1,1,1]
    position = [[0.30926,0.30926,0.30926],[2.17510,0,0],[0,2.17510,0],[0,0,2.17510]]
    orbital_coe = [[1.1714326, 0.73671471, 0.1167376]]
    orbital_exp = [[3.7805, 0.8785, 0.2857]]
    orbital_pos = [[2.17510,0,0]]
    Energy = 0
    for i in range(len(orbital_coe[0])):
        for j in range(len(orbital_coe[0])):
            Energy = Energy+orbital_coe[0][i]*orbital_coe[0][j]*H1(0.30926,0.30926,0.30926,1,0,0,0.30926,0.30926,0.30926,1,0,0,orbital_exp[0][i],orbital_exp[0][j],position,charge)
    print(Energy)
    return 0

def test5(charge,position,orbital_coe,orbital_exp,orbital_pos,orbital_typ): #This function will be treated as a part of the whole program.
    H_matrix = []
    for i in range(len(orbital_coe)):
        H_row = []
        for j in range(len(orbital_coe)):
            Energy = 0
            for k in range(len(orbital_exp[i])):
                for l in range(len(orbital_exp[j])):
                    Energy = Energy + orbital_coe[i][k] * orbital_coe[j][l] * H1(orbital_pos[i][0], orbital_pos[i][1], orbital_pos[i][2], orbital_typ[i][0], orbital_typ[i][1], orbital_typ[i][2], orbital_pos[j][0],
                                                                        orbital_pos[j][1], orbital_pos[j][2], orbital_typ[j][0], orbital_typ[j][1], orbital_typ[j][2], orbital_exp[i][k],
                                                                        orbital_exp[j][l], position, charge)
            H_row.append(Energy) #to calculate the matrix element H_{ij}
        H_matrix.append(H_row)
        print(H_row) # Thus, the calculation of H matrix is completed!!! And this function has passed the test.
    return H_matrix
def read_orbital(): #This function is merely used to store the information of atomic orbitals!
    charge = [7, 1, 1, 1]
    position = [[0.30926, 0.30926, 0.30926], [2.17510, 0, 0], [0, 2.17510, 0], [0, 0, 2.17510]]
    orbital_coe = [[3.47756711, 3.3638252, 1.04832234], [-0.193172031, 0.25835665, 0.19497785], [1.1714326, 0.73671471, 0.1167376],
                   [1.1714326, 0.73671471, 0.11673761], [1.1714326, 0.73671471, 0.11673761], [0.26359069, 0.25496339, 0.07945002],
                   [0.26359069, 0.25496339, 0.07945002], [0.26359069, 0.25496339, 0.07945002]]
    orbital_exp = [[99.9997, 18.2151, 4.9297], [3.7805, 0.8785, 0.2857], [3.7805, 0.8785, 0.2857],
                   [3.7805, 0.8785, 0.2857], [3.7805, 0.8785, 0.2857], [3.2078, 0.5843, 0.1581],
                   [3.2078, 0.5843, 0.1581], [3.2078, 0.5843, 0.1581]]
    orbital_pos = [[0.30926, 0.30926, 0.30926], [0.30926, 0.30926, 0.30926], [0.30926, 0.30926, 0.30926],
                   [0.30926, 0.30926, 0.30926], [0.30926, 0.30926, 0.30926], [2.17510, 0, 0], [0, 2.17510, 0],
                   [0, 0, 2.17510]]
    orbital_typ = [[0, 0, 0], [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 0], [0, 0, 0],
                   [0, 0, 0]]  # Attention please! All Gaussian functions of an atomic orbital share the same angular part!
    return charge,position,orbital_coe,orbital_exp,orbital_pos,orbital_typ
def overlap_matrix(orbital_coe,orbital_exp,orbital_pos,orbital_typ):
    # This function is designed to calculate the overlap integrals!
    S_matrix = []
    for i in range(len(orbital_coe)):
        S_row = []
        for j in range(len(orbital_coe)):
            S0 = 0
            for k in range(len(orbital_exp[i])):
                for l in range(len(orbital_exp[j])):
                    S0 = S0+orbital_coe[i][k]*orbital_coe[j][l]*overlap(orbital_pos[i][0],orbital_pos[j][0],orbital_typ[i][0],orbital_typ[j][0],orbital_exp[i][k],orbital_exp[j][l],Judge=True)*\
                    overlap(orbital_pos[i][1], orbital_pos[j][1], orbital_typ[i][1], orbital_typ[j][1],orbital_exp[i][k], orbital_exp[j][l], Judge=True)*overlap(orbital_pos[i][2],orbital_pos[j][2],orbital_typ[i][2],orbital_typ[j][2],orbital_exp[i][k],orbital_exp[j][l],Judge=True)
            S_row.append(S0)
        S_matrix.append(S_row)
        print(S_row)
    return S_matrix
def double_electron_matrix(orbital_coe,orbital_exp,orbital_pos,orbital_typ): #to calculate the double electron integral of ijkl
    DE_matrix = []
    for i in range(len(orbital_coe)):
        DE_1 = []
        for j in range(len(orbital_coe)):
            DE_2 = []
            for k in range(len(orbital_coe)):
                DE_3 = []
                for l in range(len(orbital_coe)):
                    DE = 0
                    for r in range(len(orbital_coe[i])):
                        for s in range(len(orbital_coe[j])):
                            for t in range(len(orbital_coe[k])):
                                for u in range(len(orbital_coe[l])):
                                    DE = DE+orbital_coe[i][r]*orbital_coe[j][s]*orbital_coe[k][t]*orbital_coe[l][u]*test2(orbital_pos[i][0],orbital_pos[j][0],orbital_typ[i][0],orbital_typ[j][0],orbital_exp[i][r],
                                         orbital_exp[j][s],orbital_pos[k][0],orbital_pos[l][0],orbital_typ[k][0],orbital_typ[l][0],orbital_exp[k][t],orbital_exp[l][u],orbital_pos[i][1],orbital_pos[j][1],orbital_pos[i][2],
                                         orbital_pos[j][2],orbital_typ[i][1],orbital_typ[j][1],orbital_typ[i][2],orbital_typ[j][2],orbital_pos[k][1],orbital_pos[l][1],orbital_pos[k][2],orbital_pos[l][2],orbital_typ[k][1],
                                         orbital_typ[l][1],orbital_typ[k][2],orbital_typ[l][2])
                                                              #x1, x2, p1, p2, alpha, beta, y1, y2, q1, q2, alpha1, beta1, x3, x4, x5, x6, p3, p4, p5, p6, y3, y4, y5, y6, q3, q4, q5, q6
                    DE_3.append(DE)
                DE_2.append(DE_3)
            DE_1.append(DE_2)
        DE_matrix.append(DE_1)
    return DE_matrix
def potential_background(charge,position):
    Potential = 0
    for i in range(len(charge)):
        for j in range(i+1,len(charge)):
            p = math.sqrt((position[i][0]-position[j][0])**2+(position[i][1]-position[j][1])**2+(position[i][2]-position[j][2])**2)
            Potential = Potential + charge[i]*charge[j]/p
    return Potential
def program():
    charge,position,orbital_coe,orbital_exp,orbital_pos,orbital_typ = read_orbital()
    H_matrix = test5(charge,position,orbital_coe,orbital_exp,orbital_pos,orbital_typ)
    S_matrix = overlap_matrix(orbital_coe,orbital_exp,orbital_pos,orbital_typ)
    DE_matrix = double_electron_matrix(orbital_coe,orbital_exp,orbital_pos,orbital_typ)
    Potential = potential_background(charge,position)
    TEST = False
    Energy = 0
    Energy1 = 0
    ei, X = LA.eigh(S_matrix)
    S_inv = LA.pinv(S_matrix)
    S_inv = np.array(S_inv)
    while TEST == False:
        Fock_ex = []
        for i in range(len(DE_matrix)):
            DE_0 = []
            for j in range(len(DE_matrix[i])):
                DE = 0
                for k in range(len(DE_matrix[i][j])):
                    for l in range(len(DE_matrix[i][j][k])):
                        for r in range(int(sum(charge)/2)):
                            DE=DE+DE_matrix[i][j][k][l]*X[k][r]*X[l][r]*2-X[k][r]*X[l][r]*DE_matrix[i][k][l][j]
                DE_0.append(DE)
            Fock_ex.append(DE_0)
        Fock_ex = np.array(Fock_ex)
        H_matrix = np.array(H_matrix)
        orb_energy,X = LA.eigh(S_inv@(H_matrix+Fock_ex))
        for i in range(int(sum(charge)/2)):
            Energy1 = Energy1+orb_energy[i]
        if abs(Energy1-Energy)<0.000001 :
            TEST = True
        Energy = Energy1
        print(Energy)
        Energy1 = 0
    print(Energy * 2 + Potential)
    print(orb_energy)
    print("Successfully completed!")
    for i in range(len(DE_matrix)):
        for j in range(len(DE_matrix[i])):
            for k in range(len(DE_matrix[i][j])):
                for l in range(len(DE_matrix[i][j][k])):
                    if abs(DE_matrix[i][j][k][l]) >0.000001 :
                        print(i,j,k,l,DE_matrix[i][j][k][l])
    print(X)
    X=[[-0.992857,-0.219941,0.000,0.000,-0.089055,-0.184147,0,0],
       [-0.034848,-0.736394,0,0,0.451106,1.27113,0,0],
       [-0.003049,0.079144,-0.415448,0.239859,-0.515084,0.324697,0.748265,-0.432010],
       [-0.003049,0.079144,0.415448,0.239859,-0.515084,0.324697,-0.748265,-0.432010],
       [-0.003049,0.079144,0,-0.479718,-0.515084,0.324697,0,0.864022],
       [0.007327,0.157274,-0.437866,0.252802,-0.124004,-0.709589,-0.871736,0.503296],
       [0.007327,0.157274,0.437866,0.252802,-0.124004,-0.709589,0.871736,0.503296],
       [0.007327,0.157274,0,-0.505604,-0.124004,-0.709589,0,-1.006595]]
    result = TEXT(DE_matrix,H_matrix,X)
    return 0

def TEXT(DE_matrix,H_matrix,X):
    Fock_ex = []
    for i in range(len(DE_matrix)):
        DE_0 = []
        for j in range(len(DE_matrix[i])):
            DE = 0
            for k in range(len(DE_matrix[i][j])):
                for l in range(len(DE_matrix[i][j][k])):
                    for r in range(5):
                        DE = DE + DE_matrix[i][j][k][l] * X[k][r] * X[l][r] * 2 - X[k][r] * X[l][r] * \
                             DE_matrix[i][k][l][j]
            DE_0.append(DE)
        Fock_ex.append(DE_0)
    Fock_ex = np.array(Fock_ex)
    H_matrix = np.matrix(H_matrix)
    print(Fock_ex+H_matrix)
    return 0
program()
#test3()
#position=[[0.30926,0.30926,0.30926],[2.17510,0.00000,0.00000],[0.00000,2.17510,0.0000],[0.0000,0.0000,2.17510]]
#wavefunction= [[[]]]

"""
3.7805 0.8785 0.2857
1.1714326 0.73671471 0.1167376
3.2078 0.5843 0.1581
0.26359069 0.25496339 0.07945002
"""
"""
3.7805 0.8785 0.2857
-0.1932 0.2584 0.1950
3.2078 0.5843 0.1581
0.26359069 0.25496339 0.07945002
"""
"""
3.2078 0.5843 0.1581
0.2635 0.2549 0.0795
3.2078 0.5843 0.1581
0.2635 0.2549 0.0795
"""
"""
99.9997 18.2151 4.9297
3.4776 3.3638 1.048
99.9997 18.2151 4.9297
3.4776 3.3638 1.048
"""
#math.sqrt(4*alpha*beta/(alpha+beta)/math.pi)