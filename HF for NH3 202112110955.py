import numpy as np
from numpy import linalg as LA
import math
import scipy.linalg
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
        return F(x1,x2,p1,p2,alpha,beta,y1,y2,p3,p4,z1,z2,p5,p6,0)
    if Judge== False and p1==0 and p2==0 :
        return overlap(y1,y2,p3,p4,alpha,beta,Judge,y1=z1,y2=z2,p3=p5,p4=p6,z1=x1,z2=x2,p5=p1,p6=p2)
    return -beta/(alpha+beta)*(x1-x2)*overlap(x1,x2,p1-1,p2,alpha,beta,Judge,y1=y1,y2=y2,p3=p3,p4=p4,z1=z1,z2=z2,p5=p5,p6=p6)\
            +(p1-1)/(2*(alpha+beta))*overlap(x1,x2,p1-2,p2,alpha,beta,Judge,y1=y1,y2=y2,p3=p3,p4=p4,z1=z1,z2=z2,p5=p5,p6=p6)\
            +p2/(2*(alpha+beta))*overlap(x1,x2,p1-1,p2-1,alpha,beta,Judge,y1=y1,y2=y2,p3=p3,p4=p4,z1=z1,z2=z2,p5=p5,p6=p6)
def overlap1(x1,x2,p1,p2,alpha,beta,y1=0,y2=0,p3=0,p4=0,z1=0,z2=0,p5=0,p6=0):
    if p1<0 or p2 <0 :
        return 0
    if p1>0 and p2 == 0:
        return 0
    if p1 == 0 and p2 >0:
        return alpha/(alpha+beta)*(x1-x2)*overlap(x1,x2,p1,p2-1,alpha,beta,y1=y1,y2=y2,p3=p3,p4=p4,z1=z1,z2=z2,p5=p5,p6=p6)\
                   +p1/(2*(alpha+beta))*overlap(x1,x2,p1-1,p2-1,alpha,beta,y1=y1,y2=y2,p3=p3,p4=p4,z1=z1,z2=z2,p5=p5,p6=p6)\
                   +(p2-1)/(2*(alpha+beta))*overlap(x1,x2,p1,p2-2,alpha,beta,y1=y1,y2=y2,p3=p3,p4=p4,z1=z1,z2=z2,p5=p5,p6=p6)
    if p1==0 and p2==0 :
        return Gauss(x1,x2,p1,p2,alpha,beta)
    if p1==0 and p2==0 and p3==0 and p4==0 and p5==0 and p6==0 :
        return F(x1,x2,p1,p2,alpha,beta,y1,y2,p3,p4,z1,z2,p5,p6)
    if p1==0 and p2==0 :
        return overlap(y1,y2,p3,p4,alpha,beta,y1=z1,y2=z2,p3=p5,p4=p6,z1=x1,z2=x2,p5=p1,p6=p2)
    return -beta/(alpha+beta)*(x1-x2)*overlap(x1,x2,p1-1,p2,alpha,beta,y1=y1,y2=y2,p3=p3,p4=p4,z1=z1,z2=z2,p5=p5,p6=p6)\
            +(p1-1)/(2*(alpha+beta))*overlap(x1,x2,p1-2,p2,alpha,beta,y1=y1,y2=y2,p3=p3,p4=p4,z1=z1,z2=z2,p5=p5,p6=p6)\
            +p2/(2*(alpha+beta))*overlap(x1,x2,p1-1,p2-1,alpha,beta,y1=y1,y2=y2,p3=p3,p4=p4,z1=z1,z2=z2,p5=p5,p6=p6)
def Gauss(x1,x2,p1,p2,alpha,beta):
    return math.sqrt(math.pi/(alpha+beta))*math.exp(-(x1-x2)**2*alpha*beta/(alpha+beta))
def F(x1,x2,p1,p2,alpha,beta,y1,y2,p3,p4,z1,z2,p5,p6,pt):
    F0=0
    l = alpha*beta/(alpha+beta)*((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
    if l <20 :
        for i in range(0,int(2*l+28+2*pt)):
            F0=F0+2*(math.pi)**2.5/alpha/beta/math.sqrt(alpha+beta)*(-l)**i/math.factorial(i)/(2*i+2*pt+1)
    else:
        if pt==0:
            F0 = math.sqrt(1/(alpha+beta)/l)*math.pi**3/alpha/beta
        else:
            F0 = math.factorial(2*pt-1)/math.factorial(pt-1)/4**pt*(math.pi)**3/l**(pt+0.5)*2/alpha/beta
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
    a = [99.9997, 18.2151, 4.9297]
    b = [3.47756711, 3.3638252, 1.04832234]
    sum1 = 0
    for j in range(len(a)):
        for l in range(len(a)):
            for m in range(len(a)):
                for n in range(len(a)):
                    sum1 = sum1+b[j]*b[l]*b[m]*b[n]*test2(0.30926,0.30926,0,0,a[j],a[l],0.30926,0.30926,0,0,a[m],a[n],0.30926,0.30926,0.30926,0.30926,0,0,0,0,0.30926,0.30926,0.30926,0.30926,0,0,0,0)
    print("complete")
    print(sum1)
    return 0
def test2(x1,x2,p1,p2,alpha,beta,y1,y2,q1,q2,alpha1,beta1,x3,x4,x5,x6,p3,p4,p5,p6,y3,y4,y5,y6,q3,q4,q5,q6,pt): # This function is used to calculate the double electron integrals!
    tot = 0
    if p1<0 or p2<0 or q1<0 or q2<0:
        return 0
    if p1== 0 and p2 ==0 and p3 == 0 and p4 == 0 and p5== 0 and p6 ==0 and q1 == 0 and q2 ==0 and q3 ==0 and q4 ==0 and q5 ==0 and q6==0:
        return math.exp(-alpha*beta/(alpha+beta)*((x1-x2)**2+(x3-x4)**2+(x5-x6)**2))*math.exp(-alpha1*beta1/(alpha1+beta1)*((y1-y2)**2+(y3-y4)**2+(y5-y6)**2))*F((alpha*x1+beta*x2)/(alpha+beta),(alpha1*y1+beta1*y2)/(alpha1+beta1),p1,p2,alpha+beta,alpha1+beta1,(alpha*x3+beta*x4)/(alpha+beta),(alpha1*y3+beta1*y4)/(alpha1+beta1),p3,p4,(alpha*x5+beta*x6)/(alpha+beta),(alpha1*y5+beta1*y6)/(alpha1+beta1),p5,p6,pt)
    if p1>0 and p2==0:
        return beta/(alpha+beta)*(x2-x1)*test2(x1,x2,p1-1,p2,alpha,beta,y1,y2,q1,q2,alpha1,beta1,x3,x4,x5,x6,p3,p4,p5,p6,y3,y4,y5,y6,q3,q4,q5,q6,pt)-(alpha1+beta1)/(alpha+beta+alpha1+beta1)*((alpha*x1+beta*x2)/(alpha+beta)-(alpha1*y1+beta1*y2)/(alpha1+beta1))*test2(x1,x2,p1-1,p2,alpha,beta,y1,y2,q1,q2,alpha1,beta1,x3,x4,x5,x6,p3,p4,p5,p6,y3,y4,y5,y6,q3,q4,q5,q6,pt+1)\
               +(p1-1)/2/(alpha+beta)*(test2(x1,x2,p1-2,p2,alpha,beta,y1,y2,q1,q2,alpha1,beta1,x3,x4,x5,x6,p3,p4,p5,p6,y3,y4,y5,y6,q3,q4,q5,q6,pt)-(alpha1+beta1)/(alpha1+beta1+alpha+beta)*test2(x1,x2,p1-2,p2,alpha,beta,y1,y2,q1,q2,alpha1,beta1,x3,x4,x5,x6,p3,p4,p5,p6,y3,y4,y5,y6,q3,q4,q5,q6,pt+1))+p2/2/(alpha+beta)*(test2(x1,x2,p1-1,p2-1,alpha,beta,y1,y2,q1,q2,alpha1,beta1,x3,x4,x5,x6,p3,p4,p5,p6,y3,y4,y5,y6,q3,q4,q5,q6,pt)\
               -(alpha1+beta1)/(alpha+beta+alpha1+beta1)*test2(x1,x2,p1-1,p2-1,alpha,beta,y1,y2,q1,q2,alpha1,beta1,x3,x4,x5,x6,p3,p4,p5,p6,y3,y4,y5,y6,q3,q4,q5,q6,pt+1))+q1/2/(alpha+beta+alpha1+beta1)*test2(x1,x2,p1-1,p2,alpha,beta,y1,y2,q1-1,q2,alpha1,beta1,x3,x4,x5,x6,p3,p4,p5,p6,y3,y4,y5,y6,q3,q4,q5,q6,pt+1)+q2/2/(alpha+beta+alpha1+beta1)*test2(x1,x2,p1-1,p2,alpha,beta,y1,y2,q1,q2-1,alpha1,beta1,x3,x4,x5,x6,p3,p4,p5,p6,y3,y4,y5,y6,q3,q4,q5,q6,pt+1)
    if p1== 0 and p2>0:
        return alpha/(alpha+beta)*(x1-x2)*test2(x1,x2,p1,p2-1,alpha,beta,y1,y2,q1,q2,alpha1,beta1,x3,x4,x5,x6,p3,p4,p5,p6,y3,y4,y5,y6,q3,q4,q5,q6,pt)-(alpha1+beta1)/(alpha+beta+alpha1+beta1)*((alpha*x1+beta*x2)/(alpha+beta)-(alpha1*y1+beta1*y2)/(alpha1+beta1))*test2(x1,x2,p1,p2-1,alpha,beta,y1,y2,q1,q2,alpha1,beta1,x3,x4,x5,x6,p3,p4,p5,p6,y3,y4,y5,y6,q3,q4,q5,q6,pt+1)\
               +(p2-1)/2/(alpha+beta)*(test2(x1,x2,p1,p2-2,alpha,beta,y1,y2,q1,q2,alpha1,beta1,x3,x4,x5,x6,p3,p4,p5,p6,y3,y4,y5,y6,q3,q4,q5,q6,pt)-(alpha1+beta1)/(alpha1+beta1+alpha+beta)*test2(x1,x2,p1,p2-2,alpha,beta,y1,y2,q1,q2,alpha1,beta1,x3,x4,x5,x6,p3,p4,p5,p6,y3,y4,y5,y6,q3,q4,q5,q6,pt+1))+p1/2/(alpha+beta)*(test2(x1,x2,p1-1,p2-1,alpha,beta,y1,y2,q1,q2,alpha1,beta1,x3,x4,x5,x6,p3,p4,p5,p6,y3,y4,y5,y6,q3,q4,q5,q6,pt)\
               -(alpha1+beta1)/(alpha+beta+alpha1+beta1)*test2(x1,x2,p1-1,p2-1,alpha,beta,y1,y2,q1,q2,alpha1,beta1,x3,x4,x5,x6,p3,p4,p5,p6,y3,y4,y5,y6,q3,q4,q5,q6,pt+1))+q1/2/(alpha+beta+alpha1+beta1)*test2(x1,x2,p1,p2-1,alpha,beta,y1,y2,q1-1,q2,alpha1,beta1,x3,x4,x5,x6,p3,p4,p5,p6,y3,y4,y5,y6,q3,q4,q5,q6,pt+1)+q2/2/(alpha+beta+alpha1+beta1)*test2(x1,x2,p1,p2-1,alpha,beta,y1,y2,q1,q2-1,alpha1,beta1,x3,x4,x5,x6,p3,p4,p5,p6,y3,y4,y5,y6,q3,q4,q5,q6,pt+1)
    if p1==0 and p2 == 0 and q1 == 0 and q2 == 0:
        return test2(x3,x4,p3,p4,alpha,beta,y3,y4,q3,q4,alpha1,beta1,x5,x6,x1,x2,p5,p6,p1,p2,y5,y6,y1,y2,q5,q6,q1,q2,pt)
    if p1==0 and p2 == 0:
        return test2(y1,y2,q1,q2,alpha1,beta1,x1,x2,p1,p2,alpha,beta,y3,y4,y5,y6,q3,q4,q5,q6,x3,x4,x5,x6,p3,p4,p5,p6,pt)
    return beta/(alpha+beta)*(x2-x1)*test2(x1,x2,p1-1,p2,alpha,beta,y1,y2,q1,q2,alpha1,beta1,x3,x4,x5,x6,p3,p4,p5,p6,y3,y4,y5,y6,q3,q4,q5,q6,pt)-(alpha1+beta1)/(alpha+beta+alpha1+beta1)*((alpha*x1+beta*x2)/(alpha+beta)-(alpha1*y1+beta1*y2)/(alpha1+beta1))*test2(x1,x2,p1-1,p2,alpha,beta,y1,y2,q1,q2,alpha1,beta1,x3,x4,x5,x6,p3,p4,p5,p6,y3,y4,y5,y6,q3,q4,q5,q6,pt+1)\
               +(p1-1)/2/(alpha+beta)*(test2(x1,x2,p1-2,p2,alpha,beta,y1,y2,q1,q2,alpha1,beta1,x3,x4,x5,x6,p3,p4,p5,p6,y3,y4,y5,y6,q3,q4,q5,q6,pt)-(alpha1+beta1)/(alpha1+beta1+alpha+beta)*test2(x1,x2,p1-2,p2,alpha,beta,y1,y2,q1,q2,alpha1,beta1,x3,x4,x5,x6,p3,p4,p5,p6,y3,y4,y5,y6,q3,q4,q5,q6,pt+1))+p2/2/(alpha+beta)*(test2(x1,x2,p1-1,p2-1,alpha,beta,y1,y2,q1,q2,alpha1,beta1,x3,x4,x5,x6,p3,p4,p5,p6,y3,y4,y5,y6,q3,q4,q5,q6,pt)\
               -(alpha1+beta1)/(alpha+beta+alpha1+beta1)*test2(x1,x2,p1-1,p2-1,alpha,beta,y1,y2,q1,q2,alpha1,beta1,x3,x4,x5,x6,p3,p4,p5,p6,y3,y4,y5,y6,q3,q4,q5,q6,pt+1))+q1/2/(alpha+beta+alpha1+beta1)*test2(x1,x2,p1-1,p2,alpha,beta,y1,y2,q1-1,q2,alpha1,beta1,x3,x4,x5,x6,p3,p4,p5,p6,y3,y4,y5,y6,q3,q4,q5,q6,pt+1)+q2/2/(alpha+beta+alpha1+beta1)*test2(x1,x2,p1-1,p2,alpha,beta,y1,y2,q1,q2-1,alpha1,beta1,x3,x4,x5,x6,p3,p4,p5,p6,y3,y4,y5,y6,q3,q4,q5,q6,pt+1)
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
def T(x1,y1,z1,p1,p2,p3,x2,y2,z2,q1,q2,q3,alpha,beta,position,charge):
    Tx = (4*alpha*beta*overlap(x1,x2,p1+1,q1+1,alpha,beta,Judge=True)+p1*(-2*beta)*overlap(x1,x2,p1-1,q1+1,alpha,beta,Judge=True)+
        p1*q1*overlap(x1,x2,p1-1,q1-1,alpha,beta,Judge=True)+(-2*alpha)*q1*overlap(x1,x2,p1+1,q1-1,alpha,beta,Judge=True))\
         *overlap(y1,y2,p2,q2,alpha,beta,Judge=True)*overlap(z1,z2,p3,q3,alpha,beta,Judge=True)
    Ty = (4*alpha*beta*overlap(y1,y2,p2+1,q2+1,alpha,beta,Judge=True)+p2*(-2*beta)*overlap(y1,y2,p2-1,q2+1,alpha,beta,Judge=True)+
        p2*q2*overlap(y1,y2,p2-1,q2-1,alpha,beta,Judge=True)+(-2*alpha)*q2*overlap(y1,y2,p2+1,q2-1,alpha,beta,Judge=True))\
         *overlap(x1,x2,p1,q1,alpha,beta,Judge=True)*overlap(z1,z2,p3,q3,alpha,beta,Judge=True)
    Tz = (4*alpha*beta*overlap(z1,z2,p3+1,q3+1,alpha,beta,Judge=True)+p3*(-2*beta)*overlap(z1,z2,p3-1,q3+1,alpha,beta,Judge=True)+
        p3*q3*overlap(z1,z2,p3-1,q3-1,alpha,beta,Judge=True)+(-2*alpha)*q3*overlap(z1,z2,p3+1,q3-1,alpha,beta,Judge=True))\
         *overlap(x1,x2,p1,q1,alpha,beta,Judge=True)*overlap(y1,y2,p2,q2,alpha,beta,Judge=True)
    return 0.5*(Tx+Ty+Tz)
def H1(x1,y1,z1,p1,p2,p3,x2,y2,z2,q1,q2,q3,alpha,beta,position,charge):
    Tx = (4*alpha*beta*overlap(x1,x2,p1+1,q1+1,alpha,beta,Judge=True)+p1*(-2*beta)*overlap(x1,x2,p1-1,q1+1,alpha,beta,Judge=True)+
        p1*q1*overlap(x1,x2,p1-1,q1-1,alpha,beta,Judge=True)+(-2*alpha)*q1*overlap(x1,x2,p1+1,q1-1,alpha,beta,Judge=True))\
         *overlap(y1,y2,p2,q2,alpha,beta,Judge=True)*overlap(z1,z2,p3,q3,alpha,beta,Judge=True)
    Ty = (4*alpha*beta*overlap(y1,y2,p2+1,q2+1,alpha,beta,Judge=True)+p2*(-2*beta)*overlap(y1,y2,p2-1,q2+1,alpha,beta,Judge=True)+
        p2*q2*overlap(y1,y2,p2-1,q2-1,alpha,beta,Judge=True)+(-2*alpha)*q2*overlap(y1,y2,p2+1,q2-1,alpha,beta,Judge=True))\
         *overlap(x1,x2,p1,q1,alpha,beta,Judge=True)*overlap(z1,z2,p3,q3,alpha,beta,Judge=True)
    Tz = (4*alpha*beta*overlap(z1,z2,p3+1,q3+1,alpha,beta,Judge=True)+p3*(-2*beta)*overlap(z1,z2,p3-1,q3+1,alpha,beta,Judge=True)+
        p3*q3*overlap(z1,z2,p3-1,q3-1,alpha,beta,Judge=True)+(-2*alpha)*q3*overlap(z1,z2,p3+1,q3-1,alpha,beta,Judge=True))\
         *overlap(x1,x2,p1,q1,alpha,beta,Judge=True)*overlap(y1,y2,p2,q2,alpha,beta,Judge=True)
    V = 0
    for i in range(len(position)):
        V = V - charge[i]*Pot1(x1,p1,alpha,x2,q1,beta,0,y1,p2,y2,q2,z1,p3,z2,q3,position[i][0],position[i][1],position[i][2])
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
        for i in range(0,int(2*l+38+2*pt)):
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
    orbital_coe = [[1.1714326, 0.73671471, 0.1167376],[0.26359069, 0.25496339, 0.07945002]]
    orbital_exp = [[3.7805, 0.8785, 0.2857],[3.2078, 0.5843, 0.1581]]
    orbital_pos = [[2.17510,0,0]]
    Energy1 = 0
    Energy2 = 0
    for i in range(len(orbital_coe[0])):
        for j in range(len(orbital_coe[0])):
            Energy1 = Energy1+orbital_coe[0][i]*orbital_coe[1][j]*T(0.30926,0.30926,0.30926,1,0,0,2.17510,0,0,0,0,0,orbital_exp[0][i],orbital_exp[1][j],position,charge)
            Energy2 = Energy2+orbital_coe[1][i]*orbital_coe[0][j]*T(2.17510,0,0,0,0,0,0.30926,0.30926,0.30926,1,0,0,orbital_exp[1][i],orbital_exp[0][j],position,charge)
    print(Energy1,Energy2)
    return 0

def test5(charge,position,orbital_coe,orbital_exp,orbital_pos,orbital_typ): #This function will be treated as a part of the whole program.
    H_matrix = []
    T_matrix = []
    for i in range(len(orbital_coe)):
        H_row = []
        T_row = []
        for j in range(len(orbital_coe)):
            Energy = 0
            T_energy = 0
            for k in range(len(orbital_exp[i])):
                for l in range(len(orbital_exp[j])):
                    Energy = Energy + orbital_coe[i][k] * orbital_coe[j][l] * H1(orbital_pos[i][0], orbital_pos[i][1], orbital_pos[i][2], orbital_typ[i][0], orbital_typ[i][1], orbital_typ[i][2], orbital_pos[j][0],
                                                                        orbital_pos[j][1], orbital_pos[j][2], orbital_typ[j][0], orbital_typ[j][1], orbital_typ[j][2], orbital_exp[i][k],
                                                                        orbital_exp[j][l], position, charge)
                    T_energy = T_energy+ orbital_coe[i][k] * orbital_coe[j][l] * T(orbital_pos[i][0], orbital_pos[i][1], orbital_pos[i][2], orbital_typ[i][0], orbital_typ[i][1], orbital_typ[i][2], orbital_pos[j][0],
                                                                        orbital_pos[j][1], orbital_pos[j][2], orbital_typ[j][0], orbital_typ[j][1], orbital_typ[j][2], orbital_exp[i][k],
                                                                        orbital_exp[j][l], position, charge)
            H_row.append(Energy) #to calculate the matrix element H_{ij}
            T_row.append(T_energy)
        H_matrix.append(H_row)
        T_matrix.append(T_row)
        print(H_row)
    return H_matrix
def read_orbital(): #This function is merely used to store the information of atomic orbitals!
    charge = [7, 1, 1, 1]
    position = [[0, 0, 0.21596924], [0, 1.775367459, -0.5039294864], [1.537513301, -0.8876837296, -0.5039294864], [-1.537513301, -0.8876837296, -0.5039294864]]
    orbital_coe = [[3.454880726, 3.341410226, 1.041371618], [-0.193163861, 0.2583722691, 0.1949974057], [1.171552502, 0.7367036516, 0.1167062357],
                   [1.171552502, 0.7367036516, 0.1167062357], [1.171552502, 0.7367036516, 0.1167062357], [0.2769343551, 0.2678388516, 0.08347367113],
                   [0.2769343551, 0.2678388516, 0.08347367113], [0.2769343551, 0.2678388516, 0.08347367113]]
    orbital_exp = [[99.10616896, 18.05231239, 4.885660238], [3.780455879, 0.8784966449, 0.2857143744], [3.780455879, 0.8784966449, 0.2857143744],
                   [3.780455879, 0.8784966449, 0.2857143744], [3.780455879, 0.8784966449, 0.2857143744], [3.425250914, 0.6239137298, 0.1688554040],
                   [3.425250914, 0.6239137298, 0.1688554040], [3.425250914, 0.6239137298, 0.1688554040]]
    orbital_pos = [[0, 0, 0.21596924], [0, 0, 0.21596924], [0, 0, 0.21596924],
                   [0, 0, 0.21596924], [0, 0, 0.21596924], [0, 1.775367459, -0.5039294864], [1.537513301, -0.8876837296, -0.5039294864],
                   [-1.537513301, -0.8876837296, -0.5039294864]]
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
                                         orbital_typ[l][1],orbital_typ[k][2],orbital_typ[l][2],0)
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
    step = 0
    ei, X = LA.eigh(S_matrix)
    S_inv = LA.pinv(S_matrix)
    #S_sqr = scipy.linalg.sqrtm(S_inv)
    A,B = LA.eigh(S_inv)
    A = np.diag(A**0.5)
    S_sqr = B@A@LA.inv(B)
    while TEST == False:
        Energy1 = 0
        step += 1
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

        orb_energy,X = LA.eigh(S_sqr.T@(H_matrix+Fock_ex)@S_sqr)
        X = S_sqr@X
        for i in range(int(sum(charge)/2)):
            Energy1 = Energy1+orb_energy[i]*2
        for i in range(len(DE_matrix)):
            for j in range(len(DE_matrix[i])):
                for s in range(int(sum(charge)/2)):
                    for k in range(len(DE_matrix[i][j])):
                        for l in range(len(DE_matrix[i][j][k])):
                            for r in range(int(sum(charge)/2)):
                                Energy1 = Energy1-(DE_matrix[i][j][k][l]*X[k][r]*X[l][r]*2-X[k][r]*X[l][r]*DE_matrix[i][k][l][j])*X[i][s]*X[j][s]
        if abs(Energy1-Energy)<0.000001 :
            TEST = True
        Energy = Energy1
        if step%5 == 0:
            print(Energy)
    print(Energy  + Potential)
    print(orb_energy)
    print(H_matrix+Fock_ex)
    print("Successfully completed!")
    print(X)
    for i in range(len(DE_matrix)):
        for j in range(len(DE_matrix[i])):
            for l in range(len(DE_matrix[i][j])):
                for m in range(len(DE_matrix[i][j][l])):
                    if abs(DE_matrix[i][j][l][m]) >0.00000001 :
                        print(i+1,j+1,l+1,m+1,DE_matrix[i][j][l][m]/2)
    return 0


if __name__ == "__main__":
    program()
#test1()
#scipy.linalg.sqrtm