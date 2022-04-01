from matplotlib.pyplot import axis
import numpy as np
import scipy
from math import *
#类和函数
class pole:
    E=0
    A=0
    L=0
    theta=0
    phi=0
    k=0
    enable=True
    def __init__(self,E,A,L,theta,phi,enable=True):
        self.theta=theta
        self.E=E
        self.A=A
        self.L=L
        self.k=E*A/L
        self.enable=enable
        self.phi=phi
def GenMatrix(poles):
    ret=np.zeros((3,3),dtype=np.float64)
    for p in poles:
        if p.enable==False:
            continue
        [xi,yi,zi]=[cos(p.phi)*cos(p.theta),cos(p.phi)*sin(p.theta),sin(p.phi)]
        ret=ret+np.array([p.k],dtype=np.float64)*np.array([[xi**2,xi*yi,xi*zi],[xi*yi,yi**2,yi*zi],[xi*zi,yi*zi,zi**2]],dtype=np.float64)
    return ret
def CalR(F,mat,error=0.01):
    mat_inv=np.linalg.inv(mat)
    cond=np.linalg.cond(mat,np.inf)
    R_error=cond*error/(1-cond*error)
    #print("cond(mat)={}".format(cond))
    ret=np.dot(mat_inv,F)
    return [ret,R_error]
##主进程
#初始化
pole1=pole(2,1,1,0,0)
pole2=pole(1,1,1,pi/2,0)
pole3=pole(1,1,1,0,pi/2)
pole4=pole(1,1,1,pi/4,pi/3,enable=True)
pole5=pole(1,1,1,pi/4,0,enable=True)

poles=(pole1,pole2,pole3,pole4,pole5)
F=np.array([1,1,1])
#打印输入
cnt=0
for p in poles:
    if p.enable==False:
        continue
    print("k{0}={1:.3f},theta{0}={2:.1f}°,phi={3:.1f}".format(cnt+1,p.k,p.theta*180/pi,p.phi*180/pi))
    cnt=cnt+1
F_r=[sqrt(F[0]**2+F[1]**2+F[2]**2),atan(F[2]/sqrt(F[0]**2+F[1]**2)),atan(F[1]/F[0])]
print("F={}".format(F))
print("F={0:.3f}，theta={2:.3f}°,phi={1:.3f}°".format(F_r[0],F_r[1]*180/pi,F_r[2]*180/pi))
#计算
mat=GenMatrix(poles)
error=0.02
[R,R_error]=CalR(F,mat,error)
R_rho=sqrt(R[0]**2+R[1]**2+R[2]**2)
R_theta=atan(R[1]/R[0])*180/pi
R_phi=atan(F[2]/sqrt(F[1]**2+F[0]**2))*180/pi
#打印输出
print("mat=\n{}".format(mat))
print("R=",end="")
print(R)
print("Rrho={0:.3f},R_theta={2:.1f}°,R_phi={1:.1f}°".format(R_rho,R_phi,R_theta))
print("relative error of F as mat {1}% error:{0:.2f}%".format(R_error*100,error*100))