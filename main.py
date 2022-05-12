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
    k=0
    enable=True
    def __init__(self,E,A,L,theta,enable=True):
        self.theta=theta
        self.E=E
        self.A=A
        self.L=L
        self.k=E*A/L
        self.enable=enable
def GenMatrix(poles):
    sum1=0.0
    sum2=0.0
    sum3=0.0
    for p in poles:
        if p.enable==False:
            continue
        sum1=sum1+p.k*cos(p.theta)**2
        sum2=sum2+p.k*cos(p.theta)*sin(p.theta)
        sum3=sum3+p.k*sin(p.theta)**2
    ret=np.array([[sum1,sum2],[sum2,sum3]],dtype=np.float64)
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
pole1=pole(1,1,1,pi/2)
pole2=pole(1,1,1,1.01*pi/4,False)
pole3=pole(1,1,1,3*pi/4,False)
pole4=pole(1,1,1,pi)
pole5=pole(1,1,1,pi/4,enable=False)

poles=(pole1,pole2,pole3,pole4,pole5)
F=np.array([cos(pi/2),sin(pi/2)])
#打印输入
cnt=0
for p in poles:
    if p.enable==False:
        continue
    print("k{2}={0:.3f},theta{2}={1:.1f}°".format(p.k,p.theta*180/pi,cnt+1))
    cnt=cnt+1
F_r=[sqrt(F[0]**2+F[1]**2),atan(F[1]/F[0])*180/pi]
print("F={0:.3f},{1:.3f}°".format(F_r[0],F_r[1]))
#计算
mat=GenMatrix(poles)
error=0.02
[R,R_error]=CalR(F,mat,error)
R_rho=sqrt(R[0]**2+R[1]**2)
R_phi=atan(R[1]/R[0])*180/pi
#打印输出
print("mat=\n{}".format(mat))
print("R=",end="")
print(R)
print("Rrho={0:.3f},R_phi={1:.1f}°".format(R_rho,R_phi))
print("relative error of F as mat {1}% error:{0:.2f}%".format(R_error*100,error*100))