import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import art3d
import numpy as np
import cmath as cm
from numpy import sin,cos,exp,sinh,cosh,pi,arctan,sqrt
from scipy.special import jv, jvp, hankel1, h1vp, hankel2, h2vp,lpmv,spherical_jn,spherical_yn
#直角坐标系
k0=50
dc = 0.3
def cartesian(ratio,a,x,V0):
    xmin, xmax = 0, a
    x1 = np.linspace(-dc-xmin,xmin,1000)
    x2 = np.linspace(xmin,xmax,1000)
    x3 = np.linspace(xmax,xmax+dc,1000)
    j=1j
    k1, k2 = k0*np.sqrt(ratio)*np.sqrt(V0), k0*cm.sqrt(1-1/ratio)
    #T = abs(k1-k2)*(k1+k2)*sin(a*k2)/(2*j*k1*k2*cos(a*k2)+(k1**2+k2**2)*sin(a*k2))
    #R = abs(2*j*exp(-j*a*k1)*k1*k2/(2*j*k1*k2*cos(a*k2)+(k1**2+k2**2)*sin(a*k2)))
    A1 = 1
    
    A2, C1 = (k1**2-k2**2)*sin(a*k2)/(2*j*k1*k2*cos(a*k2)+(k1**2+k2**2)*sin(a*k2)), 2*j*exp(-j*a*k1)*k1*k2/(2*j*k1*k2*cos(a*k2)+(k1**2+k2**2)*sin(k2*a))
    B1, B2 = -2*k1*(k1+k2)/(exp(2*j*a*k2)*(k1-k2)**2-(k1+k2)**2), 2*exp(2*j*a*k2)*k1*(k1-k2)/(exp(2*j*a*k2)*(k1-k2)**2-(k1+k2)**2)
    psai11 = A1*np.exp(j*k1*x)
    psai12 = A2*np.exp(-j*k1*x)
    psai21 = B1*np.exp(j*k2*x)
    psai22 = B2*np.exp(-j*k2*x)
    psai31 = C1*np.exp(j*k1*x)
    psai32 = 0
    
    if x in x1:
        return psai11, psai12
    elif x in x2:
        return psai21, psai22
    else:
        return psai31, psai32
    
def polar(ratio,a,r,fai,V0,m):
    rmax = a
    r2 = np.linspace(0,rmax,100)
    j=1j
    k1, k2 = k0*np.sqrt(ratio)*np.sqrt(V0), k0*cm.sqrt(1-1/ratio)
    A1, D1 = -(-k1*hankel2(m,a*k1)*h1vp(m,a*k1)+k1*hankel1(m,a*k1)*h2vp(m,a*k1))/(k2*hankel2(m,a*k1)*jvp(m,a*k2)-k1*jv(m,a*k2)*h2vp(m,a*k1)), -(k2*hankel1(m,a*k1)*jvp(m,a*k2)-k1*jv(m,a*k2)*h1vp(m,a*k1))/(k2*hankel2(m,a*k1)*jvp(m,a*k2)-k1*jv(m,a*k2)*h2vp(m,a*k1))
    psai1 = (1*hankel1(m,k1*r)+D1*hankel2(m,k1*r))*exp(j*m*fai)
    psai2 = A1*jv(m,k2*r)*exp(j*m*fai)
    if r in r2:
        return psai2
    else:
        return psai1
    
def cylindrical(ratio,a,h,r,fai,V0,m):
    zmax = h 
    rmax = a
    z = np.linspace(-dc,zmax+dc,100)
    z1 = np.linspace(-dc-zmax,0,100)
    z2 = np.linspace(0,zmax,100)
    z3 = np.linspace(zmax,zmax+dc,100)
    r2 = np.linspace(0,rmax,100)
    psair = polar(ratio,a,r,fai,V0,m)
    psaiz1, psaiz2 = cartesian(ratio,a,h,V0)
    
    p1, p2 = psair*psaiz1, psair*psaiz2
    
    return p1, p2

def Y(m,l,theta,fai):
    j=1j
    Y = lpmv(m,l,cos(theta))*exp(j*m*fai)
    return Y

sj = spherical_jn
sy = spherical_yn
def sh1(l,x):
    j=1j
    return sj(l,x)+j*sy(l,x)
def sh2(l,x):
    j=1j
    return sj(1,x)-j*sy(l,x)
def dsj(l,x):
    return sj(l,x,derivative=True)
def dsy(l,x):
    return sy(l,x,derivative=True)
def dsh1(l,x):
    j=1j
    return dsj(l,x)+j*dsy(l,x)
def dsh2(l,x):
    j=1j
    return dsj(l,x)-j*dsy(l,x)
def spherical(ratio,a,r,theta,fai,V0,l,m):
    rmax = a
    r2 = np.linspace(0,rmax,100)
    j=1j
    k1, k2 = k1, k2 = k0*np.sqrt(ratio)*np.sqrt(V0), k0*cm.sqrt(1-1/ratio)
    A,C = k1*(dsh2(l,a*k1)*sh1(l,a*k1)-dsh1(l,a*k1)*sh2(l,a*k1))/(k1*dsh2(l,a*k1)*sj(l,a*k2)-k2*dsj(l,a*k2)*sh2(l,a*k1)), (-k1*dsh1(l,a*k1)*sj(l,a*k2)+k2*dsj(l,a*k2)*sh1(l,a*k1))/(k1*dsh2(l,a*k1)*sj(l,a*k2)-k2*dsj(l,a*k2)*sh2(l,a*k1))
    psai1 = (1*sh1(l,k1*r)+C*sh2(l,k1*r))*Y(m,l,theta,fai)
    psai2 = A*sj(l,k2*r)*Y(m,l,theta,fai)
    if r in r2:
        return psai2
    else:
        return psai1

def transmission_matrix(ratio,V0,a):
    j=1j
    k1, k2 = k0*np.sqrt(ratio)*np.sqrt(V0), k0*cm.sqrt(1-1/ratio)
    B,F = (k1-k2)/(k1+k2)*exp(-2j*k2*a), 2*k1/(k1+k2)*exp(-j*k2*a)
    M = np.array([[F,-B],[0,1]])
    return M

def multi_car(ratio,a,V):
    j=1j
    k1, k2 = k0*np.sqrt(ratio)*np.sqrt(V[0]), k0*cm.sqrt(1-1/ratio)
    ncar = len(V)
    M1 = transmission_matrix(ratio,V[0],a)
    F1 = 2*j*exp(-j*a*k1)*k1*k2/(2*j*k1*k2*cos(a*k2)+(k1**2+k2**2)*sin(k2*a))
    M_total = np.linalg.matrix_power(M1,ncar)

    F_total = M_total[0,0]
    R = []
    for i in range(ncar):
        R.append(M1[0,1]*(F_total)**(i+1)/(M1[0,0])**(i+1))
    return R

#plot
def car_plot(V0,ratio,a):
    V0, ratio, a =float(V0), float(ratio), float(a)
    xmin, xmax = 0, a
    x = np.linspace(-xmin-dc,xmax+dc, 300)
    x1 = np.linspace(-dc,xmin,1000)
    x2 = np.linspace(xmin,xmax,1000)
    x3 = np.linspace(xmax,xmax+dc,1000)
    y11, y12 = cartesian(ratio,a,x1,V0)
    y21, y22 = cartesian(ratio,a,x2,V0)
    y31, y32 = cartesian(ratio,a,x3,V0)
    y1, y2 = [], []
    for i in range(len(y11)):
        y1.append(y11[i]+y12[i])
        y2.append(y21[i]+y22[i])
    yV0 = np.zeros(len(x))
    yV0[(x>xmin) & (x<xmax)] = V0

    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)

    plt.rcParams['figure.figsize'] = (8.0,6.0)
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    ax1.plot(x1,y11,label='入射波')
    ax1.plot(x1,y12,label='反射波')
    ax1.plot(x1,y1,label='复合波1')
    ax1.plot(x2,y2,label='复合波2')
    ax1.plot(x3,y31,label='透射波')
    ax1.plot(x,yV0,label='势')
    ax1.set_xticks(np.arange(-xmin-dc,xmax+dc,0.2))
    ax1.legend()
    plt.title('一维情况下波函数幅值分布')
    return fig

def polar_plot(V0,ratio,a,m):
    #V0, ratio, a, m =1,0.2,0.2,1
    V0, ratio, a, m =float(V0), float(ratio), float(a), float(m)
    rmin, rmax, faimax = 0, a, 2*pi

    Z = np.zeros((200,200))
    XX = np.linspace(-rmax-dc,rmax+dc,200)
    x = list(XX)
    YY = np.linspace(-rmax-dc,rmax+dc,200)
    y = list(YY)
    for i in range(len(x)):
        if x[i] > 0:
            for j in range(len(y)):
                fai = arctan(y[j]/x[i])
                r = sqrt(x[i]**2+y[j]**2)
                Z[j][i] = polar(ratio,a,r,fai,V0,m)
        elif x[i]<0:
            for j in range(len(y)):
                fai = arctan(y[j]/x[i])
                r = sqrt(x[i]**2+y[j]**2)
                Z[j][i] = polar(ratio,a,r,fai,V0,m) 
        else:
            for j in range(len(y)):
                if y[j] >= 0:
                    fai = pi/2
                else:
                    fai = -pi/2
                r = sqrt(x[i]**2+y[j]**2)
                Z[j][i] = polar(ratio,a,r,fai,V0,m)
    X,Y = np.meshgrid(XX,YY)
    # 设置中文字体
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False

    # 创建图形和坐标轴
    fig = plt.figure(figsize=(10, 5))  # 调整图形大小

    # 使用 gridspec_kw 自定义子图布局
    gs = fig.add_gridspec(1, 2, width_ratios=[2, 1.2])  # 调整比例

    # 默认视角
    ax1 = fig.add_subplot(gs[0], projection='3d')  # 占据两列
    u = np.linspace(0, 2 * np.pi, 50)
    h = np.linspace(0, V0, 20)
    xp = np.outer(rmax * np.cos(u), np.ones(len(h)))
    yp = np.outer(rmax * np.sin(u), np.ones(len(h)))
    zp = np.outer(np.ones(len(u)), h)

    norm1 = Normalize(-1,1)
    # 使用更柔和的颜色
    ax1.plot_surface(xp, yp, zp, color='black', alpha=0.2)  # 调整圆柱颜色和透明度
    circle = Circle((0, 0), rmax, color='black', alpha=0.2)
    ax1.add_patch(circle)
    art3d.pathpatch_2d_to_3d(circle, z=V0, zdir="z")
    ax1.plot_surface(X, Y, Z, cmap='coolwarm',norm=norm1)  # 更改 colormap
    ax1.view_init(elev=37, azim=41)  # 仰角 37 度，方位角 41 度
    ax1.set_zlim(-1, 1)
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_title("默认视角")

    # 俯视
    ax2 = fig.add_subplot(gs[1], projection='3d')  # 位于第二行第二列
    ax2.plot_surface(xp, yp, zp, color='black', alpha=0.2)  # 调整圆柱颜色和透明度
    circle = Circle((0,0), rmax, color='black', alpha=0.2)
    ax2.add_patch(circle)
    art3d.pathpatch_2d_to_3d(circle, z=V0, zdir="z")
    ax2.plot_surface(X, Y, Z, cmap='coolwarm',norm=norm1)  # 更改 colormap
    ax2.view_init(elev=90, azim=0)  # 设置俯视视角
    ax2.set_zlim(-1, 1)
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_xticks([])  # 取消 x 轴刻度
    ax2.set_yticks([])  # 取消 y 轴刻度
    ax2.set_zticks([])  # 取消 z 轴刻度
    ax2.set_title("俯视")

    plt.tight_layout()  # 调整子图布局
    plt.title('极坐标系下波函数幅值分布')
    return fig

def cylind_color(ratio,a,z,x,y,V0,m):
    r = sqrt(x**2+y**2)
    if x > 0:
        fai = arctan(y/x)
    elif x<0:
        fai = arctan(y/x)
    else:
        fai = pi/2
    psai1, psai2 = cylindrical(ratio,a,z,r,fai,V0,m)
    Normpsai = (abs(psai1+psai2))**2
    return Normpsai

def cylind_plot(V0,ratio,a,h,m):
    #V0, ratio, a, h, m = 1, 0.2, 0.2, 1, 1
    V0, ratio, a, h, m = float(V0), float(ratio), float(a), float(h), float(m)

    rmax = a
    zmax = h
    XX = np.linspace(-rmax-dc,rmax+dc,20)
    YY = np.linspace(-rmax-dc,rmax+dc,20)
    ZZ = np.linspace(-zmax-dc,zmax+dc,20)
    x,y,z = list(XX),list(YY),list(ZZ)

    fig = plt.figure(figsize=(12, 8))  # 设置图形大小
    ax = fig.add_subplot(111, projection='3d')
    plt.rcParams['figure.figsize'] = (8.0,6.0)
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    # 绘制散点图
    values = []
    absvalues = []
    for i in x:
        for j in y:
            for k in z:
                values.append(cylind_color(ratio, a, k, i, j, V0, m))
                absvalues.append(abs(cylind_color(ratio,a,k,i,j,V0,m)))
    vmid = (min(values)+max(values))*(ratio/7.5)
    # 使用归一化对象对颜色映射进行调整
    norm1 = Normalize(vmin=min(values), vmax=vmid)  # 设置第一个颜色条的范围
    norm2 = Normalize(vmin=vmid, vmax=max(values))  # 设置第二个颜色条的范围
    for i in x:
        for j in y:
            for k in z:
                # 使用绝对值计算透明度
                alpha = 0.5+0.5*abs(cylind_color(ratio, a, k, i, j, V0, m)) / max(absvalues)  
                if cylind_color(ratio,a,k,i,j,V0,m)<=vmid:
                    ax.scatter(i, j, k, c=cylind_color(ratio, a, k, i, j, V0, m), cmap='coolwarm', norm=norm1, alpha=alpha)
                else:
                    ax.scatter(i, j, k, c=cylind_color(ratio, a, k, i, j, V0, m), cmap='viridis', norm=norm2, alpha=alpha) 

    # 绘制圆柱
    u = np.linspace(0, 2 * np.pi, 50)
    h = np.linspace(0, zmax, 20)
    xp = np.outer(rmax * np.cos(u), np.ones(len(h)))
    yp = np.outer(rmax * np.sin(u), np.ones(len(h)))
    zp = np.outer(np.ones(len(u)), h)
    ax.plot_surface(xp, yp, zp, color='black', alpha=0.2)

    # 给圆柱加个‘盖子’
    circle = Circle((0, 0), rmax, color='black', alpha=0.2)
    ax.add_patch(circle)
    art3d.pathpatch_2d_to_3d(circle, z=zmax, zdir="z")

    # 设置坐标轴标签
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # 添加颜色条
    m1 = plt.cm.ScalarMappable(norm=norm1, cmap='coolwarm')
    cb1 = plt.colorbar(m1, ax=ax, shrink=0.5, pad=0.1, location="left")  # 添加第一个颜色条
    cb1.set_ticks([min(values), vmid])  # 设置刻度
    cb1.set_ticklabels(['0', str(vmid/(max(values)+min(values)))])  # 设置刻度标签
    m2 = plt.cm.ScalarMappable(norm=norm2, cmap='viridis')
    cb2 = plt.colorbar(m2, ax=ax, shrink=0.5, pad=0.1, location="right")  # 添加第二个颜色条
    cb2.set_ticks([vmid, max(values)])  # 设置刻度
    cb2.set_ticklabels([str(vmid/(max(values)+min(values))), '1'])  # 设置刻度标签
    plt.title('柱坐标系下波函数模方分布')

    return fig

def spher_color(ratio,a,z,x,y,V0,l,m):
    r = sqrt(x**2+y**2+z**2)
    if x> 0:
        fai = arctan(y/x)
    elif x<0:
        fai = arctan(y/x)
    else:
        fai = pi/2
    if x==0 and y==0:
        theta = pi/2
    else:
        theta = arctan(z/sqrt(x**2+y**2))
    psai = spherical(ratio,a,r,theta,fai,V0,l,m)
    normpsai = (abs(psai))**2
    return normpsai

def spher_plot(V0,ratio,a,m,l):
    #V0, ratio, a, m, l =1,0.2,0.2,1,1
    V0, ratio, a, m, l =float(V0), float(ratio), float(a), float(m), float(l)
    rmax = a
    XX = np.linspace(-rmax-dc,rmax+dc,15)
    YY = np.linspace(-rmax-dc,rmax+dc,15)
    ZZ = np.linspace(-rmax-dc,rmax+dc,15)
    x,y,z = list(XX),list(YY),list(ZZ)

    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(111,projection='3d')
    plt.rcParams['figure.figsize'] = (8.0,6.0)
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    values = []
    absvalues = []
    for i in x:
        for j in y:
            for k in z:
                values.append(spher_color(ratio,a,k,i,j,V0,l,m))
                absvalues.append(abs(spher_color(ratio,a,k,i,j,V0,l,m)))

    vmid = (min(values)+max(values))*(ratio/7.5)
    # 创建第一个颜色条
    norm1 = Normalize(vmin=min(values), vmax=vmid)  # 设置第一个颜色条的范围


    # 创建第二个颜色条
    norm2 = Normalize(vmin=vmid, vmax=max(values))  # 设置第二个颜色条的范围


    for i in x:
        for j in y:
            for k in z:
                alpha = 0.5 +0.5*abs(spher_color(ratio,a,k,i,j,V0,l,m))/max(absvalues)
                if spher_color(ratio,a,k,i,j,V0,l,m)<vmid:
                    ax.scatter(i,j,k,c=spher_color(ratio,a,k,i,j,V0,l,m),cmap='coolwarm',norm=norm1,alpha=alpha)
                else:
                    ax.scatter(i,j,k,c=spher_color(ratio,a,k,i,j,V0,l,m),cmap='viridis',norm=norm2,alpha=alpha)
    #绘制球体
    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, np.pi, 50)
    x = rmax * np.outer(np.cos(u), np.sin(v))
    y = rmax * np.outer(np.sin(u), np.sin(v))
    z = rmax * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, color='black', alpha=0.1)

    #设置坐标轴标签
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    m1 = plt.cm.ScalarMappable(norm=norm1, cmap='coolwarm')  # 使用 coolwarm colormap
    cb1 = plt.colorbar(m1, ax=ax, shrink=0.5, pad=0.1, location="left")  # 添加第一个颜色条
    cb1.set_ticks([min(values),vmid])
    cb1.set_ticklabels(['0',str(vmid/(max(values)+min(values)))])

    m2 = plt.cm.ScalarMappable(norm=norm2, cmap='viridis')  # 使用 viridis colormap
    cb2 = plt.colorbar(m2, ax=ax, shrink=0.5, pad=0.1, location="right")  # 添加第二个颜色条
    cb2.set_ticks([vmid,max(values)])
    cb2.set_ticklabels([str(vmid/(max(values)+min(values))),'1'])
    plt.title('球坐标系下波函数模方分布')
    return fig

def multicar_plot(V,ratio,a):
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    plt.rcParams['figure.figsize'] = (8.0,6.0)
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False

    xmin, xmax = 0, dc ##############  方便计算
    ncar = len(V)
    x1 = np.linspace(-dc,xmin,1000)
    x2 = np.linspace(xmin,xmax,1000)
    x3 = np.linspace(xmax+dc*(ncar-1),xmax+dc*ncar,1000)
    x = [x1,x2]
    for i in range(ncar):
        xi = np.linspace(xmax+dc*(i),xmax+dc*(i+1),1000)
        x.append(xi)

    

#animation
def time_factor(V0,ratio,t):
    j=1j
    hbar=16/k0**2
    E = exp(j*ratio*V0*t)
    return E

def ani_car_plot(V0,ratio,a):
    #V0,ratio,a = 1,0.2 ,0.2
    V0, ratio, a =float(V0), float(ratio), float(a)
    xmin, xmax = 0, a
    x = np.linspace(-xmin-dc,xmax+dc, 300)
    x1 = np.linspace(-dc,xmin,1000)
    x2 = np.linspace(xmin,xmax,1000)
    x3 = np.linspace(xmax,xmax+dc,1000)



    y11, y12 = cartesian(ratio,a,x1,V0)
    y21, y22 = cartesian(ratio,a,x2,V0)
    y31, y32 = cartesian(ratio,a,x3,V0)
    y1, y2 = [], []
    for i in range(len(y11)):
        y1.append(y11[i]+y12[i])
    for i in range(len(y21)):
        y2.append(y21[i]+y22[i])
    yV0 = np.zeros(len(x))
    yV0[(x>xmin) & (x<xmax)] = V0

    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)

    plt.rcParams['figure.figsize'] = (8.0,6.0)
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    line11,=ax1.plot(x1,y11,label='入射波')
    line12,=ax1.plot(x1,y12,label='反射波')
    line1, =ax1.plot(x1,y1,label='复合波1')
    line2 ,=ax1.plot(x2,y2,label='复合波2')
    line31,=ax1.plot(x3,y31,label='透射波')
    def init():
        line11.set_ydata([np.nan]*len(x1))
        line12.set_ydata([np.nan]*len(x1))
        line1.set_ydata([np.nan]*len(x1))
        line2.set_ydata([np.nan]*len(x2))
        line31.set_ydata([np.nan]*len(x3))
        return line1, ; line2, ; line11, ; line12, ; line31,

    def update(frame):
        for i in range(len(y11)):
            y11[i] = (y11[i]*time_factor(V0,ratio,float(frame)/500))
            y12[i] = (y12[i]*time_factor(V0,ratio,float(frame)/500))
            y1[i] = (y1[i]*time_factor(V0,ratio,float(frame)/500))
        for i in range(len(y2)):
            y2[i] = (y2[i]*time_factor(V0,ratio,float(frame)/500))
        for i in range(len(y31)):
            y31[i] = (y31[i]*time_factor(V0,ratio,float(frame)/500))
        line11.set_ydata(y11)
        line12.set_ydata(y12)
        line1.set_ydata(y1)
        line2.set_ydata(y2)
        line31.set_ydata(y31)
        return line1,  line2,  line11,  line12,  line31,

    ax1.plot(x,yV0,label='势')
    ax1.set_xticks(np.arange(-xmin-dc,xmax+dc,0.2))
    ax1.legend()
    plt.title('一维情况下波函数幅值分布 E/V0={}'.format(ratio))
    ani = FuncAnimation(fig,update,init_func=init,frames=2000,interval=50,blit=True)
    return ani

