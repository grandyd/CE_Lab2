import matplotlib.pyplot as plt
import numpy as np

def ACH_func(v,h,F):
   A=F*((k*(K-M*v**2)-(m*v**2)*(K+k-M*v**2))**2+((h*v)**2*(K-(M+m)*v**2)**2))**(-1/2)*((k-m*v**2)**2+(h*v)**2)**(1/2)
   ACH = A / F
   return ACH


#* Функция рисовки АЧХ
def plot_ACH(h,F):
   v=np.linspace(0,3)
   ACH_H=ACH_func(v,h,F)
   ACH_0=ACH_func(v,0,F)


   plt.plot(t, ACH_H, 'o-', linewidth = 2.0, label ='АЧХ с трением h = '+str(h))
   plt.plot(t, ACH_0, 'x-', linewidth = 2.0, label ='АЧХ без трения h = 0')
   plt.xlabel('v')
   plt.ylabel('A(v)')
   plt.legend()
   plt.grid()
   plt.savefig('ACH.png', bbox_inches='tight')
   plt.show()

def amplitude(t,h):
   v = (k/m)**(1/2)
   A = (F*(h*1j*v + k - m*(v**2)))/((K + h*1j*v + k - M*(v**2))*(h*1j*v + k - m*(v**2)) - ((h*1j*v + k)**2))
   A = abs(A)
   x = A*np.exp(1j*v*t)
   B = (F*(h*1j*v + k))/((K + h*1j*v + k - M*(v**2))*(h*1j*v + k - m*(v**2)) - ((h*1j*v + k)**2))
   y = B*np.exp(1j*v*t)
   return [t, x, y]

def plot_amplitude(t,h):
   t, x, y = amplitude(t,h)
   plt.plot(t, x, 'o-', linewidth = 2.0, label ='x(t) - Колебания основного тела')
   plt.plot(t, y, 'x-', linewidth = 2.0, label ='y(t) - Колебания демпфера')
   plt.xlabel('t')
   plt.ylabel('x(t)     y(t)')
   plt.legend()
   plt.grid()
   if h==0:
      plt.savefig('Amplitude_without_h.png', bbox_inches='tight')
   else:
      plt.savefig('Amplitude_with_h.png', bbox_inches='tight')
   plt.show()

# демпфер
def damper_func(t,h,m,k):
   v = (k/m)**(1/2)
   B = (F*(h*1j*v + k))/((K + h*1j*v + k - M*(v**2))*(h*1j*v + k - m*(v**2)) - ((h*1j*v + k)**2))
   y = abs(B)*np.exp(1j*v*t)
   return y
def plot_damper(t,h,m,k,M,K,name):
   y = damper_func(t,h,m,k)
   plt.plot(t, y, 'o-', linewidth = 2.0)
   plt.xlabel('t')
   plt.ylabel('y(t)')
   plt.grid()
   plt.savefig(name+'.png', bbox_inches='tight')
   plt.show()

# основное тело
def main_func(t,h,m,k):
   v = (k/m)**(1/2)
   A = (F*(h*1j*v + k - m*(v**2)))/((K + h*1j*v + k - M*(v**2))*(h*1j*v + k - m*(v**2)) - ((h*1j*v + k)**2))
   x = abs(A)*np.exp(1j*v*t)
   return x
def plot_main(t,h,m,k,M,K,name):
   x = main_func(t,h,m,k)
   plt.plot(t, x, 'o-', linewidth = 2.0)
   plt.xlabel('t')
   plt.ylabel('x(t)')
   plt.grid()
   plt.savefig(name+'.png', bbox_inches='tight')
   plt.show()

if __name__=='__main__':
   K = 5
   k = 2
   M = 13
   m = 10
   F = 7
   h =1
   t=np.linspace(0,40)

   plot_ACH(h,F)

   plot_amplitude(t,h)
   plot_amplitude(t,0)

   t=np.linspace(0,40)

   name='MmoreK'
   plot_main(t,h,m,k,M,K,name)
   name='MlessK'
   K = 15
   plot_main(t, h, m, k, M, K, name)
   name = 'MequalK'
   K = M
   plot_main(t, h, m, k, M, K, name)

   t = np.linspace(0, 20)
   M=13
   K=5
   name='1mMOREk'
   plot_damper(t,h,m,k,M,K,name)
   name='1mLESSk'
   k=11
   plot_damper(t,h,m,k,M,K,name)
   name = '1mEQUALk'
   m=k
   plot_damper(t, h, m, k, M, K, name)

