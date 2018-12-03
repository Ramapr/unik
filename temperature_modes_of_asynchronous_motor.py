"""
Данный скрипт позволит рассчитать переходные 
процесссы нагревания асинхронного двигателя 
при разных режимах работы (S1, S2, S3).
"""
%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt

'зададим все значения'
s1, s2, s3 = 180*60, 90*60, 0.6
t_dop = 110
t = 40
n = 750
etha = 0.74
P_n = 1100
m = 22.3
C = 381 
tmax = 70
ttime = np.linspace(0, s1/60, s1, endpoint=True)
time = np.linspace(0, s1, s1, endpoint=True)
T = 5 # число циклов для режима S3
alpha = 0
if n < 1000:
    alpha = 0.5
else:
    alpha = 0.7
xtick = s1 // 60 + 10
ytick = t_dop + t + 10

def ps1(Pn, t_oc, t_dop, alpha=0.5):
    """
	Формула для вычисления полезной механической 
	мощности двигателя при режиме S1 (темп. охл. среды отличной от 40С).
	Pn - номинаяльная мощность
	t_oc - температура фактич. охлаждающей среды
	t_dop - допускаемое превышение температуры
	alpha - коэф. потерь
    """
	
    dt = 40 - t_oc 
    return Pn*(1 + dt*(1 + alpha) / t_dop)**0.5
    
def ps2(Pn, tp, Tn, alpha=0.5):
	"""
	Мех. мощность двигателя при режиме S2.
	Pn - номинаяльная мощность
	tp - продолжительность работы двигателя [c]
	Tn - постоянная времени нагревания
	alpha - коэф. потерь
	"""
    
    return Pn*(((1+alpha)/(1 - np.exp(-1*tp/Tn))) - alpha) **0.5
    
def ps3(Pn, csi, alpha=0.5):
    """
	Мех. мощность двигателя при режиме S3.
	Pn - номинаяльная мощность
	csi - относительная продолжительность включения
	alpha - коэф. потерь
	"""
	
    return Pn / (csi/(csi + (1 + alpha)*(1 - csi)))**0.5
	
def warm_tnom(t, t_ust, t_start, Tn):
	"""
	Уравнение ПП при нагревании, когда температура двигателя 
	не равна температуре окружающей среды.
	t - время 
	t_ust - установившееся значение температуры
	t_start - начальное превышение температуры (в 0 момень времени)
	Tn - постоянная времени нагревания
	"""
	
    return t_ust*(1 - np.exp(-1 * t / Tn)) + t_start * np.exp(-1 * t / Tn)

def warm_t(t, t_ust, Tn):
    """
	Уравнение ПП при нагревании, когда температура
	двигателя равна температуре окружающей среды.
	t - время 
	t_ust - установившееся значение температуры
	Tn - постоянная времени нагревания
	"""
	
    return t_ust * (1 - np.exp(-1 * t / Tn))

def cold(t, t_ust, T_cold):
	"""
	Уравнение ПП при охлаждении.
	t - время 
	t_ust - установившееся значение температуры
	T_cold - постоянная времени охлаждения
	"""
	
    return t_ust * np.exp(-1 * t / T_cold)
	
	
'рассчитаем режим S1'
Ps1_40 = ps1(P_n, 40, t_dop, alpha) # (13)
Ps1_ = ps1(P_n, 30, t_dop)

delta_Pt_40 = Ps1_40 * (1 - etha) / etha # (7)
delta_Pt = Ps1_ * (1 - etha) / etha

A = P_n * (1 - etha) / etha / t_dop # (6)
T_n = C * m / A # (5)

t_ust_40 = delta_Pt_40 / A
t_ust = delta_Pt / A

y_40 = np.array(warm_t(time, t_ust_40, T_n))
y_ = np.array(warm_t(time, t_ust, T_n))

plt.figure(0)
plt.grid(True)
plt.plot(ttime, y_40 + t, 'b', label='t=40°C')
plt.plot(ttime, y_ + 30, 'r', label='t=30°C')
plt.legend(loc=4, borderaxespad=1)
plt.title('Режим S1')
plt.xticks(np.arange(0, xtick, 20))
plt.yticks(np.arange(0, ytick, 10))
plt.xlabel('t, min')
plt.ylabel(r'$\tau$, °C')
#plt.savefig('fig1', fmt='png')


'рассчитаем режим S2'
Ps2 = ps2(P_n, s2, T_n, alpha)
delta_Pt_2 = Ps2 * (1 - etha) / etha
t_ust_s2 = delta_Pt_2 / A

y_s2w = np.array(warm_t(time[:s2], t_ust_s2, T_n))
y_s2c = np.array(cold(time[:s2], y_s2w[len(y_s2w) - 1], T_n))
y_s2 = np.hstack((y_s2w, y_s2c))

plt.figure(1)
plt.grid(True)
plt.plot(ttime, y_s2 + t, 'b')
plt.xticks(np.arange(0, xtick, 20))
plt.yticks(np.arange(0, ytick, 10))
plt.title('Режим S2')
plt.xlabel('t, min')
plt.ylabel(r'$\tau$, °C')
#plt.savefig('fig2', fmt='png')


'рассчитаем режим S3'
Ps3 = ps3(P_n, s3, alpha)
delta_Pt_3 = Ps3 * (1 - etha) / etha
t_ust_s3 = delta_Pt_3 / A

tt = len(time) // T
warm_time = time[:int(tt * s3)]
cold_time = time[:int(tt * (1 - s3))] 
y_s3 = np.array(warm_t(warm_time, t_ust_s3, T_n))

for i in range(T):
    sub_1 = np.array(cold(cold_time, y_s3[len(y_s3) - 1], T_n))
    y_s3 = np.hstack((y_s3, sub_1))
    if i == 4:
        break
    sub_2 = np.array(warm_tnom(warm_time, t_ust_s3, y_s3[len(y_s3) - 1], T_n))
    y_s3 = np.hstack((y_s3, sub_2))

plt.figure(2)
plt.grid(True)
plt.plot(ttime, y_s3 + t, 'b')
plt.xticks(np.arange(0, 190, 20))
plt.yticks(np.arange(0, 160, 10))
plt.xlabel('t, min')
plt.ylabel(r'$\tau$, °C')
plt.title('Режим S3')
#plt.savefig('fig3', fmt='png')
