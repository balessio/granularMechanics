import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib
from scipy import optimize

# CONSTANTS
################################
################################
gouge_thickness = 6e-2 # m
gouge_length = 16e-2 # m
gouge_depth = 1e-3 # m
A = gouge_length*gouge_depth # m^2
d_c = 0.13*gouge_thickness
sig_n = 50e6 # Pa
a_minus_b = -6.44e-3 / np.log(10)
################################
################################
def calculations(k_spring, slope):
	G = slope*sig_n*gouge_thickness
	k_eff = slope*sig_n*A 
	k_H = k_spring*k_eff/(k_spring-k_eff)
	k_eff_critical = -sig_n*a_minus_b*A/d_c
	k_spring_critical = k_H*k_eff_critical/(k_H-k_eff_critical)
	printstr = 'calculations for $k_{sp}$=1e7N/m measurement:'
	printstr += '\nshear modulus G = '+str(G)+' Pa'
	printstr += '\ngouge stiffnes k_H = '+str(k_H)+' N/m'
	printstr+='\ncritical pulling spring stiffness k_sp_crit = '+str(k_spring_critical)+' N/m'
	print(printstr)

v = 10 # cm/s

k_sp = 7e5 # N/m
cutoffs = [[0.05,0.4]]
files = ['res99k7e5.txt']
dts = [9.33e-9]
labels = ['$k_{sp}$=7e5N/m']

def line_func(LPD, slope_1, offset):
	return slope_1*LPD + offset

lpd = [] # load point displacement = v*t
friction = []
sub_fric = []
sub_lpd = []
for file, dt, cutoff in zip(files,dts,cutoffs):
	sim_data = np.loadtxt(file)
	sim_length = len(sim_data)
	lpd_here = np.array((sim_data[0:sim_length,0]-sim_data[0,0]))*v*dt
	lpd.append(lpd_here)
	friction_here = sim_data[0:sim_length,2]/sim_data[0:sim_length,3]
	friction.append(friction_here)
	'''
	lpd_diff = []
	for t_a, t_b in zip(lpd_here[:200], lpd_here[1:201]):
		lpd_diff.append(t_b-t_a)
	step = np.mean(lpd_diff)

	f_s = 1/step
	lpd_i, lpd_f = int(cutoff[0]*f_s), int(cutoff[1]*f_s)
	'''
	lpd_i, lpd_f = 0, 0
	for i in range(lpd_here.size):
		if (lpd_here[i] < cutoff[0]):
			continue
		elif (lpd_i == 0):
			lpd_i = i 
		elif (lpd_here[i] > cutoff[1]):
			lpd_f = i 
			break

	sub_friction = friction_here[lpd_i:lpd_f]
	sub_LPD = lpd_here[lpd_i:lpd_f]
	
	slope_guess = (np.mean(sub_friction[-3:])-np.mean(sub_friction[:3]))/(
		lpd_here[lpd_f]-lpd_here[lpd_i])

	guess_not = [slope_guess, 0.3]
	params, params_cov = optimize.curve_fit(
		line_func, sub_LPD, sub_friction,
		p0=guess_not, sigma=np.array([1e-5]*len(sub_LPD)))
	slope_err = params_cov[0][0]**0.5
	printstring = labels[0]+': '+'slope: '+str(np.round(params[0],5))
	printstring += '+/-'+str(np.round(slope_err,5))
	printstring += ', '+'y-offset: '+str(np.round(params[1],5))
	print(printstring)
	sub_fric.append(sub_friction)
	sub_lpd.append(sub_LPD)

plt.rc('font', family='serif', serif = "cmr10", size=20)
plt.rcParams['mathtext.fontset']='cm'
plt.figure(figsize=(10,6))
matplotlib.rcParams['axes.unicode_minus'] = False

for x, y in zip(sub_lpd, sub_fric):
	plt.errorbar(x=x, y=y, linewidth=1, fmt='', label='data')
	plt.errorbar(x=x, y=line_func(x,params[0],params[1]),
		linewidth=1, fmt='', label='linear fit')

#plt.legend(loc='best')
plt.xlabel('load point displacement (cm)')
plt.ylabel('initial loading friction $\mu_{ss}$')
plt.title('Initial loading friction with linear fit\n$P$=50MPa, $v_{lp}$=10cm/s\n'+printstring)
plt.tight_layout(pad=0)
plt.grid(linestyle='--')
plt.show()

calculations(k_sp,params[0])