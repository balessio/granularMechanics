import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib
from scipy import integrate
from scipy import signal as sig 

# output format in the log file
# 0.Step 1.Atoms 2.f_brlc[1] 3.f_brlc[2] 4.f_brlc[3] 5.f_trlc[1]
# 6.f_trlc[2] 7.f_trlc[3] 8.c_3[1] 9.c_3[2] 10.c_3[3] 11.c_4[1]
# 12.c_4[2] 13.c_4[3]

#c_3[1], c_3[2], and c_3[3] are the x-, y- and z-dir position of the center...
# ...of mass of the bottom (shearing) block, respectively,
#c_4[1], c_4[2], and c_4[3] are the x-, y- and z-dir position of the center...
# ...of mass of the top (fixed) block, respectively,

print('started...')

v_sp = 10 # cm/s
delta_out = 1000
num_bins = 150
sigma_n = 1e6 # Pa
slip_threshold = 1e-5 # m
interval = (0.0,0.0)

files = ['gouge_shear_2D_hertz_Lx5_50MPa_log.txt',
	'gouge_shear_2D_hooke_Lx5_50MPa_log.txt',
	'gouge_shear_3D_hertz_thinLy_Lx5_50MPa_log.txt']
dts = [8e-09,8e-09,2e-08]
labels = ['2D hertz','2D hooke','3D hertz (thin $\ell_y$)']

which_file = 0
file, dt = files[which_file], dts[which_file]

sim_data = np.loadtxt(file)
sim_length = len(sim_data)

# number of time steps between outputs
num_dt = sim_data[1,0] - sim_data[0,0] 
# x position array of the center of mass of the sheared (bottom) block 
sheared_block_x = np.array(sim_data[0:sim_length,8])
time = np.arange(0,sim_length)*dt*delta_out

next_window, d_slip, index = True, 0, 0
velocities, time_ragged, slips = [], [], []
for i in range(0,len(sheared_block_x)):
	if next_window:
		next_window = False
		index = i
	d_slip = sheared_block_x[i]-sheared_block_x[index]
	if d_slip > slip_threshold:
		time_ragged.append((time[i]+time[index])/2)
		slips.append(d_slip*100)
		velocities.append(d_slip/(time[i]-time[index]))
		next_window = True

vFitLbound, vFitRbound = 0, 0
for i in range(0,len(time_ragged)):
	if time_ragged[i] >= interval[0]:
		vFitLbound = i 
		break
if not interval[1] == 0.0:
	for i in range(vFitLbound,len(time_ragged)):
		if time_ragged[i] >= interval[1]:
			vFitRbound = i 
			break

if not vFitRbound==0:
	subV = velocities[vFitLbound:vFitRbound]
else:
	subV = velocities[vFitLbound:]

hist, bins = np.histogram(subV, bins=np.logspace(-8, 5, num_bins))
hist = np.array(hist).astype(np.float32)
centers = []
for edge1, edge2 in zip(bins[:-1], bins[1:]):
	centers += [edge1 + ((edge2-edge1)/2)]

slip_hist, avg_tracker = hist*0, hist*0
for i in range(0,len(velocities)):
	index = (np.abs(centers - velocities[i])).argmin()
	slip_hist[index] += slips[i]
	avg_tracker[index] += 1
for i in range(0,len(slip_hist)):
	if avg_tracker[i] > 0:
		slip_hist[i] /= avg_tracker[i]

D_fromIntegratedHist = np.sum(slip_hist*hist) # need to convert from v to d
print('total displacement='+str(D_fromIntegratedHist)+' m')
'''
plt.rc('font', family='serif', serif = "cmr10", size=20)
plt.rcParams['mathtext.fontset']='cm'
plt.figure(figsize=(10,6))
matplotlib.rcParams['axes.unicode_minus'] = False

plt.stem(centers, hist, 'r', use_line_collection=True)

plt.xscale('log')
plt.yscale('log')
#plt.xlim(left=1e-4)
#plt.legend(loc='best')
plt.xlabel('plate velocity (m/s)\n')
plt.ylabel('\noccurences')
plt.title('Histogram of plate velocities, bins='+str(
	num_bins)+'\n$v_{sp}$=10cm/s,\ntotal plate displacement='+str(
		np.round(D_fromIntegratedHist,4))+' m\n'\
	+'constanct slip distance of '+str(slip_threshold)+' cm')
plt.tight_layout(pad=0)
plt.grid(linestyle='--')
plt.show()
'''

plt.rc('font', family='serif', serif = "cmr10", size=20)
plt.rcParams['mathtext.fontset']='cm'
plt.figure(figsize=(10,6))
matplotlib.rcParams['axes.unicode_minus'] = False

plt.errorbar(x=time_ragged, y=np.array(velocities)*100, linewidth=2)

#plt.xlim(left=1e-4)
#plt.legend(loc='best')
plt.xlabel('time (s)')
plt.ylabel('Plate velocity (cm/s)')
plt.title('plate velocities for $v_{sp}$=10 cm/s\n$<v_{plate}>$='\
	+str(np.round(np.mean(np.array(velocities)*100),4))+' cm/s, '\
	+'$\sigma$='+str(np.round(np.std(np.array(velocities)*100),6))+' cm/s\n'\
	+'calculated for constant slip distance of '+str(slip_threshold)+' m')
plt.tight_layout(pad=0)
plt.grid(linestyle='--')
plt.show()

plt.rc('font', family='serif', serif = "cmr10", size=20)
plt.rcParams['mathtext.fontset']='cm'
plt.figure(figsize=(10,6))
matplotlib.rcParams['axes.unicode_minus'] = False

plt.errorbar(x=time_ragged, y=slips, linewidth=2)

#plt.xlim(left=1e-4)
#plt.legend(loc='best')
plt.xlabel('time (s)')
plt.ylabel('slip distance (cm)')
plt.title('slip distance, expected to be constant '+str(slip_threshold)+' m')
plt.tight_layout(pad=0)
plt.grid(linestyle='--')
plt.show()
