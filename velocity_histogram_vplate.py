import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib

# calculates how many lines to skip when calculating velocity
def dn(dLPD, dT, v_SP, dOUT, ALPHA=10):
	denominator = dT*v_SP*dOUT*alpha
	return int(dLPD/denominator)

# output format in the log file
# 0.Step 1.Atoms 2.f_brlc[1] 3.f_brlc[2] 4.f_brlc[3] 5.f_trlc[1]
# 6.f_trlc[2] 7.f_trlc[3] 8.c_3[1] 9.c_3[2] 10.c_3[3] 11.c_4[1]
# 12.c_4[2] 13.c_4[3] 14.f_pull[1] 15.f_pull[2] 16.f_pull[3]
# 17.f_pull[4] 18.f_pull[5] 19.f_pull[6] 20.f_pull[7] 

#c_3[1], c_3[2], and c_3[3] are the x-, y- and z-dir position of the center...
# ...of mass of the bottom (shearing) block, respectively,
#c_4[1], c_4[2], and c_4[3] are the x-, y- and z-dir position of the center...
# ...of mass of the top (fixed) block, respectively,
#f_pull[1] is the x-component of the spring force (the spring is 1-D in x-dir, so...
# ...it only has x-components of force and position), and f_pull[5] is the x-dir...
# ....load point position (the end of the spring not attached to the block).

print('started...')

k = '1e8'
delta_lpd = 2.5e-3 # cm
v_sp = 10 # cm/s
dt = 2.3325e-9 # s
delta_out = 400
alpha = 25
derivative_steps = dn(delta_lpd,dt,v_sp,delta_out,alpha)
num_bins = 150

output_textfile = 'velocity_histograms/k'+k+'_veloctity_histogram.txt'
file = 'k'+k+'.txt'

try: sim_data = np.loadtxt(file)
except: sim_data = np.loadtxt('res99'+file)
sim_length = len(sim_data)

# number of time steps between outputs
num_dt = sim_data[1,0] - sim_data[0,0] 
# x position array of the center of mass of the sheared (bottom) block 
sheared_block_x = np.array(sim_data[0:sim_length,8])

# caluculate the instantaneous derivative as a function of time
############
time_window = derivative_steps*dt
velocities = []
for i, j in zip(range(0, sheared_block_x.size-derivative_steps, derivative_steps),
	range(derivative_steps, sheared_block_x.size, derivative_steps)):
	instantaneous_velocity = (sheared_block_x[j]-sheared_block_x[i])/time_window
	velocities.append(instantaneous_velocity)

# time array with centered derivative increments
time = np.arange(
	0, sheared_block_x.size-derivative_steps, derivative_steps) + derivative_steps/2
# convert time to mega seconds
time /= 1e6
print('max velocity: '+str(np.max(velocities)))

plt.rc('font', family='serif', serif = "cmr10", size=20)
plt.rcParams['mathtext.fontset']='cm'
plt.figure(figsize=(10,6))
matplotlib.rcParams['axes.unicode_minus'] = False

plt.errorbar(x=time, y=velocities, linewidth=1, fmt='')

#plt.legend(loc='lower right', fontsize=12)
plt.xlabel('time $t$ (Ms)')
plt.ylabel('\ninstantaneous x-velocity $v_x$ (m/s)')
plt.title('\nnumber of log file lines between $v$ calculations = '+str(
	derivative_steps)+'\n$v_{sp}$=10cm/s, $k_{sp}=$'+k+'N/m\n')
plt.tight_layout(pad=0)
plt.grid(linestyle='--')
#plt.ylim(bottom=0)
#plt.xscale('log')
plt.show()

hist, bins = np.histogram(velocities, bins=np.logspace(-6, 7, num_bins))
hist = np.array(hist).astype(np.float32)
centers = []
for edge1, edge2 in zip(bins[:-1], bins[1:]):
	centers += [edge1 + ((edge2-edge1)/2)]

with open(output_textfile, 'w') as ttt:
	for x,y in zip(centers,hist):
		ttt.write(str(x)+' '+str(y)+'\n')

plt.rc('font', family='serif', serif = "cmr10", size=20)
plt.rcParams['mathtext.fontset']='cm'
plt.figure(figsize=(10,6))
matplotlib.rcParams['axes.unicode_minus'] = False

plt.stem(centers, hist, 'r', use_line_collection=True)

plt.xscale('log')
plt.yscale('log')
#plt.xlim(left=1e-4)
#plt.legend(loc='best')
plt.xlabel('instantaneous velocity (m/s)\n')
plt.ylabel('\noccurences')
plt.title('Histogram of instantaneous velocities, bins='+str(
	num_bins)+'\nnumber of log file lines between $v$ calculations = '+str(
	derivative_steps)+'\n$v_{sp}$=10cm/s, $k_{sp}=$'+k+'N/m\n')
plt.tight_layout(pad=0)
plt.grid(linestyle='--')
plt.show()
