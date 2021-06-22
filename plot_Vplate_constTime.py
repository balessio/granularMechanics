import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib

print('started...')

v = 10 # cm/s
delta_out = 1000

files = ['gouge_shear_2D_hertz_Lx5_50MPa_log.txt',
	'gouge_shear_2D_hooke_Lx5_50MPa_log.txt',
	'gouge_shear_3D_hertz_thinLy_Lx5_50MPa_log.txt']
dts = [8e-09,8e-09,2e-08]
labels = ['2D hertz','2D hooke','3D hertz (thin $\ell_y$)']
derivative_steps = 10

'''
files = ['gouge_shear_2D_hertz_Lx5_50MPa_log.txt',
	'gouge_shear_2D_hooke_Lx5_50MPa_log.txt']
dts = [8e-09,8e-09]
labels = ['2D hertz','2D hooke']
'''

times = []
velocities = []
for file, dt in zip(files,dts):
	sim_data = np.loadtxt(file)
	sim_length = len(sim_data)
	# x position array of the center of mass of the sheared (bottom) block 
	sheared_block_x = np.array(sim_data[0:sim_length,8])
	# caluculate the instantaneous derivative as a function of time
	############
	time_window = derivative_steps*dt*delta_out
	velocity = []
	for i, j in zip(range(0, sheared_block_x.size-derivative_steps, derivative_steps),
		range(derivative_steps, sheared_block_x.size, derivative_steps)):
		instantaneous_velocity = (sheared_block_x[j]-sheared_block_x[i])/time_window
		velocity.append(instantaneous_velocity)
	print(min(np.abs(velocity)))
	# time array with centered derivative increments
	time = np.arange(
		0, sheared_block_x.size-derivative_steps, derivative_steps) + derivative_steps/2
	# convert time to mega seconds
	print('max velocity: '+str(np.max(velocity)))
	times.append(time/1e6)
	velocities.append(np.array(velocity)*100)

plt.rc('font', family='serif', serif = "cmr10", size=20)
plt.rcParams['mathtext.fontset']='cm'
plt.figure(figsize=(10,6))
matplotlib.rcParams['axes.unicode_minus'] = False

for x, y, label in zip(times, velocities, labels):
	avg_string = '$<V_{plate}>$='+str(round(np.mean(y),4))+' cm/s'
	std_string = '$\sigma$='+str(round(np.std(y),4))+' cm/s'
	plt.errorbar(x=x, y=y, linewidth=1, fmt='', label=label+': '+avg_string+', '+std_string)
#plt.errorbar(x=(sim_data[lines_until_step-1,0]-sim_data[0,0]
#	)*v*dt, y=0.29, yerr=0.05,linewidth=3)

plt.legend(loc='lower right', fontsize=20)
plt.xlabel('time (s)\n')
plt.ylabel('\nshearing plate velocity (cm/s)')
plt.title('$V_{slide}=10$ cm/s, varying contact law\n')
plt.tight_layout(pad=0)
plt.grid(linestyle='--')
plt.ylim(bottom=0)
#plt.xscale('log')
plt.show()