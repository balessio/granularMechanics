import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib

print('started...')

files_beforeStep = ['gouge_shear_2D_hertz_Lx5_50MPa_log.txt',
	'gouge_shear_2D_hooke_Lx5_50MPa_log.txt']
vs_beforeStep = [10,10]
dts_beforeStep = [8e-09]*len(files_beforeStep)


files_afterStep = ['log_hertz_vStep_1e0_noText.txt','log_hertz_vStep_1en2_noText.txt',
'log_hertz_vStep_1en3_noText.txt','log_hooke_vStep_1e0_noText.txt',
'log_hooke_vStep_1en2_noText.txt','log_hooke_vStep_1en3_noText.txt']
vs_afterStep = [100,1,0.1,100,1,0.1]
dts_afterStep = [8e-09]*len(files_afterStep)
line_of_step = [1060000000,1060000000,1060000000,980000000+50000000,980000000+50000000,980000000+50000000]
labels = ['hertz $V=100$cm/s','hertz $V=1$cm/s','hertz $V=0.1$cm/s',
'hooke $V=100$cm/s','hooke $V=1$cm/s','hooke $V=0.1$cm/s']

#############
which_one = 5
files_afterStep = [files_afterStep[which_one]]
vs_afterStep = [vs_afterStep[which_one]]
dts_afterStep = [8e-09]*len(files_afterStep)
line_of_step = [line_of_step[which_one]]
labels = [labels[which_one]]


lpd = [] # load point displacement = v*t
friction = [] # friction = f_x/f_z
vStep_loc = []
hertz_or_hooke = 0
for file, vStep, dt, lineStep in zip(files_afterStep,vs_afterStep,dts_afterStep,line_of_step):
	sim_data = np.loadtxt(file)
	sim_length = len(sim_data)
	if 'hertz' in file: hertz_or_hooke = 0
	if 'hooke' in file: hertz_or_hooke = 1
	previousData = np.loadtxt(files_beforeStep[hertz_or_hooke])
	previousData = previousData[0:np.where(previousData[0:,0]==lineStep)[0][0],0:]
	vPrev = vs_beforeStep[hertz_or_hooke]
	dtPrev = dts_beforeStep[hertz_or_hooke]
	lpd_here = np.array((sim_data[0:sim_length,0]-sim_data[0,0]))*vStep*dt
	lpdPrev = np.array((previousData[0:,0]-previousData[0,0]))*vPrev*dtPrev
	lpd_combined = np.array(list(lpdPrev)+list(lpd_here+lpdPrev[-1]))
	lpd.append(lpd_combined)
	friction_here = sim_data[0:sim_length,2]/sim_data[0:sim_length,3]
	friction_prev = previousData[0:,2]/previousData[0:,3]
	#if '3D' in file: friction_here = sim_data[0:sim_length,2]/sim_data[0:sim_length,4]
	friction_combined = np.array(list(friction_prev)+list(friction_here))
	friction.append(friction_combined)

plt.rc('font', family='serif', serif = "cmr10", size=20)
plt.rcParams['mathtext.fontset']='cm'
plt.figure(figsize=(10,6))
matplotlib.rcParams['axes.unicode_minus'] = False

for x, y, label in zip(lpd, friction, labels):
	plt.errorbar(x=x, y=y, linewidth=1, fmt='', label=label)

plt.legend(loc='lower right', fontsize=12)
plt.xlabel('Load-point displacement (cm)')
plt.ylabel('normalized shear stress ($f_x/f_z$)')
plt.title('velocity step from $V_{slide}=10$cm/s, varying contact law\n')
plt.tight_layout(pad=0)
plt.grid(linestyle='--')
plt.ylim(bottom=0)
#plt.xscale('log')
plt.show()
