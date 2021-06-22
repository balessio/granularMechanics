#filenames = ['../hertz/log_vStep_1e0','../hertz/log_vStep_1en2','../hertz/log_vStep_1en3',
#'../hooke/log_vStep_1e0','../hooke/log_vStep_1en2','../hooke/log_vStep_1en3']

#filenames = ['hooke/log_vStep_1e0','hooke/log_vStep_1en2','hooke/log_vStep_1en3']
#filenames = ['hertz/continue/log_continue_1en2','hertz/continue/log_continue_1en3']
filenames = ['quasi 2D hertz/gouge_shear_3D/log_vStep_1e0',
'quasi 2D hertz/gouge_shear_3D/log_vStep_1en2',
'quasi 2D hertz/gouge_shear_3D/log_vStep_1en3']


for fileName in filenames:
	with open(fileName+'.txt', 'r') as log:
		lines = log.readlines()
	with open(fileName+'_noText.txt', 'w') as fix:
		for line in lines:
			try:
				vals = line.split()
				a = float(vals[13])
				if len(vals)==14:
					fix.write(line)
			except:	continue
