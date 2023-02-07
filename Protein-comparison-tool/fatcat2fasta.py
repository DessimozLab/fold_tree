
def FatcatToFasta(fatout, filename):
	print fatout 
	print filename
	alns = ['','']
	files = []
	BlockNext = False
	Blocknum = 0
	with open( fatout , 'r' ) as fatfile:
		for line in fatfile:
			try:
				print line
				if 'file from local' in line:
					files.append(line.split('/')[-1])	
				#format  = Chain 1:   16 PHC---------------SKTPIVRAQTSQNAMS-------RGMQMQFSIGLHT---------AVC----
				if 'Chain 1:' in line:
					alns[0] += line.split(':')[1].strip().split(' ')[1]
				if 'Chain 2:' in line:
					alns[1] += line.split(':')[1].strip().split(' ')[1]
			except:
				print 'line error:'
				print line
		alnstr1 = '>'+files[0][:-1] + 'Fatcat_block' +'\n' +alns[0] + '\n'
		alnstr2 = '>'+files[1][:-1] + 'Fatcat_block' +'\n' +alns[1] + '\n'
		print alnstr1
		print alnstr2
		handle = open(filename  , 'w')
		handle.write(alnstr1 + '\n')
		handle.write(alnstr2)
		handle.close()
		return filename
