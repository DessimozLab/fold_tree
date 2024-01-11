import subprocess
cmd = snakemake.params.cmd
print('running' , cmd)
p = subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
stdout, stderr= p.communicate()
if p.returncode == 0:
    print('OK')
else:
    # Analyze exit code and stderr and decide what to do next
    print(p.returncode)
    print(stderr.decode())
    #write empty file
    open( snakemake.output[0],  'a').close()
    print('empty output file', snakemake.output[0])

