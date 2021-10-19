

import subprocess, datetime

doy  = [i for i in range(38,39)]
day  = [i for i in range(7,14)]
month= [2 for i in range(7)]
year = [2021 for i in range(7)]

CMD_EXE= r'Cube.exe -conf %s -ts %s 00:00:00 -te %s 23:59:30'
CONF   = r'conf_ppp.conf'

for i in range(len(doy)):
    t1 = datetime.datetime.now()
    
    ts = '%04d/%02d/%02d' % (year[i], month[i], day[i])
    cmd = CMD_EXE % (CONF, ts, ts)
    print(cmd)

    subprocess.run(cmd)
    
    t2 = datetime.datetime.now()
    print(f'({t1})-({t2})={t2-t1}\n')

    
