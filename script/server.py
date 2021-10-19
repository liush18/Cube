

import subprocess, datetime


def network():
    for i in range(len(dow)):
        t1 = datetime.datetime.now()
        cmd = CMD_EXE % (week[i], dow[i])
        print(cmd)
        subprocess.run(cmd)
        t2 = datetime.datetime.now()
        print(f'({t1})-({t2})={t2-t1}\n')

    
dow  = [0]
week = [2144 for i in range(7)]

EXE = r'Cube.exe -conf conf_net.conf '
CMD_EXE = EXE + r'-week %04d -dow %d'
network()




