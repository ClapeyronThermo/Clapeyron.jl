import numpy as np, json

mmin, mmax = 1, 64
ymin, ymax = 1/mmax, 1/mmin # y=1/m
expsL, expsV = [], []
N = 16 # degree in 1/m direction
k = np.arange(0,N+1,1)
ynodes = (ymax-ymin)/2*np.cos(k*np.pi/N) + (ymin+ymax)/2
mnodes = 1/ynodes

coefL = []
coefV = []
xmin = []
xmax = []

for m in mnodes:
    Ex = json.load(open(f'output/PCSAFT_VLE_m{m:0.12e}_expansions.json'))
    AL = []
    AV = []
    min = []
    max = []

    i = 0
    for ex in Ex:
        AL.append(ex['coefL'])
        AV.append(ex['coefV'])
        min.append(ex['xmin'])
        max.append(ex['xmax'])
        i+=1
    coefL.append(AL)
    coefV.append(AV)
    xmin.append(min)
    xmax.append(max)
print(xmin)
