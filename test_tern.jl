model = UNIFAC(["water",("dichloromethane",["CH2CL2"=>1]),("acetone",["CH3"=>1,"CH3CO"=>1])])

N = 500
T = 298.15
p = 1e5
z0 = [0.5,0.5,1e-10]

K0 = [1000.,0.0001,0.001]
x_LLE = zeros(N,6)
idxend = N
(x,n,G) = tp_flash(model,p,T,z0,RRTPFlash(K0=K0,equilibrium=:lle))
x_LLE[1,1:3] = x[1,:]
x_LLE[1,4:6] = x[2,:]

K0 = x[1,:]./x[2,:]
z0[1:2] = x[1,1:2]/2+x[2,1:2]/2
z0[3] += 2/N
for i in 2:N
    global z0, K0, x_LLE
    (x,n,G) = tp_flash(model,p,T,z0,RRTPFlash(K0=K0,equilibrium=:lle))
    x_LLE[i,1:3] = x[1,:]
    x_LLE[i,4:6] = x[2,:]
    K0 = x[1,:]./x[2,:]

    z0 = (x_LLE[i,1:3]+x_LLE[i,4:6])-(x_LLE[i-1,1:3]+x_LLE[i-1,4:6])/2
    if abs(x[1,1]-x[2,1]) < 1e-3
        idxend=i-1
        break
    end
end

x_LLE = vcat(x_LLE[1:idxend,1:3],reverse(x_LLE[1:idxend,4:6],dims=1));