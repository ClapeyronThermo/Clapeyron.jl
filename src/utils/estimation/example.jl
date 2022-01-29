using Clapeyron, BlackBoxOptim

model = SAFTVRMie(["helium"])

e = Estimation(model)

function saturation_p(model,T)
    sat = saturation_pressure(model,T)
    return sat[1]
end

function saturation_ρL(model,T)
    sat = saturation_pressure(model,T)
    return 1/sat[2]
end

function crit_T(model)
    sat = crit_pure(model)
    return sat[1]
end

function error(model,property,cond,data)
    if size(cond,2)==3
        pred = property.(model,cond[:,1],cond[:,2],cond[:,3])
    elseif size(cond,2)==2
        pred = property.(model,cond[:,1],cond[:,2])
    elseif size(cond,2)==1
        pred = property.(model,cond[:,1])
    else
        pred = property.(model)
    end
    return sum(abs.((pred.-data)./data))/length(data)
end

function estimate(estimator,params,guesses,properties,cond,data)
    model = return_model(estimator, params, guesses)
    
    return sum([properties[i][1]*error(model,properties[i][2],cond[i],data[i]) for i ∈ 1:length(properties)])
end

Exp = [2.1768	5039.3	36480
2.2768	6382.1	36433
2.3768	7951.9	36346
2.4768	9768.8	36228
2.5768	11852	36086
2.6768	14221	35923
2.7768	16896	35741
2.8768	19894	35543
2.9768	23235	35329
3.0768	26938	35100
3.1768	31021	34855
3.2768	35502	34594
3.3768	40401	34318
3.4768	45736	34024
3.5768	51526	33713
3.6768	57789	33383
3.7768	64546	33032
3.8768	71816	32659
3.9768	79620	32260
4.0768	87980	31834
4.1768	96917	31375
4.2768	106460	30880
4.3768	116620	30341
4.4768	127450	29750
4.5768	138960	29094
4.6768	151190	28357
4.7768	164190	27512
4.8768	178010	26514
4.9768	192690	25281
5.0768	208330	23612
5.1768	225070	20441];

properties=[(2.,saturation_p),(2.,saturation_ρL),(1.,crit_T)]
cond = [Exp[:,1],Exp[:,1],[0.,0.,0.,0.]']
data = [Exp[:,2],Exp[:,3],[5.1953]]
params = [:epsilon,:sigma,:lambda_r]
guesses = [[3.],[3.5],[16.]]

f(x) = estimate(e,params,[[x[1]],[x[2]*1e-10],[x[3]]],properties,cond,data)

x=bboptimize(f; SearchRange = [(3.7, 5.), (3.3,3.8), (12.,16)],MaxTime=1000);