![JuliaSAFT_logo](./JuliaSAFT_logo.jpg)

# Methods

### The problem 

This document aims to outline all of the various tools used to obtain the relevant properties from a SAFT-type equation of state. In short, SAFT equations of state provide the Helmholtz free energy of a system at a given composition $\mathbf{z}$, volume $V$ and temperature $T$:
$$
A=A(\mathbf{z},V,T)
$$
 Taking derivatives of this function (within the JuliaSAFT module, this is done using automatic differentiation) can give us a wide range of properties which are given in the appendix. However, it is more common that we are interested in the state of a system at certain conditions ($\mathbf{z}_0$, $p_0$ , $T_0$). The answer to this can be determined from the following, deceptively simple, minimisation of the Gibbs free energy:
$$
\min G(\mathbf{z}_0,p_0,T_0)
$$
In the case of SAFT-type equations of state, this can be expressed as:
$$
\min A(\mathbf{z}_0,V,T_0)-p_0V
$$
What isn't obvious in this formulation of the problem is what the variables to be optimised are. Re-expressing this problem in greater detail:
$$
\min \sum_{i=1}^{n_\mathrm{phase}}\phi_i(A(\mathbf{z}_i,V_i,T_0)-p_0V_i)
$$

$$
\mathrm{s.t.} \left(\sum_{i=1}^{n_\mathrm{phases}}\phi_iz_{j,i}\right)-z_{j,0}=0\quad\forall j \in [1,n_\mathrm{species}]
$$

where the subscript $i$ denotes properties related to a phase $i$ and $\phi_i$ is the mole fraction of phase $i$. One can already see the difficulties behind solving such a problem as we do not often know before-hand how many phases there may be at the conditions ($\mathbf{z}_0$, $p_0$ , $T_0$) and thus, we won't know what variables to optimise for. In addition, we will want the global minimum and, particularly in systems with many components, there may be many local minima which we will need to eliminate. 

Nevertheless, if we know certain things about the system before-hand, we can reduce the problem to one that is easier to solve. 

### Pressure solvers

