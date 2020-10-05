# Methods

## The problem

This document aims to outline all of the various tools used to obtain the relevant properties from a SAFT-type equation of state. In short, SAFT equations of state provide the Helmholtz free energy of a system at a given composition $\mathbf{z}$, volume $V$ and temperature $T$:

``A=A(\mathbf{z},V,T)``

 Taking derivatives of this function (within the OpenSAFT module, this is done using automatic differentiation) can give us a wide range of properties which are given in the appendix. However, it is more common that we are interested in the state of a system at certain conditions ($\mathbf{z}_0$, $p_0$ , $T_0$). The answer to this can be determined from the following, deceptively simple, minimisation of the Gibbs free energy:

``\min G(\mathbf{z}_0,p_0,T_0)``

In the case of SAFT-type equations of state, this can be expressed as:

``\min A(\mathbf{z}_0,V,T_0)+p_0V``

What isn't obvious in this formulation of the problem is what the variables to be optimised are. Re-expressing this problem in greater detail:

``\min \sum_{i=1}^{n_\mathrm{phase}}\phi_i(A(\mathbf{z}_i,V_i,T_0)+p_0V_i)``

``\mathrm{s.t.} \left(\sum_{i=1}^{n_\mathrm{phases}}\phi_iz_{j,i}\right)-z_{j,0}=0\quad\forall j \in [1,n_\mathrm{species}]``

where the subscript $i$ denotes properties related to a phase $i$ and $\phi_i$ is the molar fraction of phase $i$. One can already see the difficulties behind solving such a problem as we do not often know before-hand how many phases there may be at the conditions ($\mathbf{z}_0$, $p_0$ , $T_0$) and thus, we won't know what variables to optimise for. In addition, we will want the global minimum and, particularly in systems with many components, there may be many local minima which we will need to eliminate.

Nevertheless, if we know certain things about the system before-hand, we can reduce the problem to one that is easier to solve.

## Pressure solvers

Let us make one simplifying assumption: we know that the system exists in a single phase. This greatly simplifies the problem to:

``\min_V A(\mathbf{z}_0,V,T_0)+p_0V``

Where, as there is no phase split, the only variable we need to optimise is the volume. We can see that this is equivalent to solving for the volume at which the pressure predicted by the equation of state equals $p_0$:

``\min_V A(\mathbf{z}_0,V,T_0)+p_0V\rightarrow\frac{\partial }{\partial V}(A(\mathbf{z}_0,V,T_0)+p_0V)=\frac{\partial A}{\partial V}(\mathbf{z}_0,V,T_0)+p_0=-p(\mathbf{z}_0,V,T_0)+p_0=0``

Effectively, we can re-word this as a root-finding problem. When using the van der Waals or engineering equations of state which are expressed as $p(\mathbf{z},V,T)$, it is easier to solve them this way. When there is only one candidate phase, there is no significant advantage between expressing the problem as either an optimisation or root-finding problem.

However, there will be a range of pressures below the critical temperature where there will be more than one candidate phase (corresponding to the vapour, liquid and unstable phases). Treating this as a root-finding problem has the added difficulty of there being an additional, unstable solution. Treating this as an optimisation problem means we never need to worry about this unstable phase (it corresponds to a local maxima).

Actually determining the values of $V$ that minimise this equation is quite straightforward, although, with a few subtleties. Within OpenSAFT, we have used the local, derivative-based method of moving assymptotes (MMA) algorithm as implemented in `NLopt.jl` module. The reason for selecting this method is because, as a local derivative-based algorithm, it will be faster than other methods. This algorithm in particular also allows us to add inequality constraints; this is particularly important as there are certain values of $V$ which will are unphysical. These can be identified through the packing fraction:

``\eta=\frac{N_\mathrm{A}\pi}{6V}\sum_ix_im_id_i^3``

Without going into significant detail about the SAFT equation itself, a packing fraction greater than or equal to one results in unphysical values. As  a result, we have a lower bound for the volume:

``V\geq\frac{N_\mathrm{A}\pi}{6}\sum_ix_im_id_i^3``

One other issue to consider when solving this problem is that, within the liquid phase, the gradients are very large which can be difficult for algorithms to handle (even when providing the exact derivatives through automatic differentiation). There are two solutions to this:

1. Good initial guesses: We can re-express the volume of a system in dimensionless units using a variable commonly used in the SAFT equations of state, the packing fraction:

   ``\eta=\frac{N_\mathrm{A}\pi}{6V}m\sigma^3``

   Where, for most fluids, we expect the packing fraction to be close to $0.9$ in the liquid phase. Thus, if we know that the system is within the liquid phase, we can use the initial guess:

   ``V_0 = \frac{N_\mathrm{A}\pi m\sigma^3}{6\times0.9}``

   It should also be mentioned that, if we know that the system is within the vapour phase which typically has a packing fraction of $10^{-3}$, we can use the initial guess of:

   ``V_0 = \frac{N_\mathrm{A}\pi m\sigma^3}{6\times10^{-3}}``

   Generally speaking, if we know which phase our system is in, we can use the initial guesses to find the volumes $V$ that correspond to that phase, if the phase exists at the given conditions.
2. Solving for $\log_{10}{V}$ rather than $V$, i.e.:

   ``\min_x A(\mathbf{z}_0,10^x,T_0)+p_010^x``

   This somewhat reduces the magnitude of the gradients in the liquid phase.

Using the above tricks, one should be able to obtain the value of $V$ that minimises the Gibbs free energy. The only question to answer now is: if there is more than one local minima, how do we identify the stable phase? In this case, we need to use a global optimisation algorithm. In the case of OpenSAFT, a tunneling algorithm has been implemented although any other such algorithms can be used; the tunneling algorithm was selected as it still relies on gradient-based methods and is generally the recommended algorithm for such problems.
