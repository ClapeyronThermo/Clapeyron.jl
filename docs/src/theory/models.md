# Models

Here, we give a high-level description of equations of state and the models provided by OpenSAFT.

## Equations of state

Equations of state provide a functional form to obtain a thermodynamic property, $F$, at given conditions $\boldsymbol{\Omega}$:

``F = f(\boldsymbol{\Omega};\boldsymbol{\Xi})``

where $f$ is the equation of state. There are many ways one can develop an equation of state, however, the most commonly used approach is through what is known as the canonical ensemble. More information on this can be found in Statistical Mechanics textbooks. This typically results in equations that determine the Helmholtz free energy, $A$, at a given temperature, $T$, system volume, $V$ and system composition, $\mathbf{N}$. It is also typical for an equation of state to require parameters, $\boldsymbol{\Xi}$, to model certain species. What these parameters are depend on the equation of state.

## Ideal model

One equation of state that most engineers and scientists should be very familiar with is the ideal gas equation, commonly expressed as:

``pV = Nk_\mathrm{B}T``

where $p$ is the pressure, $N$ is the total number of particles and $k_\mathrm{B}$ is the Boltzmann constant. It is possible to derive this equation by assuming that species can be modelled as infinitesimally  small particles. This is surprisingly appropriate for a variety of species in the gas phase at high temperature and low pressure.

However, if we wish to determine other thermodynamic properties, we would need to integrate the above equation with respect to volume to determine the Helmholtz free energy:

``A_\mathrm{ideal} =- \int p\,\mathrm{d}V =- Nk_\mathrm{B}T\ln{V}+c(T,\mathbf{N})``

We can see that, in just using the original equation, we've lost a temperature and composition dependence in the Helmholtz free energy. If we follow the derivation from statistical or quantum mechanics, we can obtain the following equation (we denote this as the `Monomer` model in OpenSAFT):

``\frac{A_\mathrm{ideal}}{Nk_\mathrm{B}T} = \left(\sum_ix_i\ln{(\rho_i\Lambda_i^3)}\right)-1``

where $x_i$, $\rho_i$, and $\Lambda_i$ are the molar composition, number density and thermal de Broglie wavelength of species $i$, respectively. For the purposes of vapour-liquid equilibrium and most thermodynamic properties, this ideal model is sufficient (one can even ignore $\Lambda_i$ in the case of the former; we use this as the default `Basic` model). In addition, some would argue that this is the only ideal model as it only considers the translational motion of infinitesimally small particles.

However, polyatomic species typically have vibrational and rotational modes of motion as well which are typically included as part of the ideal term. These can also be derived from statistical mechanics giving:

``\frac{A_{\mathrm{ideal}}}{Nk_{\mathrm{B}}T}=\sum_{i}x_{i}\bigg[\ln\left(\rho_{i}\Lambda_{i}^{3}\right)-\frac{N_{\mathrm{rot},i}}{2} \ln \frac{T}{\theta_{\mathrm{rot},i}}+\sum^{N_{\mathrm{vib},i}}_{\mathrm{v}}g_{i,\mathrm{v}}\left[\frac{\theta_{\mathrm{vib},i,\mathrm{v}}}{2T}+\ln\left(1-\exp{-(\theta_{\mathrm{vib},i,\mathrm{v}}/T)}\right)\right]-1\bigg]``

where $N_{\mathrm{rot},i}$, $\theta_{\mathrm{rot},i}$ and $N_{\mathrm{vib},i}$ are the number of rotations, vibrations and rotational temperature of a species $i$, respectively. $g_{i,\mathrm{v}}$ and $\theta_{\mathrm{vib},i,\mathrm{v}}$ are the degeneracy and vibrational temperature of a vibrational mode $\mathrm{v}$ on species $i$, respectively. The `Walker` model provides the necessary parameters to use such an equation. However, the more-commonly used approach is through the use of ideal isobaric heat-capacity, $C_{p,i}^0$, correlations, such as the `Reid`, `Wilhoit` and `Aly and Lee` models. With the ideal isobaric heat-capacity, it is possible to determine the ideal Helmholtz free energy using the following equation:

``\frac{A_{\mathrm{ideal}}}{Nk_\mathrm{B}T} = \sum_{i=1}^{N_{\mathrm{Component}}} x_i\left[\ln{\frac{\rho_i}{\rho_0}}
    + \frac{1}{Nk_\mathrm{B}T} \int_{T_0}^T \!\!C_{p,i}^0 dT + \frac{H_{0,i}}{Nk_\mathrm{B}T}- \frac{1}{Nk_{B}}\!\!\int_{T_0}^T \frac{C_{p,i}^0}{T} dT -\ln{\frac{T}{T_0}}-\frac{S_{0,i}}{Nk_\mathrm{B}} - 1\right]``

Note that the reference states, $\rho_0$, $H_{0,i}$ and $S_{0,i}$, can typically be neglected as these will not impact or contribute to most thermodynamic properties of interest.



## Cubic models

Some of the more-popular equations of state have been the engineering cubic equations. The first of these is the van der Waals (`vdW`) equation of state, written as:

``p = \frac{Nk_\mathrm{B}T}{V-Nb}-\frac{N^2a}{V^2}``

where $a$ and $b$ are the model parameters which can be related to the critical temperature and pressure of a species. Although this equation was originally empirical, it is possible to derive this equation from statistical thermodynamics where $b$ corresponds to the excluded volume of a single species and $a$ quantifies the magnitude of attraction between species. As a result, the first term typically accounts for the repulsive interactions between species and the second accounts for attractive interactions. Although its simple functional form makes calculations quite straight-forward, this model is inadequate for modelling the liquid phase and vapour-liquid equilibrium properties. 

As result, wanting to keep with the van der Waals equation's simple form, a few engineering cubic equations have been developed. The first noteworthy one of these is the Redlich-Kwong (`RK`) equation:

``p = \frac{Nk_\mathrm{B}T}{V-Nb}-\frac{N^2a}{\sqrt{T}V(V+Nb)}``

There is no physical justification for the change in the second term, however, it was found to give improved modelling of the liquid phase. This equation was subsequently improved upon by Soave, resulting in the `SRK` equation of state:

 ``p = \frac{Nk_\mathrm{B}T}{V-Nb}-\frac{N^2\alpha(T;\omega)}{V(V+Nb)}``

The $\alpha$ function requires an additional parameter, the acentricity factor, which is effectively a measure of the location of the saturation pressure when $T/T_c=0.7$. The idea behind this is, if you can capture both the critical point and another point along the vapour curve, you will improve the accuracy of your equation of state. This is indeed what happened. Further improvements were also made by Peng and Robinson who introduced their own equation of state (`PR`):

``p = \frac{Nk_\mathrm{B}T}{V-Nb}-\frac{N^2\alpha(T;\omega)}{V^2+2NbV+b^2N^2}``

Both the SRK and PR equations of state are comparable in performance, although the latter generally models liquid densities to a greater degree of accuracy. However, when it comes to modelling complex species such as chains or associating species, both models tend to perform badly. We do note that, within OpenSAFT,  the cubic plus association (`CPA`) equation of state has been labelled as a cubic equation; whilst this is the case, we will describe it in greater detail when discussing the SAFT-type models as it main improvement over other cubics is borrowed from the SAFT theory.

Something that may be apparent in all these equations is the fact that these are all functions that give a pressure and, thus, must be integrated to obtain the Helmholtz free energy. Like the ideal gas equation, there will be missing temperature and compositional dependences which need to be included.

One may also wonder how to model mixtures using such equations; this can be achieved using _mixing rules_. Although there are many variants, one of the more-popular ones is the van der Waals one-fluid mixing rules which treats the mixture as having the same parameters $\bar{a}$ and $\bar{b}$ which can be determined from:

``\bar{a}=\sum_i\sum_jx_ix_ja_{ij}``

``\bar{b}=\sum_i\sum_jx_ix_jb_{ij}``

More-complicated mixing rules do exist (such as the Wong-Sandler mixing rule) which will be made available in OpenSAFT. When $i=j$, $a$ and $b$ are just the normal van der Waals parameters for the pure. However, when $i\neq j$, these parameter characterise the unlike interactions between $i$ and $j$. We typically need to use _combining rules_ (not to be confused with _mixing rules_) to determine the unlike parameters. Examples of these include:

``b_{ij}=\frac{b_i+b_j}{2}``

``a_{ij} = (1-k_{ij})\sqrt{a_ia_j}``

where $k_{ij}$ can be set to 0 but, using either more-advanced combining rules or regression to experimental data, can be tuned to improve the effectiveness of the combining rule. Further details on this will be given for the SAFT models.

## SAFT models

In comparison to the cubic equations of state, equations based on the Statistical Associating Fluid Theory (SAFT) take a more-theoretical approach. As mentioned earlier, the van der Waals equation can be derived from statistical mechanics where the resultant Helmholtz free energy is given by:

``\frac{A}{Nk_\mathrm{B}T} = \frac{A_\mathrm{ideal}}{Nk_\mathrm{B}T}+\frac{A_\mathrm{HS}}{Nk_\mathrm{B}T}+\frac{A_\mathrm{1}}{(Nk_\mathrm{B}T)^2}``

where the ideal and hard-sphere (HS) terms combine to give the repulsive term whilst the $A_1$ term results in the attractive term. We can see that, in the van der Waals equation, species are effectively modelled as hard-spheres with dispersive interactions (we sometimes can these London dispersion interactions). The last two terms can be merged into what is referred to as the monomer or segment term. Whilst this is a step up from the ideal term, most species can't be modelled effectively as single spheres and, in cases like water, experiences interactions more complex that simple dispersion (dipoles and hydrogen bonding). 

Using Wertheim's TPT1 theory of association, it is possible to model species as chains which can interact through both dispersive and associative interactions. The latter is described as interactions through associations sites on the segments which are strong and highly directional (like hydrogen bonding and dipole interactions). This results in the addition of two extra contributions to the Helmholtz free energy (note that the HS and dispersive terms have been merged into a monomer term):

`` \frac{A}{Nk_\mathrm{B}T} = \frac{A_\mathrm{ideal}}{Nk_\mathrm{B}T}+\frac{A_\mathrm{mono.}}{Nk_\mathrm{B}T}+\frac{A_\mathrm{chain}}{Nk_\mathrm{B}T}+\frac{A_\mathrm{assoc.}}{Nk_\mathrm{B}T}``

The chain term accounts for the formation of chains of spherical segments and is generally expressed as:

``\frac{A_\mathrm{chain}}{Nk_\mathrm{B}T}=-\sum_ix_i(m_i-1)\ln{g_{ii}(d_{ii})}``

where $g_{ij}(r_{ij})$ is the pair-distribution function (i.e. the likelihood of a segment of species $i$ being present at a distance $r$ from another segment of species $j$) for a hard-sphere. Many SAFT equations differ in how to express this pair-distribution function. We note here the introduction of the Barker-Henderson hard-sphere diameter, $d$ which is given by:

``d = \int_0^\sigma (1-\exp{-\beta\phi(r)})dr``

where $\phi(r)$ is our _effective_ pair potential and $\beta=1/(k_\mathrm{B}T)$. This effectively gives a temperature dependence to the size of our segment and accounts for our segment becoming _softer_ as temperature rises.

The association term accounts for the highly-directional associative interactions. For most SAFT equations of state, it is expressed as:

``\frac{A_\mathrm{assoc.}}{Nk_\mathrm{B}T}=\sum_ix_i\left(\sum_a\left(\ln{X_{i,a}}-\frac{X_{i,a}}{2}\right)+\frac{M_i}{2}\right)``

where $X_{i,a}$ is the fraction of association sites $a$ on species $i$ not bonded to another and $M_i$ is the number of association sites on species $i$.  $X_{i,a}$ can be solved for using the following system of equations:

``X_{i,a} = (1+\rho\sum_jx_j\sum_bX_{j,b}\Delta_{ij,ab})^{-1}``

An important aspect of the association term is that the above equation results in a system of equations that typically needs to be solved iteratively; this greatly increases the computational cost of the SAFT equations. $\Delta_{ij,ab}$ is the association strength between site $a$ on species $i$ with site $b$ on species $j$; this is also an aspect where SAFT equations usually differ but can all be written generally as:

``\Delta_{ij,ab} = F_{ij,ab}K_{ij,ab}I_{ij,ab}``

Where $F_{ij,ab}$ is Mayer's function given by:

``F = \exp{-\beta\epsilon^\mathrm{assoc.}}-1``

where $\epsilon^\mathrm{assoc.}$ is the potential depth of the association interaction. $K$ and $I$ differ between equations but, generally, these represent the length scale of the interaction and the likelihood that the sites are correctly orientated such that they overlap, respectively.

Surprisingly, the monomer term is one of the aspects that most distinguishes the different SAFT equations where no two variants use the same equation. However, in general, the monomer term is composed of more than one term:

`` \frac{A_\mathrm{mono.}}{Nk_\mathrm{B}T}=\frac{A_\mathrm{HS}}{Nk_\mathrm{B}T}+\frac{A_\mathrm{1}}{(Nk_\mathrm{B}T)^2}+\frac{A_\mathrm{2}}{(Nk_\mathrm{B}T)^3}+\frac{A_\mathrm{3}}{(Nk_\mathrm{B}T)^4}+...``

This is known as a Barker-Henderson perturbative expansion; the $N^\mathrm{th}$ order terms account for interactions between $N$ segments. Most SAFT equations truncate this sum at just the second order term. These terms generally account for the dispersive interactions between segments.

#### Parameters

Although different SAFT equations use different parameters, most share a common set. These include the parameters that characterise the dispersive interactions (which are usually modelled as pair potentials): the potential depth $\epsilon$ and the segment size $\sigma$. We point out here that this potential (and its parameters) is not a _bare_ pair potential which only accounts for the interactions of two species; it is an _effective_ pair potential which accounts for the effects of other species being around the interacting pair, in some cases quantum effects and, if associative interactions are not modelled separately, account for non-dispersive interactions.

As species can now be modelled as chains of segments, the number of segments, $m$, also becomes a parameter. One thing to point out about this parameter is it need not be an integer (despite what its name suggest); non-integer values of $m$ can usually be interpreted as segments merging within the chain. 

Only required for associating species, SAFT equations usually require a parameter for the potential well depth of the association, $\epsilon^\mathrm{assoc.}$ and a parameter characterising the length-scale of the interaction (either a bonding volume, $\kappa^\mathrm{assoc.}$ or length, $r_c^\mathrm{assoc.}$). In the case of the dispersive and associative interaction parameters, there will also be the equivalent parameters characterising unlike interactions between species in a mixture (which can also be obtained from combining rules).

Unfortunately, due to the complex function form of SAFT equations, it is impossible to directly relate these parameters to critical properties like in the engineering cubics. These parameters are typically obtained by regression using experimental data (typical pure-component saturation pressure and saturated liquid density data).

We will next go through each of the variants of the SAFT equation available in OpenSAFT and what makes these unique.

### Original SAFT

### CK-SAFT

### SAFT-VR SW

### soft-SAFT

### PC-SAFT

#### sPC-SAFT

### SAFT-VR Mie

#### SAFT-VRQ Mie

### SAFT-$\gamma$ Mie

