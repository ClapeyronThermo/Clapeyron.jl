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

Derived by Chapman _et al._ (1990), this is the first variant of the SAFT equation of state. This equation can be seen as a `proof of concept' as not many parameters are available (none for mixtures). Nevertheless, some noteworthy features of this equation is its use of a semi-empirical equation to obtain the hard-sphere diameter that depends on the number of segments of a species (no other SAFT variant does this). The chain term uses the hard-sphere pair-distribution function, which has a much-simpler analytical form than what some other SAFT equations choose to use. The association strength, $\Delta$ is evaluated in a unique way as well:

``\Delta_{ij,ab}=d_{ij}^3g_{ij}^\mathrm{HS}F_{ij,ab}\kappa_{ij,ab}``

where $\kappa_{ij,ab}$ is dimensionless. Unfortunately, the implementation of `ogSAFT` in `OpenSAFT` cannot yet replicate the figures from the original paper. The reason for this is that the monomer / segment term presented in the paper is not the one used to generate the results. The actual term used is developed by Twu _et al._ (1980) and we are currently attempting to implement this within `OpenSAFT` but it is not clear, as of yet, how it was implemented within the original equation.

### CK-SAFT

If the SAFT equation derived by Chapman _et al._ was the prototype, the variant developed by Huang and Radosz (1990) was the first usable SAFT equation, with over a 100 pure-component parameters and many unlike parameters available. `CKSAFT` effectively simplifies many of the computationally-intensive parts of `ogSAFT`, using a simpler equation to obtain the hard-sphere diameter and actually providing the correct monomer term within the paper. The chain term between the two equations is identical. Similarly, the association strength only has a minor change:

``\Delta_{ij,ab}=\sigma_{ij}^3g_{ij}^\mathrm{HS}F_{ij,ab}\kappa_{ij,ab}``

which slightly reduces the computational cost. However, the most-noteworthy simplification came with the association term. As mentioned earlier, the association fraction needs to be solved for iteratively. However, Huang and Radosz proposed approximations of the association fraction that could be used to solve for the association term explicitly, greatly reducing the computational intensity of these calculations. These approximations have not been implemented within `OpenSAFT` as of yet, but these only impact calculations for species other than alcohols and carboxylic acids. We also point out that Huang and Radosz introduced the concept of association schemes which helps classify species based on how they interaction through association.

### SAFT-VR SW

Gil-Villegas _et al._ (1997) developed a new class of SAFT equations known as SAFT variable range. Here, more emphasis was placed on the potentials used to characterise dispersion interactions where a new parameter was introduced through the potential shape. Whilst many versions of SAFT-VR are proposed, each using different underlying potentials, the one that was chosen as the default was SAFT-VR square-well (SW) with the potential shape parameter $\lambda$ (characterising the width of the potential well). Within this framework, novel expressions for the monomer and chain terms were proposed, both being based on the SW potential. The association term remained largely unchanged, with the association term having the most-noteworthy modification:

``\Delta_{ij,ab}=g_{ij}^\mathrm{SW}F_{ij,ab}\kappa_{ij,ab}``

Here, $\kappa_{ij,ab}$ now carries units of volume. Not many parameters are available for this equation of state, primarily being use to model alkanes and perfluoro-alkanes. However, compared to most other SAFT variants, SAFT-VR SW has possibly seen the most extensions, having a group-contribution alternative (SAFT-$\gamma$ SW), electrolyte (SAFT-VRE SW) and cross-over theory (SAFT-VRX SW). 

### soft-SAFT

Developed by Blas and Vega (2001), whereas SAFT equations up until now have used a hard-sphere reference from which to build the equation of state, soft-SAFT chooses to use a Lennard-Jones reference instead. Because of this, compared to all other SAFT equations, soft-SAFT relies heavily on correlations obtained from molecular-dynamic simulations to obtain the monomer term, pair-distribution function and association strength. Like SAFT-VR SW, soft-SAFT does not have an extensive database of parameters, but has been extended multiple times (cross over theory being the more-noteworthy extension).

### PC-SAFT

Possibly the most-popular variant of the SAFT equation, Perturbed-Chain (not polymer-chain) SAFT was developed by Gross and Sadowski (2001) and, like soft-SAFT, chooses a different reference state than previous SAFT equations. This time, we start from the hard-chain (HC), not hard-sphere, expressing the SAFT equation as:

`` \frac{A}{Nk_\mathrm{B}T} = \frac{A_\mathrm{ideal}}{Nk_\mathrm{B}T}+\frac{A_\mathrm{HC}}{Nk_\mathrm{B}T}+\frac{A_\mathrm{disp.}}{Nk_\mathrm{B}T}+\frac{A_\mathrm{assoc.}}{Nk_\mathrm{B}T}``

This isn't as significant a change as one might initially think as, effectively, the hard-sphere and chain terms (which uses a hard-sphere pair distribution function like CK-SAFT) are combined into the hard chain term. The dispersion term is then simply another correlation, only this time depends on the number of segments as well. It carries many similarities with CK-SAFT, using the same expression for the hard-sphere diameter, pair-distribution function and association term.

The primary reason behind PC-SAFT's popularity is three-fold. For one, the code for PC-SAFT is available open-source. Secondly, there is an abundance of parameters available (over 250), including unlike parameters. Finally, many variants of the PC-SAFT equation have been developed. These include:

* Polar PC-SAFT (PPC-SAFT)
* PC-Polar SAFT (PCP-SAFT); yes, these are distinct equations
* Electrolyte PC-SAFT (ePC-SAFT)
* Electrolyte PPC-SAFT (ePPC-SAFT)
* Critical-point based PC-SAFT (CP-PC-SAFT)
* Critical-point based PPC-SAFT (CP-PPC-SAFT)
* Group-contribution PC-SAFT (GC-PC-SAFT)
* Group contribution PPC-SAFT (GC-PPC-SAFT)

We will aim to provide some of these variants at a later date.

#### sPC-SAFT

Nevertheless, we do provide one of these variants, being the simplified PC-SAFT equation (developed by Von Solms _et al._ (2003)). Here, the only modifications are to the hard-chain and association terms where, instead of using the generalised expressions for the hard-sphere term and hard-sphere pair-distribution function, by averaging the hard-sphere diameter (effectively treating mixtures as being made-up of identically-sized segments), the pure-component versions of these properties are used instead. The benefit of this is that pure-component parameters determined for PC-SAFT can still be used here, and only the unlike parameters need to be modified.

Similar to PC-SAFT, variants of the sPC-SAFT equation also exist, although no-where near as extensive. Most notably, a significant group-contribution method is available.

### SAFT-VR Mie

One of the most-novel SAFT equation of state, derived by Lafitte _et al._ (2013), this equation is effectively an extension of the SAFT-VR framework developed by Gil-Villegas _et al._ (1997), with further improvements. First of these is extending the Barker-Henderson perturbative expansion to third order instead of second order:

`` \frac{A_\mathrm{mono.}}{Nk_\mathrm{B}T}=\frac{A_\mathrm{HS}}{Nk_\mathrm{B}T}+\frac{A_\mathrm{1}}{(Nk_\mathrm{B}T)^2}+\frac{A_\mathrm{2}}{(Nk_\mathrm{B}T)^3}+\frac{A_\mathrm{3}}{(Nk_\mathrm{B}T)^4}``

We do point out that, whilst the first two terms are developed following the SAFT-VR framework, the third order term is more akin to a correlation regressed using molecular dynamic simulations of Mie fluids. This third order term resulted in significant improvements in the modelling of properties near the critical point (without using cross-over theory). The chain term also received further improvements as a result. This is also the only SAFT equation which evaluates the hard-sphere diameter analytically, although numerical approximations are needed (we note that the original SAFT-VR Mie equation used 10-point Gauss-Legendre quadrature, whilst the newer version uses 5-point Gauss-Laguerre quadrature).

However, three different versions of the association strength have been developed:

* Hard-sphere kernel:

  ``\Delta_{ij,ab}=\sigma_{ij}^3g_{ij}^\mathrm{HS}F_{ij,ab}K_{ij,ab}``

* Lennard-Jones kernel:

  ``\Delta_{ij,ab}=F_{ij,ab}K_{ij,ab}I_{ij,ab}(\epsilon_{ij},\sigma_{ij})``

* Mie kernel:

  ``\Delta_{ij,ab}=F_{ij,ab}K_{ij,ab}I_{ij,ab}(\epsilon_{ij},\sigma_{ij},\lambda_{ij})``

Unfortunately, it seems that there have been inconsistencies between which of these kernels is used in different publications. The current 'default' SAFT-VR Mie equation uses the Lennard-Jones kernel, as such, this is the one used in `OpenSAFT`. We do intend to provide the option to switch between these kernels.

As it uses a Mie potential is characterised by two shape parameters, $\lambda_\mathrm{a}$ (characterising the attractive part) and $\lambda_\mathrm{r}$ (characterising the repulsive part), both of these have become parameters for each species (although $\lambda_\mathrm{a}$ is usually set to 6). An interesting aesthetic change is with the number of segments where this is now separated into the shape factor, $S$, and the number of segments $v^*$. The latter must now be an integer and the former is a direct measure of how 'fused' the segments are. As we have different association terms, we also have different sets of parameters where the only difference is the length-scale. In the Lennard-Jones and Mie kernels, $K_{ij,ab}$ is the 'bonding volume', whereas, in the hard-sphere kernel, it is a 'bonding length', $r_{ij,ab}^c$.

The SAFT-VR Mie does not have a significantly large repository of parameters (compensated by its group-contribution variant) and has only been extended to electrolytes (SAFT-VRE Mie and eSAFT-VR Mie). 

#### SAFT-VRQ Mie

A very recent extension of the SAFT-VR Mie equation is the SAFT-VRQ Mie equation developed by Aasen _et al._ (2019) which modifies the underlying Mie potential using a Feynman-Hibbs potential, which means that a single species is represented by a sum of three Mie potentials. This method attempts to classically account for quantum effects present in small species such as helium, hydrogen and neon. Unfortunately, this equation is limited to just the monomer term and, even then, it is very computationally intensive. We do note that the current implementation in `OpenSAFT` can only model pure-component properties, but we will extend this to mixture in future versions.

### SAFT-$\gamma$ Mie

The group-contribution version of SAFT-VR Mie, developed by Papaioannou _et al._ (2014), the SAFT-$\gamma$ Mie equation uses the same general framework as SAFT-VR Mie, although, as it is a group-contribution method, we are able to model heterogenous chains (in previous SAFT equations, all segments in a chain were the same size). The group-contribution methodology is based on that developed by Lymperiadis _et al._ (2008). 37 groups are currently available for this equation. A noteworthy advantage of using groups is that unlike parameters between groups can be estimated from pure-component data; these can then be readily extended to mixtures without further regression. 

This equation has also been extended to electrolytes through SAFT-$\gamma$E Mie.