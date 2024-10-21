## Models

Here, we give a high-level description of equations of state and the models provided by `Clapeyron`.

### Equations of state

An equation of state is a functional form, $f$ (say), that allows us to obtain a thermodynamic property, $F$, at given conditions $\boldsymbol{\Omega}$:

``F = f(\boldsymbol{\Omega};\boldsymbol{\Xi})``.

There are many ways one can develop an equation of state, however, these must respect the constraints on the number of properties we can specify, as required by the Gibbs Phase Rule:

``\mathrm{DoF} = N_\mathrm{species} - N_\mathrm{phase} + 2``

Here, ``\mathrm{DoF}`` means "Degrees of Freedom"; this is the number of so-called intensive state properties (or, in other words, those that are independent of system size) that we can specify. As we can see, the largest number of degrees of freedom we can have is $N_\mathrm{species}+1$; system size itself is not included so, in practice, this represents one more variable that we can specify, giving $N_\mathrm{species}+2$. Thus, taking the simple case of a single species, we can specify at most three conditions in our system. For a traditional equation of state, we specify volume, $V$, temperature, $T$, and the size of the system – for example the number of particles, $N$, or moles, $n$; the equation of state then returns the pressure, $p$. Many modern equations of state are derived using what is known as the canonical ensemble (more information can be found in Statistical Mechanics textbooks) and, accordingly, the three variable chosen are again usually $T$, $V$ and $N$; the output of these equations is usually the Helmholtz free energy, $A$.

Many equations of state are based on an underlying molecular model. Consequently, it is also typical for an equation of state to require parameters, $\boldsymbol{\Xi}$, to model certain species. The nature of these parameters depends on the equation of state.

### Ideal-gas equation of state

One equation of state that most engineers and scientists should be very familiar with is the ideal-gas equation, commonly expressed as:

``pV = Nk_\mathrm{B}T``

where $p$ is the pressure, $N$ is the total number of particles and $k_\mathrm{B}$ is the Boltzmann constant. Ideal-gas molecules are "invisible" to each other; there is zero interaction between ideal-gas molecules (this is, in essence, the ideal-gas model), so you can think of them as being infinitely small. This equation was first written down (although in slightly different form) in 1834 by Émile Clapeyron (in whose honour `Clapeyron` is named). Most (if not all) subsequent equations of state are descended from Clapeyron's equation, which is surprisingly appropriate for a wide variety of species in the gas phase at high-enough temperature and low-enough pressure.

Unfortunately, if we wish to determine other thermodynamic properties this equation is a little inconvenient. It would be much easier if it was expressed in the form of the free energy, from which other properties can then be evaluated using standard thermodynamic relationships. Accordingly, we could first integrate the above equation with respect to volume to determine the Helmholtz free energy:

``A_\mathrm{ideal} =- \int p\,dV =- Nk_\mathrm{B}T\ln{V}+c(T,N)``

This is still a little inconvenient, however, since we have to deal with a tricky constant of integration. Fortunately, we can instead derive ``A_\mathrm{ideal}`` from statistical mechanics (using just a few well-known results from quantum mechanics). Following this route, we obtain (for a pure component (_i.e._, a single species))

``\frac{A_\mathrm{ideal}}{Nk_\mathrm{B}T} = ln{(\rho\Lambda^3)}-1``,

where $\rho = N/V$ is the number density, and $\Lambda$ is the thermal de Broglie wavelength, which introduces the kinetic contributions to the free energy (strictly speaking, with this notation only translations are included). We can generalise this as a sum over species $i$ for a multicomponent mixture:

``\frac{A_\mathrm{ideal}}{Nk_\mathrm{B}T} = \left(\sum_ix_i\ln{(\rho_i\Lambda_i^3)}\right)-1``,

where $x_i$ is the molar composition, and the subscript denote that the variable so decorated relates to species $i$. This equation represents the `MonomerIdeal` form in `Clapeyron`. For the purposes of vapour–liquid-equilibrium properties, one can even ignore $\Lambda_i$ (since it cancels out in solving the phase equilibrium); we therefore use this as the default `BasicIdeal` model.

The kinetic energy of polyatomic species includes contributions from vibrational and rotational modes of motion, as well as translational; we must also account for these in the ideal free energy. The statistical-mechanical derivation of the ideal free energy becomes a little more complicated but can still be done, resulting in the following expression:

``\frac{A_{\mathrm{ideal}}}{Nk_{\mathrm{B}}T}=\sum_{i}x_{i}\bigg[\ln\left(\rho_{i}\Lambda_{i}^{3}\right)-\frac{N_{\mathrm{rot},i}}{2} \ln \frac{T}{\theta_{\mathrm{rot},i}}+\sum^{N_{\mathrm{vib},i}}_{\mathrm{v}}g_{i,\mathrm{v}}\left[\frac{\theta_{\mathrm{vib},i,\mathrm{v}}}{2T}+\ln\left(1-\exp{-(\theta_{\mathrm{vib},i,\mathrm{v}}/T)}\right)\right]-1\bigg]``.

Here $N_{\mathrm{rot},i}$, $\theta_{\mathrm{rot},i}$ and $N_{\mathrm{vib},i}$ are the number of rotations, the number of vibrations and the rotational temperature of a species $i$, respectively; $g_{i,\mathrm{v}}$ and $\theta_{\mathrm{vib},i,\mathrm{v}}$ represent the degeneracy and vibrational temperature of a vibrational mode $\mathrm{v}$ of species $i$. The `WalkerIdeal` model provides the necessary parameters to use such an equation. However, the more-commonly used approach is through the use of correlations of the ideal isobaric heat capacity, $C_{p,i}^0$, such as the `ReidIdeal`, `WilhoitIdeal` and `AlyLeeIdeal` models. From the ideal isobaric heat capacity, it is possible to determine the ideal Helmholtz free energy using the following equation:

``\frac{A_{\mathrm{ideal}}}{Nk_\mathrm{B}T} = \sum_{i=1}^{N_{\mathrm{Component}}} x_i\left[\ln{\frac{\rho_i}{\rho_0}}
    + \frac{1}{Nk_\mathrm{B}T} \int_{T_0}^T \!\!C_{p,i}^0 dT + \frac{H_{0,i}}{Nk_\mathrm{B}T}- \frac{1}{Nk_{B}}\!\!\int_{T_0}^T \frac{C_{p,i}^0}{T} dT -\ln{\frac{T}{T_0}}-\frac{S_{0,i}}{Nk_\mathrm{B}} - 1\right]``

Note that the reference states, $\rho_0$, $H_{0,i}$ and $S_{0,i}$, can typically be neglected as these will not impact or contribute to most thermodynamic properties of interest.

### Cubic equations of state

This is the most-popular class of equations of state. The progenitor of these is the van der Waals (`vdW`) equation of state, published in 1873, which can be written as:

``p = \frac{Nk_\mathrm{B}T}{V-Nb}-\frac{N^2a}{V^2}``

where $a$ and $b$ are the model parameters. Although the vdW equation was phenomenological in origin, it, too, can be derived from statistical thermodynamics. Strictly speaking, $b$ accounts for the space taken up by the molecules themselves (it corresponds to the excluded volume per molecule) and $a$ quantifies the magnitude of attraction between species. As a result, the first term is often thought of as accounting for the repulsive interactions between molecules, while the second accounts for attractive interactions. In principle, therefore, one could obtain values of $a$ and $b$ for a particular species from (for example) spectroscopic information. However, since $a$ and $b$ can be related to the critical temperature and pressure of the vdW fluid, to relate the equation of state to a particular species, it is conventional to use the critical temperature and pressure of the species to obtain working values of the parameters. 

Unfortunately, although its simple functional form makes calculations quite straight-forward, the vdW equation is inadequate for quantitative modelling, particularly for volumetric properties, and is most useful only for providing a qualitative description of the thermodynamic properties of the fluid. 
As a result, many other engineering cubic equations have been developed, retaining (as far as possible) the simple mathematical form of van der Waals' equation. The first noteworthy one of these is the Redlich-Kwong (`RK`) equation:

``p = \frac{Nk_\mathrm{B}T}{V-Nb}-\frac{N^2a}{\sqrt{T}V(V+Nb)}``

There is no physical justification for the change in the second term; its origin is entirely empirical. The authors made the modification so that the equation would provide better gas-phase fugacities. This equation was subsequently improved upon by Soave, resulting in the `SRK` equation of state:

``p = \frac{Nk_\mathrm{B}T}{V-Nb}-\frac{N^2\alpha(T;\omega)}{V(V+Nb)}``

The $\alpha$ function requires an additional parameter, the acentric factor (or acentricity), which is effectively a measure of the location of the saturation pressure when $T/T_\mathrm{c}=0.7$, where $T_\mathrm{c}$ is the critical temperature. The idea behind this is, if you can capture both the critical point and another point along the vapour-pressure curve, you will improve the accuracy of your equation of state. This is indeed what happened. Although Soave described his equation as a "modified Redlich–Kwong equation", in truth it is more than this. The introduction of the $\alpha$ function represents a giant step forwards; the inclusion of a similar $\alpha$ function is probably the key feature in the equation of Peng and Robinson, who introduced their equation of state (`PR`) to provide improved liquid-phase volumetric properties. In addition to the inclusion of an $\alpha$ function, Peng and Robinson further revised the attractive term:

``p = \frac{Nk_\mathrm{B}T}{V-Nb}-\frac{N^2\alpha(T;\omega)}{V^2+2NbV+b^2N^2}``

The SRK and PR equations of state are comparable in performance, although the latter generally provides liquid densities with a greater degree of accuracy, while the former usually provides better fugacities. However, when it comes to modelling complex species such as polymers (macromolecules), or associating species, both equations struggle to perform well. This is unsurprising, since the underlying molecular model remains, in essence, a "van der Waalsian sphere" – in other words,a hard spherical core surrounded by a region of attraction. A more-sophisticated molecular model is required to account well for the increased molecular complexities of these species.

Before moving on from cubic equations of state we note that, within `Clapeyron`, the cubic plus association (`CPA`) equation of state is supported. A CPA equation is the amalgamation of a cubic equation (usually SRK, as in `Clapeyron`, or PR) with the association term from the SAFT equation, which we will meet later. Strictly speaking, it is neither a cubic nor a SAFT equation of state but, rather, occupies a middle ground between these two classes of equation. 

Something that may be apparent in all these equations is the fact that these are all functions that return the pressure and, thus, must be integrated to obtain the Helmholtz free energy. Like the ideal-gas equation, there will be missing temperature and compositional dependencies which need to be included.

#### Mixtures with cubic equations of state

One may wonder how to model mixtures using such equations. This can be achieved using _mixing rules_, in conjunction with _combining rules_. Although there are many variants, one of the more-popular mixing rules is the van der Waals one-fluid mixing rule: the mixture is treated as a hypothetical pure fluid, characterised by parameters $\bar{a}$ and $\bar{b}$ that are given by

``\bar{a}=\sum_i\sum_jx_ix_ja_{ij}``

``\bar{b}=\sum_i\sum_jx_ix_jb_{ij}``

When $i=j$, $a$ and $b$ are just the normal van der Waals parameters for the pure. However, when $i\neq j$, these parameter characterise the unlike interactions between $i$ and $j$. We typically need to use _combining rules_ (not to be confused with _mixing rules_) to determine the unlike parameters. Examples of these include:

``b_{ij}=\frac{b_i+b_j}{2}``

``a_{ij} = (1-k_{ij})\sqrt{a_ia_j}``

where $k_{ij}$ can be set to 0 but, using either more-advanced combining rules or regression to experimental data, can be tuned to improve the effectiveness of the combining rule. Further details on this will be given for the SAFT models.

More-complicated mixing rules (such as the Wong-Sandler mixing rule) are available and implemented in `Clapeyron`. 

### SAFT equations of state

In comparison to the cubic equations of state, equations based on the Statistical Associating Fluid Theory (SAFT) are based on a more-theoretical approach, although still can be considered as descendants of van der Waals' equation. As mentioned earlier, the van der Waals equation can be derived from statistical mechanics, whereby the Helmholtz free energy of the van der Waals fluid is obtained as

``\frac{A}{Nk_\mathrm{B}T} = \frac{A_\mathrm{ideal}}{Nk_\mathrm{B}T}+\frac{A_\mathrm{HS}}{Nk_\mathrm{B}T}+\frac{A_\mathrm{1}}{(Nk_\mathrm{B}T)^2}``;

here the ideal and hard-sphere (HS) terms combine to give the repulsive term (of the pressure form of the equation) whilst the $A_1$ term results in the attractive term. We can see from this that, using the van der Waals equation, species are effectively modelled as hard spheres with dispersive interactions (we sometimes call these London dispersion interactions). The latter two terms can be merged into what is referred to as the monomer or segment term. 

Whilst, as already noted, this is clearly a step up from the ideal gas, most species can't be modelled effectively as single spheres; they may be highly non-spherical in shape (as is usually the case with large molecules), or they may experience interactions that are more complex than simple dispersion. A classic example of the latter is water; although the water molecule is small and (at first glance) may appear simple, the behaviour of water is very strongly influenced by hydrogen-bonding interactions.  

Using Wertheim's TPT1 theory of association, it is possible to model molecules as chains of spheres; the shape of the model molecule can thereby be tailored to represent that of the real molecule far better than a single sphere. Wertheim's TPT1 theory can also be used to account for intermolecular association interactions (such as dipole–dipole interactions, or hydrogen bonding), which are strongly directional. These are described using associations sites that are located on one or more of the spherical segments comprising the chain molecule. This results in the addition of two extra contributions to the Helmholtz free energy (note that the HS and dispersive terms have been merged into a monomer term):

``\frac{A}{Nk_\mathrm{B}T} = \frac{A_\mathrm{ideal}}{Nk_\mathrm{B}T}+\frac{A_\mathrm{mono.}}{Nk_\mathrm{B}T}+\frac{A_\mathrm{chain}}{Nk_\mathrm{B}T}+\frac{A_\mathrm{assoc.}}{Nk_\mathrm{B}T}``

The chain term accounts for the formation of chains of spherical segments and is generally expressed as

``\frac{A_\mathrm{chain}}{Nk_\mathrm{B}T}=-\sum_ix_i(m_i-1)\ln{g_{ii}(d_{ii})}``,

where $g_{ii}(r_{ii})$ is the pair distribution function for species $i$ (which carries information about the structure of the fluid; it expresses how likely it is that a segment of species $i$ is present at a distance $r_{ii}$ from another segment of species $i$). Many SAFT equations differ in how this pair distribution function is expressed. We note here the introduction of the Barker-Henderson hard-sphere diameter, $d_{ii}$ which is given (dropping the subscripts for clarity) by

``d = \int_0^\sigma (1-\exp{-\beta\phi(r)})dr``;

here $\phi(r)$ is our _effective_ pair potential and $\beta=1/(k_\mathrm{B}T)$. This effectively gives a temperature dependence to the size of our segment and accounts for our segment becoming _softer_ as temperature rises.

The association term accounts for the highly-directional associative interactions (for example, hydrogen bonding). For most SAFT equations of state, it is expressed as:

``\frac{A_\mathrm{assoc.}}{Nk_\mathrm{B}T}=\sum_ix_i\left(\sum_a\left(\ln{X_{i,a}}-\frac{X_{i,a}}{2}\right)+\frac{M_i}{2}\right)``

where $X_{i,a}$ is the fraction of association sites $a$ on species $i$ _not_ bonded to another and $M_i$ is the number of association sites on species $i$. $X_{i,a}$ can be obtained by solving the following system of equations:

``X_{i,a} = (1+\rho\sum_jx_j\sum_bX_{j,b}\Delta_{ij,ab})^{-1}``

An important aspect of the association term is that the above system of equations typically needs to be solved iteratively; this greatly increases the computational cost of the SAFT equations when modelling associating species (compared to modelling non-associating species, for example, or to using cubic equations of state). $\Delta_{ij,ab}$ is the association strength between site $a$ on species $i$ with site $b$ on species $j$; this is also an aspect where SAFT equations usually differ but can all be written generally as

``\Delta_{ij,ab} = F_{ij,ab}K_{ij,ab}I_{ij,ab}``

where $F_{ij,ab}$ is Mayer's function, given by

``F = \exp{-\beta\epsilon^\mathrm{assoc.}}-1``,

where $\epsilon^\mathrm{assoc.}$ is the potential depth of the association interaction. $K$ and $I$ differ between equations but, generally, these represent the length scale of the interaction and the likelihood that the sites are correctly orientated such that they overlap, respectively.

Surprisingly, the monomer term is one of the aspects that most distinguishes the different SAFT equations; no two variants use the same equation. However, in general, the monomer term is composed of more than one term:

``\frac{A_\mathrm{mono.}}{Nk_\mathrm{B}T}=\frac{A_\mathrm{HS}}{Nk_\mathrm{B}T}+\frac{A_\mathrm{1}}{(Nk_\mathrm{B}T)^2}+\frac{A_\mathrm{2}}{(Nk_\mathrm{B}T)^3}+\frac{A_\mathrm{3}}{(Nk_\mathrm{B}T)^4}+...``

This expression is known as a Barker–Henderson perturbative expansion. These terms generally account for the dispersive interactions between segments; the $n^\mathrm{th}$ order term account for interactions between $n$ segments. In most SAFT equations, this expansion is truncated at just the second-order term.

#### Parameters

Although different SAFT equations use different parameters, most share a common set. These include the parameters that characterise the dispersive interactions (which are usually modelled as pair potentials): the potential depth $\epsilon$ (usually expressed as $\epsilon / k_{\mathrm{B}}$, in Kelvin) and the segment size $\sigma$ (in Angstrom). We point out here that this potential (and its parameters) is not a _bare_ pair potential, which accounts only for the interactions of two species (in vacuum); it is an _effective_ pair potential, which accounts for the effects of other species being around the interacting pair, in some cases quantum effects and, if associative interactions are not modelled separately, account for non-dispersive interactions.

As species can now be modelled as chains of segments, the number of segments, $m$, also becomes a parameter. One thing to point out about this parameter is that it need not be an integer (despite what its name suggests); non-integer values of $m$ can usually be interpreted as segments merging within the chain. 

In the case of associating species, SAFT equations usually require a parameter for the potential well depth of the association, $\epsilon^\mathrm{assoc.}$ (analogously to $\epsilon$, this is usually expressed as $\epsilon^\mathrm{assoc.} / k_{\mathrm{B}}$, in Kelvin) and a parameter characterising the length-scale of the interaction (either a bonding volume, $\kappa^\mathrm{assoc.}$, usually expressed in Å$^3$, or a length, $r_c^\mathrm{assoc.}$, either in meters or dimensionless (reduced by the segment size)). In the case of the dispersive and associative interaction parameters, there will also be the equivalent parameters characterising unlike interactions between species in a mixture (which can also be obtained from combining rules).

Unfortunately, due to the complex function form of SAFT equations, it is impossible to directly relate these parameters to critical properties, as is done with the engineering cubics. Instead, these parameters are typically obtained by regression using experimental data (typically pure-component saturation-pressure and saturated-liquid-density data).

We will next go through each of the variants of the SAFT equation available in `Clapeyron` and what makes these unique.

#### Original SAFT

Derived by Chapman _et al._ (1990), this is the first variant of the SAFT equation of state. This equation can be seen as a `proof of concept' as not many parameters are available (none for mixtures). Nevertheless, a noteworthy feature of this equation is the use of a semi-empirical equation to obtain the hard-sphere diameter that depends on the number of segments of a species (in no other SAFT variant is this done). The hard-sphere pair-distribution is used in the chain term; this has a much-simpler analytical form than what is chosen for use in some other SAFT equations. The association strength, $\Delta$ is evaluated in a unique way as well:

``\Delta_{ij,ab}=d_{ij}^3g_{ij}^\mathrm{HS}F_{ij,ab}\kappa_{ij,ab}``

where $\kappa_{ij,ab}$ is dimensionless. Unfortunately, the implementation of `ogSAFT` in `Clapeyron` cannot yet replicate the figures from the original paper. The reason for this is that the monomer / segment term presented in the paper is not the one used to generate the results. The actual term used is developed by Twu _et al._ (1980) and we are currently attempting to implement this within `Clapeyron` but it is not clear, as of yet, how it was implemented within the original equation.

#### CK-SAFT

If the SAFT equation derived by Chapman _et al._ was the prototype, the variant developed by Huang and Radosz (1990) was the first usable SAFT equation, with over 100 pure-component parameters and many unlike parameters available. Many of the computationally-intensive parts of `ogSAFT` are simplified in `CKSAFT`; a simpler equation is used to obtain the hard-sphere diameter, and the monomer term provided within the paper is the correct one. The chain term is identical in the two equations. Similarly, the association strength only has a minor change:

``\Delta_{ij,ab}=\sigma_{ij}^3g_{ij}^\mathrm{HS}F_{ij,ab}\kappa_{ij,ab}``;

this slightly reduces the computational cost. However, the most-noteworthy simplification came with the association term. As mentioned earlier, the association fraction needs to be solved for iteratively. However, Huang and Radosz proposed approximations of the association fraction that could be used to solve for the association term explicitly, greatly reducing the computational intensity of these calculations. These approximations have not been implemented within `Clapeyron` as of yet, but these only impact calculations for species other than alcohols and carboxylic acids. We also point out that Huang and Radosz introduced the concept of association schemes, which helps classify species based on how they interaction through association.

#### SAFT-VR SW

Gil-Villegas _et al._ (1997) developed a new class of SAFT equations known as SAFT variable range. Here, more emphasis was placed on the potentials used to characterise dispersion interactions, and a new parameter was introduced through the potential shape. Whilst many versions of SAFT-VR are proposed, each using different underlying potentials, the one that was chosen as the default was SAFT-VR square-well (SW) with the additional potential-shape-parameter $\lambda$ (characterising the width of the potential well). Within this framework, novel expressions for the monomer and chain terms were proposed, both being based on the SW potential. The association term remained largely unchanged, with the association strength having the most-noteworthy modification:

``\Delta_{ij,ab}=g_{ij}^\mathrm{SW}F_{ij,ab}\kappa_{ij,ab}``.

Here, $\kappa_{ij,ab}$ now carries units of volume. Not many parameters are available for this equation of state, primarily being used to model alkanes and perfluoro-alkanes. However, compared to most other SAFT variants, SAFT-VR SW has possibly seen the most extensions, having a group-contribution alternative (SAFT-$\gamma$ SW), electrolyte (SAFT-VRE SW) and cross-over theory (SAFT-VRX SW). 

#### soft-SAFT

Developed by Blas and Vega (2001), whereas the SAFT equations discussed up until now had been based on a hard-sphere reference from which to build the equation of state, in soft-SAFT a Lennard-Jones reference is used instead. Because of this, compared to all other SAFT equations, soft-SAFT relies heavily on correlations obtained from molecular-dynamics simulations to obtain the monomer term, pair distribution function and association strength. Like SAFT-VR SW, soft-SAFT does not have an extensive database of parameters, but has been extended multiple times (cross-over theory being the more-noteworthy extension).

#### PC-SAFT

Possibly the most-popular variant of the SAFT equation, Perturbed-Chain (not polymer-chain) SAFT was developed by Gross and Sadowski (2001) and, like soft-SAFT, a different reference state is chosen, as compared with previous SAFT equations. This time, we start from the hard-chain (HC), not the hard-sphere, expressing the SAFT equation as:

``\frac{A}{Nk_\mathrm{B}T} = \frac{A_\mathrm{ideal}}{Nk_\mathrm{B}T}+\frac{A_\mathrm{HC}}{Nk_\mathrm{B}T}+\frac{A_\mathrm{disp.}}{Nk_\mathrm{B}T}+\frac{A_\mathrm{assoc.}}{Nk_\mathrm{B}T}``

This isn't as significant a change as one might initially think as, effectively, the hard-sphere and chain terms (which uses a hard-sphere pair distribution function like CK-SAFT) are combined into the hard-chain term. The dispersion term is then simply another correlation, only this time depending on the number of segments as well. It carries many similarities with CK-SAFT, using the same expression for the hard-sphere diameter, pair distribution function and association term.

The primary reasons behind PC-SAFT's popularity are three-fold. For one, the code for PC-SAFT is available open-source. Secondly, there is an abundance of parameters available (over 250), including unlike parameters. Finally, many variants of the PC-SAFT equation have been developed. These include:

* Polar PC-SAFT (PPC-SAFT)
* PC-Polar SAFT (PCP-SAFT); yes, these are distinct equations
* Electrolyte PC-SAFT (ePC-SAFT)
* Electrolyte PPC-SAFT (ePPC-SAFT)
* Polyelectrolyte ePC-SAFT (epPC-SAFT)
* Critical-point based PC-SAFT (CP-PC-SAFT)
* Critical-point based PPC-SAFT (CP-PPC-SAFT)
* Group-contribution PC-SAFT (GC-PC-SAFT)
* Group contribution PPC-SAFT (GC-PPC-SAFT)

We will aim to provide some of these variants at a later date.

#### sPC-SAFT

We do already provide one of the PC-SAFT variants, namely the simplified PC-SAFT equation (developed by Von Solms _et al._ (2003)). Here, the only modifications are to the hard-chain and association terms where, instead of using the generalised expressions for the hard-sphere term and hard-sphere pair distribution function, by averaging the hard-sphere diameter (effectively treating mixtures as being made up of identically sized segments), the pure-component versions of these properties are used instead. The benefit of this is that pure-component parameters determined for PC-SAFT can still be used here, and only the unlike parameters need to be modified.

Similar to PC-SAFT, variants of the sPC-SAFT equation also exist, although nowhere near as extensive. Most notably, a significant group-contribution method is available.

#### SAFT-VR Mie

One of the most-novel SAFT equations of state, derived by Lafitte _et al._ (2013), this equation is effectively an extension of the SAFT-VR framework developed by Gil-Villegas _et al._ (1997), with further improvements. First of these is extending the Barker–Henderson perturbative expansion to third order instead of second order:

``\frac{A_\mathrm{mono.}}{Nk_\mathrm{B}T}=\frac{A_\mathrm{HS}}{Nk_\mathrm{B}T}+\frac{A_\mathrm{1}}{(Nk_\mathrm{B}T)^2}+\frac{A_\mathrm{2}}{(Nk_\mathrm{B}T)^3}+\frac{A_\mathrm{3}}{(Nk_\mathrm{B}T)^4}``

We do point out that, whilst the first two terms are developed following the SAFT-VR framework, the third-order term is more akin to a correlation regressed using molecular-dynamics simulation data for Mie fluids. This third-order term resulted in significant improvements in the modelling of properties near the critical point (without using cross-over theory). The chain term also received further improvements as a result. This is also the only SAFT equation in which the hard-sphere diameter is evaluated analytically, although numerical approximations are needed (we note that the original SAFT-VR Mie equation used 10-point Gauss-Legendre quadrature, whilst the newer version uses 5-point Gauss-Laguerre quadrature).

However, three different versions of the association strength have been developed:

* Hard-sphere kernel:

  ``\Delta_{ij,ab}=\sigma_{ij}^3g_{ij}^\mathrm{HS}F_{ij,ab}K_{ij,ab}``

* Lennard-Jones kernel:

  ``\Delta_{ij,ab}=F_{ij,ab}K_{ij,ab}I_{ij,ab}(\epsilon_{ij},\sigma_{ij})``

* Mie kernel:

  ``\Delta_{ij,ab}=F_{ij,ab}K_{ij,ab}I_{ij,ab}(\epsilon_{ij},\sigma_{ij},\lambda_{ij})``

Unfortunately, it seems that there have been inconsistencies between which of these kernels is used in different publications. In the current 'default' SAFT-VR Mie equation the Lennard-Jones kernel is used; as such, this is the one used in `Clapeyron`. We do intend to provide the option to switch between these kernels.

As a Mie potential is characterised by two shape parameters, $\lambda_\mathrm{a}$ (characterising the attractive part) and $\lambda_\mathrm{r}$ (characterising the repulsive part), both of these have become parameters for each species (although $\lambda_\mathrm{a}$ is usually set to 6). As we have different association terms, we also have different sets of parameters where the only difference is the length-scale. In the Lennard-Jones and Mie kernels, $K_{ij,ab}$ is the 'bonding volume', whereas, in the hard-sphere kernel, it is a 'bonding length', $r_{ij,ab}^c$.

The SAFT-VR Mie does not have a significantly large repository of parameters (compensated by its group-contribution variant) and has only been extended to electrolytes (SAFT-VRE Mie and eSAFT-VR Mie). 

#### SAFT-VRQ Mie

A very recent extension of the SAFT-VR Mie equation is the SAFT-VRQ Mie equation developed by Aasen _et al._ (2019) in which the underlying Mie potential is modified using a Feynman-Hibbs potential, which means that a single species is represented by a sum of three Mie potentials. This method attempts to classically account for quantum effects present in small species such as helium, hydrogen and neon. Unfortunately, this equation is limited to just the monomer term and, even then, it is very computationally intensive. We do note that the current implementation in `Clapeyron` can only model pure-component properties, but we will extend this to mixtures in future versions.

#### SAFT-$\gamma$ Mie

The group-contribution version of SAFT-VR Mie, developed by Papaioannou _et al._ (2014), the SAFT-$\gamma$ Mie equation rests on the same general framework as SAFT-VR Mie, although, as it is a group-contribution method, we are able to model heterogeneous chains (in SAFT equations discussed previously, all segments in a chain were the same size). The group-contribution methodology is based on that developed by Lymperiadis _et al._ (2008). An interesting aesthetic change is with the number of segments where this is now separated into the shape factor, $S$, and the number of segments $v^*$. The latter must now be an integer and the former is a direct measure of how 'fused' the segments are. Approximately 60 groups are currently available for this equation. A noteworthy advantage of using groups is that unlike parameters between groups can be estimated from pure-component data; these can then be readily extended to mixtures without further regression.

This equation has also been extended to electrolytes through SAFT-$\gamma$E Mie.

## Methods

### The problem

The aim of this document is to outline all of the various tools used to obtain the relevant properties from a SAFT-type equation of state. In short, SAFT equations of state provide the Helmholtz free energy of a system at a given composition $\mathbf{z}$, volume $V$ and temperature $T$:

``A=A(\mathbf{z},V,T)``

Taking derivatives of this function (within the `Clapeyron` module, this is done using automatic differentiation) can give us a wide range of properties which are given in the appendix. However, it is more common that we are interested in the state of a system at certain conditions ($\mathbf{z}_0$, $p_0$ , $T_0$). The answer to this can be determined from the following, deceptively simple, minimisation of the Gibbs free energy:

``\min G(\mathbf{z}_0,p_0,T_0)``

In the case of SAFT-type equations of state, this can be expressed as:

``\min A(\mathbf{z}_0,V,T_0)+p_0V``

What isn't obvious in this formulation of the problem is how to identify the variables that are to be optimised. Re-expressing this problem in greater detail:

``\min \sum_{i=1}^{n_\mathrm{phase}}\phi_i(A(\mathbf{z}_i,V_i,T_0)+p_0V_i)``

``\mathrm{s.t.} \left(\sum_{i=1}^{n_\mathrm{phases}}\phi_iz_{j,i}\right)-z_{j,0}=0\quad\forall j \in [1,n_\mathrm{species}]``

where the subscript $i$ denotes properties related to a phase $i$, and $\phi_i$ is the molar fraction of phase $i$. One can already see the difficulties behind solving such a problem as we do not often know beforehand how many phases there may be at the conditions ($\mathbf{z}_0$, $p_0$ , $T_0$) and thus, we won't know which variables to optimise. In addition, we will want the global minimum and, particularly in systems with many components, there may be many local minima that we will need to eliminate.

Nevertheless, if we know certain things about the system beforehand, we can reduce the problem to one that is easier to solve.

### Volume solvers

Let us make one simplifying assumption: we know that the system exists in a single phase. This greatly simplifies the problem to:

``\min_V A(\mathbf{z}_0,V,T_0)+p_0V``

where, as there is no phase split, the only variable we need to optimise is the volume. We can see that this is equivalent to solving for the volume at which the pressure predicted by the equation of state equals $p_0$:

``\min_V A(\mathbf{z}_0,V,T_0)+p_0V\rightarrow\frac{\partial }{\partial V}(A(\mathbf{z}_0,V,T_0)+p_0V)=\frac{\partial A}{\partial V}(\mathbf{z}_0,V,T_0)+p_0=-p(\mathbf{z}_0,V,T_0)+p_0=0``

Effectively, we can re-word this as a root-finding problem. One slight issue with this is that there is often more than one root (there can actually be up to five, even in SAFT-type equations). The true root will be the one that minimises the Gibbs free energy; thus we must first find the candidate phases and determine their Gibbs free energy before reporting the volume.

For the cubics, this problem is quite straightforward given that (as the name suggests) all these equations can be rearranged as a cubic equation in $V$:

``p_0=\frac{RT_0}{V-b}-\frac{a}{(V-c_1)(V-c_2)}\rightarrow a_0+a_1V+a_2V^2+a_3V^3=0``.

Thus, it is very easy to solve for all the roots in a cubic using analytical expressions. However, for other equations of state, we must use non-linear root-finding algorithms. In order to avoid the unstable phases, we try to use initial guesses close to what will be the 'true' phases:

1. Vapour: Since we can use automatic differentiation to obtain the virial coefficient for any model, we can actually obtain an initial guess extremely close to the final solution using

   ``\frac{p_0}{RT_0} = \frac{n_0}{V}+\frac{n_0}{V^2}B(T)\rightarrow V_0=\frac{RT_0}{p_0}\frac{-1+\sqrt{1+4p_0B(T_0)/(RT_0)}}{2}``.

2. Liquid: The best we can do here is to obtain the volume corresponding to a large packing fraction (we typically pick 0.6-0.8):
   ``V_0 = \frac{N_\mathrm{A}\pi}{6\times 0.8}m\sigma^3``.
   We are still looking for ways to improve this but the volume function is quite reliable as of now.

One other issue to consider when solving this problem is that, within the liquid phase, the gradients are very large; this can be difficult for algorithms to handle (even when providing the exact derivatives through automatic differentiation). We try to reduce magnitude of these derivatives by solving for the logarithm of the volume instead. 
