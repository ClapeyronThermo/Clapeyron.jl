---
title: 'OpenSAFT: A Julia implementation of SAFT-type equations of state'
tags:
  - Julia
  - thermodynamics
  - molecular modelling of fluids
  - SAFT
  - equations of state
authors:
  - name: Hon-Wa Yew^[Corresponding author]
    affiliation: 1
  - name: Andrés Riedemann
    affiliation: 3
  - name: Pierre J. Walker
    affiliation: 2
affiliations:
 - name: Independent researcher
   index: 1
 - name: Imperial College London, Department of Chemical Engineering, South Kensington Campus, SW7 2AZ, London, U.K.
   index: 2
 - name: University of Concepción, Department of Chemical Engineering, Concepción, Región del Bio Bio, Chile
   index: 3
date: 01 October 2020
---

# Statement of need
Thermodynamic equations of state (EoS) remain an important tool in a chemical engineer's toolbox, playing a key role in process and molecular modelling software. Engineering cubics, such as the Soave-Redlich-Kwong and Peng-Robinson EoS, are often the go-to EoS for chemical engineers. However, with the development of more advanced processes, particularly in the area of carbon-capture, such equations of state are inadequate due to the complexity of the interactions between the molecules involved. It is for this reason that statistical association fluid theory (SAFT) equations of state have grown in popularity, with the ability to explicitly model complex phenomena such as hydrogen bonding and the flexibility to be extended to model other systems (such as polymers and electrolytes).

Nevertheless, in spite of this, SAFT-type equations of state remain largely inaccessible to most chemical engineers, partly due to the complexity of the equations themselves and the methods required to use them, but also due to the lack of open-source software providing these tools. This is worsened by the abundance of SAFT variants, from the original SAFT equation, to the popular PC-SAFT and to the most-recent SAFT-$\gamma$ Mie, adding a barrier-to-entry in the usage of these equations. Most commercial tools that provide SAFT equations of state are typically only compatible with one or two versions of the equation (examples include ASPEN with PC-SAFT and gPROMS with SAFT-$\gamma$ Mie).

It is for the reasons previously listed that OpenSAFT was developed, in the hope of not only providing the variants of the SAFT equation of state and the required methods to use it, but to also increase the understanding and accessibility of this equation amongst chemical engineers.

# Summary

# Examples
probably make a nice VLE figure comparing all of the available EoS.

# Acknowledgements

The authors thank doctors Andrew J. Haslam, Spiros Kornopoulos, Felipe A. Hurtado, Sara Febra and professors Amparo Galindo, Georgios M. Kontogeorgis, Xiaodong Liang and George Jackson for introducing us to the topic of thermodynamic equations of state and their continued help in the development of this package.
