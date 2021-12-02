"""
tp_flash(model, p, T, n, numphases;  MaxSteps = 1e4*(numphases-1), PopulationSize = 50, TraceMode = :silent)
        MaxSteps = 1e4*(numphases-1), PopulationSize = 50, TraceMode = :silent)

Routine to solve non-reactive multicomponent flash problem by finding global minimum of Gibbs Free Energy via Differential Evolution. User must assume a number of phases, numphases. If true number of phases is smaller than numphases, model should predict either (a) identical composition in two or more phases, or (b) one phase with negligible total number of moles. If true number of phases is larger than numphases, a thermodynamically unstable solution will be predicted.

Inputs: 
 - T, Temperature
 - p, Pressure
 - n, vector of number of moles of each species
 - numphases, assumed number of phases  

Outputs - Tuple containing: 
 - x_ij, Array of mole fractions of species j in phase i
 - n_ij, Array of mole numbers of species j in phase i, mol
 - G, Gibbs Free Energy of Equilibrium Mixture
"""

function tp_flash(model::EoSModel, p, T, n, numphases; 
        MaxSteps = 1e4*(numphases-1), PopulationSize = 50, TraceMode = :silent)
    
    numspecies = length(model.components)
    
    if numspecies != length(n)
        error("There are ", numspecies,
            " species in the model, but the number of mole numbers specified is ", 
            length(n))
    end
    
    if numphases == 1
        return (n, n / sum(n), gibbs_free_energy(model, p, T, n))
    end
    
    #=
    Function to calculate Gibbs Free Energy for given partition of moles between phases.
    This is a little tricky. 
    We must find a way of uniquely converting a vector of numbers,
    each in (0, 1), to a partition. We must be careful that 
    the mapping is 1-to-1, not many-to-1, as if many inputs
    map to the same physical state in a redundant way, there
    will be multiple global optima, and the global optimization 
    will perform poorly.
    Our approach is to specify (numphases-1) numbers in (0,1) for
    each species. We then scale these numbers systematically in order to partition
    the species between the phases. Each set of (numphases - 1) numbers
    will result in a unique partition of the species into the numphases
    phases.
    =#
    
    function GibbsFreeEnergy(dividers)
    
        dividers = reshape(dividers, 
            (numphases - 1, numspecies))
        
        #Initialize arrays xij and nvalsij, 
        #where i in 1..numphases, j in 1..numspecies
        #xij is mole fraction of j in phase i.
        #nvals is mole numbers of j in phase i.
        x = zeros(numphases, numspecies)
        nvals = zeros(numphases, numspecies)
        
        #Calculate partition of species into phases
        for j = 1:numspecies
            nvals[1,j] = dividers[1,j]*n[j]
            for i = 2:numphases - 1
                nvals[i,j] = dividers[i,j]*(n[j] - sum(nvals[1:(i-1), j]))
            end
            nvals[numphases, j] = n[j] - sum(nvals[1:(numphases-1), j])
        end
        
        #Calculate mole fractions xij
        for i = 1:numphases
            x[i, :] = nvals[i, :] / sum(nvals[i, :])
        end      
        
        #Calculate Overall Gibbs Free Energy (J/mol)
        #If any errors are encountered, return 1e9, ensuring point is discarded
        #by DE Algorithm
        G = 0.0
        for i = 1:numphases
            try
                Gphase = gibbs_free_energy(model, p, T, x[i, :]) * sum(nvals[i, :])
                if Gphase == NaN
                    return 1e9
                end
                G += Gphase
            catch
                return 1e9
            end
        end
        
        return G
    end
    
    #Minimize Gibbs Free Energy
    result = bboptimize(GibbsFreeEnergy; SearchRange = (0.0, 1.0), 
        NumDimensions = numspecies*(numphases-1), MaxSteps=MaxSteps, PopulationSize = PopulationSize, 
        TraceMode = TraceMode)
    
   dividers = reshape(best_candidate(result), 
            (numphases - 1, numspecies))
        
    #Initialize arrays xij and nvalsij, 
    #where i in 1..numphases, j in 1..numspecies
    #xij is mole fraction of j in phase i.
    #nvals is mole numbers of j in phase i.
    x = zeros(numphases, numspecies)
    nvals = zeros(numphases, numspecies)

    #Calculate partition of species into phases
    for j = 1:numspecies
        nvals[1,j] = dividers[1,j]*n[j]
        for i = 2:numphases - 1
             nvals[i,j] = dividers[i,j]*(n[j] - sum(nvals[1:(i-1), j]))
        end
        nvals[numphases, j] = n[j] - sum(nvals[1:(numphases-1), j])
    end

    #Calculate mole fractions xij
    for i = 1:numphases
        x[i, :] = nvals[i, :] / sum(nvals[i, :])
    end      
    
    return (x, nvals, best_fitness(result))
end

export tp_flash
