
a_res(model::IdealModel,V,T,z) = zero(Base.promote_eltype(model,V,T,z))

function volume_impl(model::IdealModel,p,T,z,phase,threaded,vol0)
    return sum(z)*R̄*T/p
end

lb_volume(model::IdealModel,z) = zero(eltype(z))

idealmodel(model::IdealModel) = model

@newmodelsingleton ZeroIdeal IdealModel
a_ideal(::ZeroIdeal,V,T,z) = zero(Base.promote_eltype(V,T,z))

#just for completion
function eos_g(model::IdealModel,p,T,z)
    R = Rgas(model)
    RT = R*T
    n = sum(z)
    V = n*RT/p
    return n*RT*(a_ideal(model,V,T,z) + 1)
end

function _update_idealuserlocations_for_GC(idealmodel::IdealModel,ideal_userlocations,Mw::SingleParam)
   type_idealmodel = typeof(idealmodel)
   return _update_idealuserlocations_for_GC(type_idealmodel,ideal_userlocations,Mw)
end

function _update_idealuserlocations_for_GC(idealmodel::Type{<:IdealModel},ideal_userlocations,Mw::SingleParam)
   if isempty(ideal_userlocations)
      if hasfield(idealmodel,:params)
         paramtype = fieldtype(idealmodel,:params)
         ref_state = hasfield(paramtype,:reference_state)
         nparams = length(fieldnames(paramtype))
         if hasfield(paramtype,:Mw) && nparams - ref_state == 1 # only models with Mw and reference state or just Mw
            ideal_userlocations = (;Mw=Mw.values)
         end
      end
   end
   return ideal_userlocations
end