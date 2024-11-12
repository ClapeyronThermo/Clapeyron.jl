function flash(specifier,model,v1,v2,z,args...;kwargs...)

end


function ph_flash(model,p,h,z,T0 = Tproperty(model,p,h,z))

end

function ps_flash(model,p,s,z,T0 = Tproperty(model,p,s,z,entropy))
    
end

