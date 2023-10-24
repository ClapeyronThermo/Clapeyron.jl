using RDKitMinimalLib, ChemicalIdentifiers

res = search_chemical(component)
mol = get_mol(res.smiles)
mol_list = get_substruct_matches(mol,mol)

queries = get_qmol.(groups[:,1])

atoms = mol_list[1]["atoms"]
group_list = []
group_id = []
group_occ_list = []
atoms_list = []
coverage_atoms = []
coverage_bonds = []
for i in 1:length(groups[:,1])
    # println(atoms_list)
    if !isempty(get_substruct_match(mol,queries[i]))
        println(i)
        smatch = get_substruct_matches(mol,queries[i])
        for j in 1:length(smatch)
            # println(smatch[j]["atoms"])
            if isempty(atoms_list)
                append!(group_list,[groups[i,1]])
                append!(group_id,i)
                append!(group_occ_list,1)
                append!(atoms_list,smatch[j]["atoms"])
                append!(coverage_atoms,[smatch[j]["atoms"]])
                append!(coverage_bonds,[smatch[j]["bonds"]])
            else
                # If no atoms covered by this group are already covered by other groups
                if sum(smatch[j]["atoms"] .∈ [atoms_list])==0 
                    append!(atoms_list,smatch[j]["atoms"])
                    if !(groups[i,1] in group_list)
                        append!(group_list,[groups[i,1]])
                        append!(group_id,i)
                        append!(group_occ_list,1)
                        append!(coverage_atoms,[smatch[j]["atoms"]])
                        append!(coverage_bonds,[smatch[j]["bonds"]])
                    else
                        group_occ_list[end] += 1
                        append!(coverage_atoms[end],smatch[j]["atoms"])
                        append!(coverage_bonds[end],smatch[j]["bonds"])
                    end
                else 
                    # Check which groups group i has an overlap with
                    id = 0
                    ng_rm = 0
                    for k in 1:length(group_id)
                        id += 1
                        # Does group 1 cover any atoms of group id
                        if sum(smatch[j]["atoms"] .∈ [coverage_atoms[id]])>0
                            # We only care if group i covers _more_ atoms than group k
                            if ((length(smatch[j]["atoms"])>length(coverage_atoms[id])) & 
                                # Also make sure that group i covers all the atoms of group k 
                                (sum(smatch[j]["atoms"] .∈ [coverage_atoms[id]]).==length(coverage_atoms[id]))) |
                                (length(smatch[j]["bonds"])>length(coverage_bonds[id]))
                                # find out which atoms are covered
                                overlap_atoms = coverage_atoms[id][coverage_atoms[id] .∈ [smatch[j]["atoms"]]]
                                id_rm = group_id[id]
                                name_rm = group_list[id]
                                bond_rm =  coverage_bonds[id][coverage_bonds[id] .∈ [smatch[j]["bonds"]]]
                                filter!(e->e∉overlap_atoms,atoms_list)
                                filter!(e->e∉overlap_atoms,coverage_atoms[id])
                                group_occ_list[id] -= 1

                                # If group k no longer covers any atoms, remove it
                                if group_occ_list[id] == 0
                                    filter!(e->e≠0,group_occ_list)
                                    filter!(e->!isempty(e),coverage_atoms)
                                    deleteat!(coverage_bonds,id)
                                    filter!(e->e≠id_rm,group_id)
                                    filter!(e->e≠name_rm,group_list)
                                    id -= 1 
                                end
                                ng_rm +=1
                            end
                        end    
                    end
                    if ng_rm > 0
                        if !(groups[i,1] in group_list)
                            append!(group_list,[groups[i,1]])
                            append!(group_id,i)
                            append!(group_occ_list,1)
                            append!(coverage_atoms,[smatch[j]["atoms"]])
                            append!(coverage_bonds,[smatch[j]["bonds"]])
                            append!(atoms_list,smatch[j]["atoms"])
                        else
                            group_occ_list[end] += 1
                            append!(atoms_list,smatch[j]["atoms"])
                            append!(coverage_atoms[end],smatch[j]["atoms"])
                            append!(coverage_bonds[end],smatch[j]["bonds"])
                        end
                    end
                end
            end
        end
    end
end

if !(sum(atoms_list .∈ [atoms])==length(atoms))
    @warn "Could not find all groups for"*component
end

return (component,[groups[group_id[i],2] => group_occ_list[i] for i in 1:length(group_id)])