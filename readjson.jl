import JSON

dict2 = Dict()
open("all_data.json", "r") do f
    global dict2
    dict2=JSON.parse(f)  # parse and transform data
end

println(dict2)
