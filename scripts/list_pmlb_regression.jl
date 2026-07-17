using PythonCall

dataset_lists = pyimport("pmlb.dataset_lists")
pmlb = pyimport("pmlb")

df = dataset_lists.df_summary
regression = df[df["task"].eq("regression")]
names         = pyconvert(Vector{String}, regression["dataset"].tolist())
n_instances   = pyconvert(Vector{Int},    regression["n_instances"].tolist())
n_features    = pyconvert(Vector{Int},    regression["n_features"].tolist())
n_continuous  = pyconvert(Vector{Int},    regression["n_continuous_features"].tolist())

order = sortperm(n_instances)
names        = names[order]
n_instances  = n_instances[order]
n_features   = n_features[order]
n_continuous = n_continuous[order]

println(rpad("dataset", 50), lpad("n_instances", 12), lpad("n_features", 12), lpad("n_continuous", 14))
println("-"^88)
for (name, ni, nf, nc) in zip(names, n_instances, n_features, n_continuous)
    # try 
    #     pmlb.fetch_data(name)
    # catch e 
    #     println("Skipping $name $e")
    #     continue
    # end
    println(rpad(name, 50), lpad(ni, 12), lpad(nf, 12), lpad(nc, 14))
end
println("\nTotal: $(length(names)) regression datasets")
