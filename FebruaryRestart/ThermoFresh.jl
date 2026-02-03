using CSV, DataFrames

tf = CSV.read(
    "FebruaryRestart/thermtol_taxonomy_final.csv",
    DataFrame;
    missingstring = "",
    ntasks = 1,
    pool = false
)
tf_nomiss = tf[.!ismissing.(tf.class), :]

tf_fish = tf_nomiss[tf_nomiss.class .== "Actinopteri", :]

species = unique(tf_fish.species)  # Vector{String}

tf_df = CSV.read(
    "FebruaryRestart/thermtol_comb_final.csv",
    DataFrame;
    missingstring = "",
    ntasks = 1,
    pool = false
)
