using CSV, DataFrames

gt = CSV.read(
    "FebruaryRestart/GlobalTherm_upload_02_11_17.csv",
    DataFrame;
    missingstring = "",
    ntasks = 1,     # avoid multithreaded inference
    pool = false    # ← CRITICAL: disable string pooling
)

for name in names(gt)
    col = gt[!, name]
    if eltype(col) <: Union{Missing, AbstractString}
        gt[!, name] = passmissing(String).(col)
    end
end
