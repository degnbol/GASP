#!/usr/bin/env julia
using DataFrames

"""
Left join where duplicate names will neither raise an error nor be suffixed, 
but rather values in df2 will replace missing entries in df1, and raise an 
error if modifications occur to non-missing entries.
Another difference to DataFrames.leftjoin is the on keyword can contain column 
names not found between df1 and df2, as long as at least one on column name is 
shared.
"""
function joinleft(df1::DataFrame, df2::DataFrame, on)
    _on = intersect(String.(on), names(df1), names(df2))
    @assert length(_on) > 0 "No on columns shared: $on"
    df = leftjoin(df1, df2; on=_on, makeunique=true)
    # makeunique will then append _1 and _2 to duplicate column names
    nms = names(df1)
    nms = nms[nms .* "_1" .∈ Ref(names(df))]
    for nm in nms
        nm1 = nm*"_1"
        ismiss = ismissing.(df[!, nm])
        df[ismiss, nm] .= df[ismiss, nm1]
        nomiss = .!ismiss .& .!ismissing.(df[!, nm1])
        if !all(df[nomiss, nm] .== df[nomiss, nm1])
            modifications = df[nomiss, [_on..., nm, nm1]] |> unique
            error("joinleft has to modify:\n" * repr(modifications))
        end
        select!(df, Not(nm1))
    end
    df
end
joinleft(df1::DataFrame, df2::DataFrame, on::Union{AbstractString,Symbol}) = joinleft(df1, df2, [on])

function joinleft(dfs::Vector{DataFrame}, on)
    df = dfs[1]
    for df2 in dfs[2:end] # handles edge case length(df) == 1.
        df = joinleft(df, df2, on)
    end
    df
end

