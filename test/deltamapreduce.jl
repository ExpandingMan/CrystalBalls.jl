#==============================================================================================
Theory here is to consider any function R[f] as a mapreduce.

To compute only part of reduce that altered, invert reducer and then re-apply
as seen below...
==============================================================================================#


function deltamapreduce(R, v::Array{T}, mapper::Function, reducer::Function, idx, x::T) where T
    reducer(inv(reducer)(R, mapper(v[idx])), mapper(x))
end
