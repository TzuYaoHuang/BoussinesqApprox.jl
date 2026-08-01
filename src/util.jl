
BC!(a::AbstractArray{T,D}, perdir, N = size(a)) where {T,D} = for j ∈ 1:D
    if j∈perdir
        @loop a[I] = a[CIj(j,I,N[j]-1)] over I ∈ slice(N,1,j)
        @loop a[I] = a[CIj(j,I,2)] over I ∈ slice(N,N[j],j)
    else
        @loop a[I] = a[CIj(j,I,2)] over I ∈ slice(N,1,j)
        @loop a[I] = a[CIj(j,I,N[j]-1)] over I ∈ slice(N,N[j],j)
    end
end