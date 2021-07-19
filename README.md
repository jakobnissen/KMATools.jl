# KMATools.jl

Package for parsing various files produced by *k-mer* alignment progran KMA. Tested on KMA 1.3.22.

See docstrings for functions `parse_spa`, `parse_res` and `parse_mat`.

The functions in this package work on `IO`-objects. To work with a file, do e.g.:
```julia
path = "my_path.spa"
spa = open(path) do io
    parse_spa(io, path)
end
```

For gzipped files, you can do:
```julia
using CodecZlib
path = "my_path.mat.gz"
mat = open(GzipDecompressorStream, path) do io
    parse_mat(io, path)
end
```