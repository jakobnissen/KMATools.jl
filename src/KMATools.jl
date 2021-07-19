"""
    KMATools

Package for parsing various files produced by KMA. Tested on KMA 1.3.22.
"""
module KMATools

using BioSymbols

imap(f) = x -> Iterators.map(f, x)
ifilter(f) = x -> Iterators.filter(f, x)

const SPA_HEADER = join(
    [
        "#Template",
        "Num",
        "Score",
        "Expected",
        "Template_length",
        "Query_Coverage",
        "Template_Coverage",
        "Depth",
        "tot_query_Coverage",
        "tot_template_Coverage",
        "tot_depth",
        "q_value",
        "p_value",
    ],
    '\t'
)

"""
    parse_spa(io::IO, path::String) -> Vector{NamedTuple}

Parse .spa file, returning a vector of `NamedTuple` with the following names:
`template, num, score, expected, tlen, qcov, tcov, depth, total_qcov, total_tcov,`
`total_depth, qval, pval`.

Coverages and identities are represented as fractions (`Float64`) in [0.0, 1.0]
"""
function parse_spa(io::IO, path::String)
    counted_lines = enumerate(eachline(io))
    header = let
        it = iterate(counted_lines)
        header = if it === nothing
            error("Missing header from file \"$path\"")
        else
            last(first(it))
        end
        if strip(header) != SPA_HEADER
            error("Malformed header in file \"$path\"")
        end
        header
    end
    counted_lines |> 
        imap(i -> (first(i), strip(last(i)))) |>
        ifilter(i -> !isempty(last(i))) |>
        imap() do (line_number, line)
            fields = split(line, '\t')
            length(fields) != 13 && error("Expected 13 fields in file \"$path\", line \"$line_number\"")
            return (;
                template    = String(fields[1]),
                num         = parse(UInt, fields[2]),
                score       = parse(UInt, fields[3]),
                expected    = parse(UInt, fields[4]),
                tlen        = parse(UInt, fields[5]),
                qcov        = round(parse(Float64, fields[6]) / 100, digits=6),
                tcov        = round(parse(Float64, fields[7]) / 100, digits=6),
                depth       = parse(Float64, fields[8]),
                total_qcov  = round(parse(Float64, fields[9]) / 100, digits=6),
                total_tcov  = round(parse(Float64, fields[10]) / 100, digits=6),
                total_depth = parse(Float64, fields[11]),
                qval        = parse(Float64, fields[12]),
                pval        = parse(Float64, fields[13]),
            )
        end |>
    collect
end

const RES_HEADER = join(
    [
        "#Template",
        "Score",
        "Expected",
        "Template_length",
        "Template_Identity",
        "Template_Coverage",
        "Query_Identity",
        "Query_Coverage",
        "Depth",
        "q_value",
        "p_value",
    ],
    '\t'
)

"""
    parse_res(io::IO, path::String) -> Vector{NamedTuple}

Parse .res file, returning a vector of `NamedTuple` with the following names:
`template, score, expected, tlen, tid, tcov, qid, qcov, depth qval, pval`.

Coverages and identities are represented as fractions (`Float64`) in [0.0, 1.0]
"""
function parse_res(io::IO, path::String)
    counted_lines = enumerate(eachline(io))
    header = let
        it = iterate(counted_lines)
        header = if it === nothing
            error("Missing header from file \"$path\"")
        else
            last(first(it))
        end
        if strip(header) != RES_HEADER
            error("Malformed header in file \"$path\"")
        end
        header
    end
    counted_lines |>
        imap(i -> (first(i), strip(last(i)))) |>
        ifilter(i -> !isempty(last(i))) |>
        imap() do (line_number, line)
            fields = split(line, '\t')
            length(fields) != 11 && error("Expected 11 fields in file \"$path\", line \"$line_number\"")
            return (;
                template = String(fields[1]),
                score    = parse(UInt, fields[2]),
                expected = parse(UInt, fields[3]),
                tlen     = parse(UInt, fields[4]),
                tid      = round(parse(Float64, fields[5]) / 100, digits=6),
                tcov     = round(parse(Float64, fields[6]) / 100, digits=6),
                qid      = round(parse(Float64, fields[7]) / 100, digits=6),
                qcov     = round(parse(Float64, fields[8]) / 100, digits=6),
                depth    = parse(Float64, fields[9]),
                qval     = parse(Float64, fields[10]),
                pval     = parse(Float64, fields[11]),
            )
        end |>
    collect
end

function split!(v::Vector{SubString{String}}, s::Union{String, SubString{String}}, sep::UInt8)
    n = 0
    start = 1
    @inbounds for i in 1:ncodeunits(s)
        if codeunit(s, i) == sep
            n += 1
            n >= length(v) && throw(BoundsError(v, n+1))
            substr = SubString(s, start, i-1)
            v[n] = substr
            start = i + 1
        end
    end
    @inbounds v[n+1] = SubString(s, start, ncodeunits(s))
    n + 1 < length(v) && error("not enough fields for split!")
    v
end

"""
    parse_map(io::IO, path::String) -> Vector{Tuple{String, Vector{Row}}}

Parse .mat file, returning a vector of sequence matrices. A sequence matrix
is a `Tuple{String, Vector{Row}}` with the sequence name as the first part. The `Row`
is `Tuple{DNA, NTuple{6, UInt32}}`, with the 6 depths being the number of
A, C, G, T, N and gap, respectively.
"""
function parse_mat(io::IO, path::String)
    rowT = Tuple{DNA, NTuple{6, UInt32}}
    result = Vector{Tuple{String, Vector{rowT}}}()
    current = Vector{rowT}()
    fields = Vector{SubString{String}}(undef, 7)
    linedepths = Vector{UInt32}(undef, 6)
    header = nothing
    for line in eachline(io) |> imap(strip) |> ifilter(!isempty)
        if startswith(line, '#')
            if !isempty(current)
                @assert header isa String
                push!(result, (header, copy(current)))
                empty!(current)
            end
            header = String(line[2:end])
            continue
        end
        isnothing(header) && error("Expected header in file \"$path\"")
        split!(fields, line, UInt8('\t'))
        if ncodeunits(first(fields)) != 1
            error("Multi-character reference nucleotide in file \"$path\"")
        end
        refnuc = DNA(first(first(line)))
        depth_tuple = ntuple(i -> parse(UInt32, @inbounds fields[i+1]), Val(6))
        @inbounds for i in 1:6
            linedepths[i] = parse(UInt32, fields[i+1])
        end
        push!(current, (refnuc, depth_tuple))
    end
    isempty(current) || push!(result, (header, current))
    return result
end

export parse_spa, parse_res, parse_mat

end # module
