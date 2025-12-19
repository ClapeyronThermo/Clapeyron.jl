using PlotlyLight
using PlotlyKaleido

_strip_quotes(s::AbstractString) = (startswith(s, "\"") && endswith(s, "\"")) ? replace(s[2:end-1], "\"\"" => "\"") : s

function _csv_split(line::AbstractString)
    fields = String[]
    buf = IOBuffer()
    in_quotes = false

    i = firstindex(line)
    while i <= lastindex(line)
        c = line[i]
        if c == '"'
            if in_quotes
                # escaped quote: "" -> "
                j = nextind(line, i)
                if j <= lastindex(line) && line[j] == '"'
                    write(buf, '"')
                    i = j
                else
                    in_quotes = false
                end
            else
                in_quotes = true
            end
        elseif c == ',' && !in_quotes
            push!(fields, String(take!(buf)))
        else
            write(buf, c)
        end
        i = nextind(line, i)
    end
    push!(fields, String(take!(buf)))
    return fields
end

function _parse_tpflash_algo_population(algo::AbstractString)
    m = match(r"(\d+)$", algo)
    m === nothing && error("cannot infer population size from algo name: $algo (expected trailing digits, e.g. rdex40)")
    pop = parse(Int, m.captures[1])
    family = replace(algo, r"(\d+)$" => "")
    return family, pop
end

"""
    plot_tpflash_summary_csv(path; outpath=nothing, title=nothing, theme=:light, bar_width=nothing, bargap=nothing, bargroupgap=nothing)

Plot a TPFlash stage2 CSV (either `*_summary.csv` or `*_per_case.csv`) as:
- x axis: population size (parsed from `algo`, e.g. `rdex40 -> 40`)
- y1 axis (left): mean U-score (`mean_u_norm`)
- y2 axis (right): overall success rate (`overall_success_rate`)

If multiple algorithm families exist in the same CSV (e.g. `rdex40`, `bbo20`),
each family is plotted as its own pair of traces (U-score + success rate).

Bar width tuning:
- `bar_width` sets the per-trace bar width. With numeric x (e.g. 20, 30, 40), the unit is the x-axis unit.
  For example, if x points are spaced by 10, `bar_width=6` gives visible gaps between groups.
- `bargap` / `bargroupgap` adjust spacing between bars (layout-level). These work best with categorical x.

If `outpath` is provided, saves the plot via PlotlyKaleido.
Returns the PlotlyLight `Plot`.
"""
function plot_tpflash_summary_csv(
    path::AbstractString;
    outpath::Union{Nothing,AbstractString}=nothing,
    title::Union{Nothing,AbstractString}=nothing,
    theme::Symbol=:light,
    bar_width::Union{Nothing,Real}=nothing,
    bargap::Union{Nothing,Real}=nothing,
    bargroupgap::Union{Nothing,Real}=nothing,
)
    header = nothing
    header_idx = Dict{String,Int}()
    rows = NamedTuple[]

    open(path, "r") do io
        for (ln, line) in enumerate(eachline(io))
            isempty(strip(line)) && continue
            if header === nothing
                header = _csv_split(strip(line))
                header_idx = Dict{String,Int}(h => i for (i, h) in enumerate(header))
                continue
            end

            cols = _csv_split(strip(line))
            length(cols) >= 2 || error("malformed CSV row at line $ln in $path: $line")

            algo_col = get(header_idx, "algo", nothing)
            algo_col === nothing && error("CSV missing required column `algo`: $path")

            algo = _strip_quotes(cols[algo_col])
            mean_u_norm = if haskey(header_idx, "mean_u_norm")
                parse(Float64, cols[header_idx["mean_u_norm"]])
            elseif haskey(header_idx, "u_norm")
                parse(Float64, cols[header_idx["u_norm"]])
            else
                error("CSV must contain `mean_u_norm` (summary) or `u_norm` (per-case): $path")
            end
            overall_success_rate = if haskey(header_idx, "overall_success_rate")
                parse(Float64, cols[header_idx["overall_success_rate"]])
            elseif haskey(header_idx, "success_rate")
                parse(Float64, cols[header_idx["success_rate"]])
            else
                error("CSV must contain `overall_success_rate` (summary) or `success_rate` (per-case): $path")
            end
            family, pop = _parse_tpflash_algo_population(algo)

            push!(rows, (; family, pop, mean_u_norm, overall_success_rate))
        end
    end

    isempty(rows) && error("no data rows found in $path")

    # If the input is a per-case CSV, aggregate to a single row per algo.
    if header !== nothing && haskey(header_idx, "u_norm") && haskey(header_idx, "success_rate") && !(haskey(header_idx, "mean_u_norm") && haskey(header_idx, "overall_success_rate"))
        acc = Dict{Tuple{String,Int},Tuple{Float64,Float64,Int}}()
        for r in rows
            k = (r.family, r.pop)
            su, ssr, n = get(acc, k, (0.0, 0.0, 0))
            acc[k] = (su + r.mean_u_norm, ssr + r.overall_success_rate, n + 1)
        end
        rows = [(; family=k[1], pop=k[2], mean_u_norm=su / n, overall_success_rate=ssr / n) for (k, (su, ssr, n)) in acc]
    end

    families = unique(getfield.(rows, :family))
    sort!(families)

    plot_title = title === nothing ? basename(path) : String(title)
    template = Config(
        font=Config(family="Arial", size=18, color=theme == :light ? "black" : "#dbdbdb"),
        plot_bgcolor=theme == :light ? "white" : "rgba(0,0,0,0)",
        paper_bgcolor=theme == :light ? "white" : "rgba(0,0,0,0)",
        margin=Config(l=70, r=15, t=40, b=60, pad=0),
        xaxis=Config(
            showline=true,
            linecolor=theme == :light ? "black" : "#dbdbdb",
            linewidth=1.5,
            mirror=true,
            ticks="inside",
            gridcolor=theme == :light ? "rgba(0,0,0,0.1)" : "rgba(255,255,255,0.2)",
            zeroline=false,
            standoff=9,
        ),
        legend=Config(
            bgcolor=theme == :light ? "white" : "rgba(0,0,0,0.5)",
            x=0.01,
            y=0.99,
            standoff=9,
        ),
    )
    template.yaxis = template.xaxis

    layout = deepcopy(template)
    layout.title = Config(text=plot_title)
    layout.barmode = "group"

    layout.xaxis = deepcopy(template.xaxis)
    layout.xaxis.title = Config(text="Population size")
    layout.xaxis.automargin = true
    layout.xaxis.tickmode = "linear"
    layout.xaxis.tick0 = 0
    layout.xaxis.dtick = 10

    layout.yaxis = deepcopy(template.xaxis)
    layout.yaxis.title = Config(text="Mean normalized U-score")
    layout.yaxis.automargin = true
    layout.yaxis.range = (0, 1)

    layout.yaxis2 = deepcopy(template.xaxis)
    layout.yaxis2.title = Config(text="Overall success rate")
    layout.yaxis2.overlaying = "y"
    layout.yaxis2.side = "right"
    layout.yaxis2.range = (0, 100)
    layout.yaxis2.ticksuffix = "%"
    layout.yaxis2.automargin = true

    bargap === nothing || (layout.bargap = Float64(bargap))
    bargroupgap === nothing || (layout.bargroupgap = Float64(bargroupgap))

    fig = nothing
    for fam in families
        fam_rows = [r for r in rows if r.family == fam]
        sort!(fam_rows; by = r -> r.pop)

        pops = getfield.(fam_rows, :pop)
        u = getfield.(fam_rows, :mean_u_norm)
        sr = 100 .* getfield.(fam_rows, :overall_success_rate)

        width = bar_width === nothing ? nothing : Float64(bar_width)
        if fig === nothing
            fig = plot.bar(; x=pops, y=u, name="$(fam) U-score", yaxis="y", layout=layout, width=width)
        else
            fig.bar(; x=pops, y=u, name="$(fam) U-score", yaxis="y", width=width)
        end
        fig.scatter(; x=pops, y=sr, mode="lines+markers", name="$(fam) success rate", yaxis="y2")
    end

    fig === nothing && error("no traces produced for $path")

    display(fig)
    if outpath !== nothing
        PlotlyKaleido.start()
        PlotlyKaleido.savefig(fig, String(outpath))
    end

    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    if isempty(ARGS)
        error("usage: julia plot_tpflash_summary.jl <path/to/*_summary.csv> [out.png]")
    end
    path = ARGS[1]
    outpath = length(ARGS) >= 2 ? ARGS[2] : nothing
    plot_tpflash_summary_csv(path; outpath=outpath)
end
