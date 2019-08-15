using CSV, Gadfly, Cairo, Fontconfig, Dates

if ispath("/moto/sscc/projects/biasedexpectations")
  root_dir = "/moto/sscc/projects/biasedexpectations"
elseif ispath("/research/hmc")
  root_dir = "/research/hmc"
else
  @error "No valid directory for root directory found"
  exit(1)
end
cd(root_dir)

set_default_plot_size(18inch, 30inch)

t = Theme(plot_padding = [1inch, 2inch, 1inch, 1inch])

function colors(rng = :)
    return ["dark red", "dark blue", "dark orange", "purple", "black"][rng]
end

set_default_plot_size(18inch, 12inch)


function plotchain(date, theme = Theme())
  p = joinpath(root_dir, "data/output/official")
  df = CSV.read(joinpath(p, "filtered_means_$(date).csv"), copycols=true)
  fig = plot(
    layer(df[1:100_000,:], x = :state_3, Geom.density, Theme(default_color = colors(1)) ), 
    layer(df[1:150_000,:], x = :state_3, Geom.density,Theme(default_color = colors(2)) ), 
    layer(df[1:200_000,:], x = :state_3, Geom.density,Theme(default_color = colors(3)) ), 
    layer(df, x = :state_3, Geom.density, Theme(default_color = colors(5)) ), 
    theme,
    Guide.ylabel("Probability"), 
    Guide.title("Distribution of Parameter Draws"), 
    Guide.manual_color_key("Legend", "1:" .* ["100,000", "150,000", "200,000","250,000"], colors([1,2,3,5])),
    )
return fig
end

d = Date(2009,06)
plotchain(d, t) |> PDF("plots/chain_mean_$(d).pdf")


