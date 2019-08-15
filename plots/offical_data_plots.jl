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

t = Theme(
  plot_padding = [1inch, 2inch, 1inch, 1inch]
)

xticks = range(DateTime(1980,1),DateTime(2019,1),step=Dates.Year(5))

function colors(rng=:)
  return ["dark red", "dark blue", "dark orange", "black"][rng]
end

function plot_mean(df)
  fig = plot(df, 
    layer(x = :date, y = :state_1_mean, Geom.line, Theme(default_color = colors(1)) ),
    layer(x = :date, y = :state_2_mean, Geom.line, Theme(default_color = colors(2)) ),
    layer(x = :date, y = :state_3_mean, Geom.line, Theme(default_color = colors(3)) ),
    layer(x = :date, y = :offic_inf, Geom.line, Theme(default_color = colors(4)) ),
    Guide.xlabel("Sample End Date"), 
    Guide.ylabel("State Mean"), 
    Guide.title("Means of Distributions"), 
    Coord.cartesian(xmin=Date(1980,01), xmax=Date(2018,03)),
    Scale.x_continuous(minvalue=Date(1980,01), maxvalue=Date(2018,03)),
    Guide.xticks(ticks=xticks),
    Guide.manual_color_key("Legend", append!(["State $(i)" for i in 1:3], ["Inflation "]), colors()),
    t
    )
  return fig
end

function plot_std(df)
  fig = plot(df, 
    layer(x = :date, y = :state_1_mean, Geom.line, Theme(default_color = colors(1)) ),
    layer(x = :date, y = :state_2_mean, Geom.line, Theme(default_color = colors(2)) ),
    layer(x = :date, y = :state_3_mean, Geom.line, Theme(default_color = colors(3)) ),
    Guide.xlabel("Sample End Date"), 
    Guide.ylabel("State Std."), 
    Guide.title("Std. of Distributions"), 
    Coord.cartesian(xmin=Date(1980,01), xmax=Date(2018,03)),
    Guide.xticks(ticks=xticks),
    Guide.manual_color_key("Legend", ["State $(i)" for i in 1:3], colors()),
    t
    )
  return fig
end

function plot_forecast(df)
  fig = plot(df, 
    layer(x = :date, y = :forecast_12_mean, Geom.line, Theme(default_color = colors(1)) ),
    layer(x = :date, y = :offic_inf, Geom.line, Theme(default_color = colors(4)) ),
    Guide.xlabel("Sample End Date"), 
    Guide.ylabel("Forecast"), 
    Guide.title("12-Month Ahead Forecasts"), 
    Coord.cartesian(xmin=Date(1980,01), xmax=Date(2018,03)),
    Guide.xticks(ticks=xticks),
    Guide.manual_color_key("Legend", ["Forecast", "Inflation"], colors()),
    t
    )
end

function plot_state_probs(df)
  fig = plot(df, 
    layer(x = :date, y = :state_1_mean, Geom.step, Theme(default_color = colors(1)) ),
    layer(x = :date, y = :state_2_mean, Geom.step, Theme(default_color = colors(2)) ),
    layer(x = :date, y = :state_3_mean, Geom.step, Theme(default_color = colors(3)) ),
    Guide.xlabel("Sample End Date"), 
    Guide.ylabel("State Probabilities"), 
    Guide.title("State Probabilities at Sample End Date"), 
    Coord.cartesian(xmin=Date(1980,01), xmax=Date(2018,03)),
    Guide.xticks(ticks=xticks),
    Guide.manual_color_key("Legend", ["State $(i)" for i in 1:3], colors()),
    t
    )
  return fig
end

p = joinpath(root_dir, "data/output/official")
df1 = CSV.read(joinpath(p, "filtered_means_summary.csv"), copycols=true)
df2 = CSV.read(joinpath(p, "filtered_variances_summary.csv"), copycols=true)
df3 = CSV.read(joinpath(p, "filtered_state_probs_summary.csv"), copycols=true)
df4 = CSV.read(joinpath(p, "forecasts_summary.csv"), copycols=true)
df5 = CSV.read(joinpath(root_dir, "data/processed/inflation.csv"))

df1 = join(df1, df5, on = :date, kind = :left)
df2 = join(df2, df5, on = :date, kind = :left)
df3 = join(df3, df5, on = :date, kind = :left)
df4 = join(df4, df5, on = :date, kind = :left)

df2[:,[:state_1_mean, :state_2_mean, :state_3_mean]] .= sqrt.(df2[:,[:state_1_mean, :state_2_mean, :state_3_mean]])

vstack(plot_mean(df1), plot_std(df2), plot_forecast(df4), plot_state_probs(df3)) |> PDF("plots/official_data.pdf")


t = Theme(
  plot_padding = [1inch, 2inch, 0.2inch, 0.2inch]
)
set_default_plot_size(12inch, 20inch)
function plot_mean(d1, d2, d3, d4; root_dir = root_dir)
  p = joinpath(root_dir, "data/output/official")
  df1 = CSV.read(joinpath(p, "filtered_means_$(d1).csv"))
  df2 = CSV.read(joinpath(p, "filtered_means_$(d2).csv"))
  df3 = CSV.read(joinpath(p, "filtered_means_$(d3).csv"))
  df4 = CSV.read(joinpath(p, "filtered_means_$(d4).csv"))
  m1 = layer(df1, x = :state_1, Geom.density, Theme(default_color = colors(1)))
  m2 = layer(df2, x = :state_1, Geom.density, Theme(default_color = colors(2)))
  m3 = layer(df3, x = :state_1, Geom.density, Theme(default_color = colors(3)))
  m4 = layer(df4, x = :state_1, Geom.density, Theme(default_color = colors(4)))
  fig1 = plot(m1, m2, m3, m4, Guide.xlabel("Inflation Rate"), Guide.ylabel("Density"), Guide.manual_color_key("Legend", ["$(d1)", "$(d2)", "$(d3)", "$(d4)"], colors()),t, Guide.title("Distribution of μ, Low Inflation State"))
  m1 = layer(df1, x = :state_2, Geom.density, Theme(default_color = "red"))
  m2 = layer(df2, x = :state_2, Geom.density, Theme(default_color = "blue"))
  m3 = layer(df3, x = :state_2, Geom.density, Theme(default_color = "green"))
  m4 = layer(df4, x = :state_2, Geom.density, Theme(default_color = "purple"))
  fig2 = plot(m1, m2, m3, m4, Guide.xlabel("Inflation Rate"), Guide.ylabel("Desnsity"), Guide.manual_color_key("Legend", ["$(d1)", "$(d2)", "$(d3)", "$(d4)"], colors()),t, Guide.title("Distribution of μ, Mid-Inflation State"))
  m1 = layer(df1, x = :state_3, Geom.density, Theme(default_color = "red"))
  m2 = layer(df2, x = :state_3, Geom.density, Theme(default_color = "blue"))
  m3 = layer(df3, x = :state_3, Geom.density, Theme(default_color = "green"))
  m4 = layer(df4, x = :state_3, Geom.density, Theme(default_color = "purple"))
  fig3 = plot(m1, m2, m3, m4, Guide.xlabel("Inflation Rate"), Guide.ylabel("Density"), Guide.manual_color_key("Legend", ["$(d1)", "$(d2)", "$(d3)", "$(d4)"], colors()),t, Guide.title("Distribution of μ, High Inflation State"))
  fig = vstack(fig1, fig2, fig3)
  fig |> PDF("plots/state_mean_distributions.pdf")
end


plot_mean(Date(2009,4), Date(2009,5), Date(2009,6), Date(2009,7))