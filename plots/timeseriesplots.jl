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


function plot_mean(df, noise)
  m1 = layer(df, x = :date, y = :state_1_mean, Geom.line, Theme(default_color = "red"))
  m2 = layer(df, x = :date, y = :state_2_mean, Geom.line, Theme(default_color = "blue"))
  m3 = layer(df, x = :date, y = :state_3_mean, Geom.line, Theme(default_color = "green"))
  set_default_plot_size(12inch, 10inch)
  fig = plot(m1, m2, m3, Guide.xlabel("Sample End Date"), Guide.ylabel("State Mean"))
  fig |> PDF("plots/state_mean_$(noise).pdf")
end

p = joinpath(root_dir, "data/output/signals_official_noise_1.0_allsignal")
df1 = CSV.read(joinpath(p, "filtered_means_dispersion.csv"))
plot_mean(df1, "1.0")

p = joinpath(root_dir, "data/output/signals_official_noise_0.1_allsignal")
df2 = CSV.read(joinpath(p, "filtered_means_dispersion.csv"))
plot_mean(df2, "0.1")

p = joinpath(root_dir, "data/output/official")
df3 = CSV.read(joinpath(p, "filtered_means_summary.csv"))
plot_mean(df3, "0.0")

p = joinpath(root_dir, "data/output/official")
df4 = CSV.read(joinpath(p, "filtered_variances_summary.csv"))


function plot_var(df, noise)
  m1 = layer(df, x = :date, y = :state_1_mean, Geom.line, Theme(default_color = "red"))
  m2 = layer(df, x = :date, y = :state_2_mean, Geom.line, Theme(default_color = "blue"))
  m3 = layer(df, x = :date, y = :state_3_mean, Geom.line, Theme(default_color = "green"))
  set_default_plot_size(12inch, 10inch)
  fig = plot(m1, m2, m3, Guide.xlabel("Sample End Date"), Guide.ylabel("State Std."),  Guide.manual_color_key("Legend", ["State 1", "State 2", "State 3"], ["red", "blue", "green"]) )
  fig |> PDF("plots/state_std_$(noise).pdf")
end

p = joinpath(root_dir, "data/output/official")
df4 = CSV.read(joinpath(p, "filtered_variances_summary.csv"), copycols=true)
df4[:,[:state_1_mean, :state_2_mean, :state_3_mean]] .= sqrt.(df4[:,[:state_1_mean, :state_2_mean, :state_3_mean]])
plot_var(df4, "0.0")


function plot_forecast(df, noise)
  m1 = layer(df, x = :date, y = :forecast_12_mean, Geom.line, Theme(default_color = "red"))
  set_default_plot_size(12inch, 10inch)
  fig = plot(m1, Guide.xlabel("Sample End Date"), Guide.ylabel("12-Month Forecast") )
  fig |> PDF("plots/state_forecast_$(noise).pdf")
end

p = joinpath(root_dir, "data/output/official")
df5 = CSV.read(joinpath(p, "forecasts_summary.csv"), copycols=true)
plot_forecast(df5, "0.0")


function plot_forecast_state(df1, df2, noise)
  m1 = layer(df1, x = :date, y = :forecast_12_mean, Geom.line, Theme(default_color = "red"))
  m2 = layer(df2, x = :date, y = :state_1_mean, Geom.line, Theme(default_color = "red"))
  m3 = layer(df2, x = :date, y = :state_2_mean, Geom.line, Theme(default_color = "blue"))
  m4 = layer(df2, x = :date, y = :state_3_mean, Geom.line, Theme(default_color = "green"))
  set_default_plot_size(12inch, 10inch)
  fig1 = plot(m1, Guide.xlabel("Sample End Date"), Guide.manual_color_key("Ignore", ["State 1", "State 2", "State 3"], ["red", "blue", "green"]), Coord.cartesian(xmin=Date(1980,01) ), Guide.title("12-Month Forecasts") )
  fig2 = plot(m2, m3, m4, Guide.xlabel("Sample End Date"), Guide.ylabel("State Probabilities"),  Guide.manual_color_key("Legend", ["State 1", "State 2", "State 3"], ["red", "blue", "green"]), Coord.cartesian(xmin=Date(1980,01) ), Guide.title("State Probabilities at Time of Forecast") )
  fig = vstack(fig1, fig2)
  fig |> PDF("plots/state_forecast_probs_$(noise).pdf")
end

p = joinpath(root_dir, "data/output/official")
df5 = CSV.read(joinpath(p, "forecasts_summary.csv"), copycols=true)
df6 = CSV.read(joinpath(p, "filtered_state_probs_summary.csv"), copycols=true)
plot_forecast_state(df5, df6, "0.0")

function plot_mean(d1, d2, d3, d4; root_dir = root_dir)
  p = joinpath(root_dir, "data/output/official")
  df1 = CSV.read(joinpath(p, "filtered_means_$(d1).csv"))
  df2 = CSV.read(joinpath(p, "filtered_means_$(d2).csv"))
  df3 = CSV.read(joinpath(p, "filtered_means_$(d3).csv"))
  df4 = CSV.read(joinpath(p, "filtered_means_$(d4).csv"))
  m1 = layer(df1, x = :state_1, Geom.density, Theme(default_color = "red"))
  m2 = layer(df2, x = :state_1, Geom.density, Theme(default_color = "blue"))
  m3 = layer(df3, x = :state_1, Geom.density, Theme(default_color = "green"))
  m4 = layer(df4, x = :state_1, Geom.density, Theme(default_color = "purple"))
  set_default_plot_size(12inch, 10inch)
  fig1 = plot(m1, m2, m3, m4, Guide.xlabel("State 1 Mean"), Guide.ylabel("Density"), Guide.manual_color_key("Legend", ["$(d1)", "$(d2)", "$(d3)", "$(d4)"], ["red", "blue", "green", "purple"]))
  m1 = layer(df1, x = :state_2, Geom.density, Theme(default_color = "red"))
  m2 = layer(df2, x = :state_2, Geom.density, Theme(default_color = "blue"))
  m3 = layer(df3, x = :state_2, Geom.density, Theme(default_color = "green"))
  m4 = layer(df4, x = :state_2, Geom.density, Theme(default_color = "purple"))
  set_default_plot_size(12inch, 10inch)
  fig2 = plot(m1, m2, m3, m4, Guide.xlabel("State 2 Mean"), Guide.ylabel("Desnsity"), Guide.manual_color_key("Legend", ["$(d1)", "$(d2)", "$(d3)", "$(d4)"], ["red", "blue", "green", "purple"]))
  m1 = layer(df1, x = :state_3, Geom.density, Theme(default_color = "red"))
  m2 = layer(df2, x = :state_3, Geom.density, Theme(default_color = "blue"))
  m3 = layer(df3, x = :state_3, Geom.density, Theme(default_color = "green"))
  m4 = layer(df4, x = :state_3, Geom.density, Theme(default_color = "purple"))
  set_default_plot_size(12inch, 10inch)
  fig3 = plot(m1, m2, m3, m4, Guide.xlabel("State 3 Mean"), Guide.ylabel("Density"), Guide.manual_color_key("Legend", ["$(d1)", "$(d2)", "$(d3)", "$(d4)"], ["red", "blue", "green", "purple"]))
  fig = vstack(fig1, fig2, fig3)
  fig |> PDF("plots/state_mean_$(d1)_$(d2)_$(d3).pdf")
end

plot_mean(Date(2009,4), Date(2009,5), Date(2009,6), Date(2009,7))


function plot_var(d1, d2, d3, d4; root_dir = root_dir)
  p = joinpath(root_dir, "data/output/official")
  df1 = CSV.read(joinpath(p, "filtered_variances_$(d1).csv"))
  df2 = CSV.read(joinpath(p, "filtered_variances_$(d2).csv"))
  df3 = CSV.read(joinpath(p, "filtered_variances_$(d3).csv"))
  df4 = CSV.read(joinpath(p, "filtered_variances_$(d4).csv"))
  m1 = layer(df1, x = :state_1, Geom.density, Theme(default_color = "red"))
  m2 = layer(df2, x = :state_1, Geom.density, Theme(default_color = "blue"))
  m3 = layer(df3, x = :state_1, Geom.density, Theme(default_color = "green"))
  m4 = layer(df4, x = :state_1, Geom.density, Theme(default_color = "purple"))
  set_default_plot_size(12inch, 10inch)
  fig1 = plot(m1, m2, m3, m4, Guide.xlabel("State 1 Variance"), Guide.ylabel("Density"), Guide.manual_color_key("Legend", ["$(d1)", "$(d2)", "$(d3)", "$(d4)"], ["red", "blue", "green", "purple"]))
  m1 = layer(df1, x = :state_2, Geom.density, Theme(default_color = "red"))
  m2 = layer(df2, x = :state_2, Geom.density, Theme(default_color = "blue"))
  m3 = layer(df3, x = :state_2, Geom.density, Theme(default_color = "green"))
  m4 = layer(df4, x = :state_2, Geom.density, Theme(default_color = "purple"))
  set_default_plot_size(12inch, 10inch)
  fig2 = plot(m1, m2, m3, m4, Guide.xlabel("State 2 Variance"), Guide.ylabel("Density"), Guide.manual_color_key("Legend", ["$(d1)", "$(d2)", "$(d3)", "$(d4)"], ["red", "blue", "green", "purple"]))
  m1 = layer(df1, x = :state_3, Geom.density, Theme(default_color = "red"))
  m2 = layer(df2, x = :state_3, Geom.density, Theme(default_color = "blue"))
  m3 = layer(df3, x = :state_3, Geom.density, Theme(default_color = "green"))
  m4 = layer(df4, x = :state_3, Geom.density, Theme(default_color = "purple"))
  set_default_plot_size(12inch, 10inch)
  fig3 = plot(m1, m2, m3, m4, Guide.xlabel("State 3 Variance"), Guide.ylabel("Density"), Guide.manual_color_key("Legend", ["$(d1)", "$(d2)", "$(d3)", "$(d4)"], ["red", "blue", "green", "purple"]))
  fig = vstack(fig1, fig2, fig3)
  fig |> PDF("plots/state_var_$(d1)_$(d2)_$(d3).pdf")
end

plot_var(Date(2009,4), Date(2009,5), Date(2009,6), Date(2009,7))