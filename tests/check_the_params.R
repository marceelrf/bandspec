
# Criar o grid de parâmetros
param_grid <- expand.grid(sigma = c(1, 5, 10),
                          gamma = c(1, 5, 10))
param_grid
# Gerar a lista de bandas
params_models <- do.call(band_params,
                         apply(param_grid, 1, function(row) {
                           band_voigtian(
                             band_id = paste0("band_", which(apply(param_grid, 1, function(r) all(r == row)))),
                             x0 = 1000,  # substitua pelo valor x0 desejado
                             sigma = row["sigma"],
                             gamma = row["gamma"],
                             amp = 1
                           )
                         })
)
viridis::inferno(n = 16)
band_preview(950:1050, params = params_models,show_sum = F,
             colors = viridis::inferno(n = 9))
