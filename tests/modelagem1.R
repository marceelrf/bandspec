library(tidyverse)
devtools::load_all("../tidyspec")

set_spec_wn(col_name = "Wavenumber")

tmp <- spec_select(.data = Region, CoHA100)

spec_smartplot(tmp)
spec_smartplotly(tmp)

params <- band_params(
  band_voigtian(band_id = "band_1",
                x0 = 1070,
                sigma = 10,
                gamma = 2,
                amp = 25),
  band_voigtian(band_id = "band_2",
                x0 = 1128,
                sigma = 4,
                gamma = 15,
                amp = 30),
  band_voigtian(band_id = "band_3",
                x0 = 1168,
                sigma = 1,
                gamma = 10,
                amp = 2),
  band_voigtian(band_id = "band_4",
                x0 = 1234,
                sigma = 10,
                gamma = 10,
                amp = 2)
)

band_preview(tmp,spectrum = "CoHA100",params = params)


band_fit(.data = tmp,
         spectrum = "CoHA100",
         params = params,
         control = nls.lm.control(maxiter = 1024))
