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
                amp = 22),
  band_voigtian(band_id = "band_2",
                x0 = 1130,
                sigma = 8,
                gamma = 12,
                amp = 25)
  # band_voigtian(band_id = "band_3",
  #               x0 = 1170,
  #               sigma = 1,
  #               gamma = 10,
  #               amp = 3)
)

band_preview(tmp,spectrum = "CoHA100",params = params)
