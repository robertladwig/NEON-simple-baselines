installed_packages <- installed.packages()
if (is.element('feasts', installed_packages) == F) {
  install.packages('feasts')
}

# Script to run forecasts
source('./Models/Physics.R')
message('Physics model submitted')
