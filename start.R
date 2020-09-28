source('dependencies.R')

runApp('main', port = 18232, launch.browser = F)

# ssh -L 12345:localhost:18232 jsicherman@nelson.msl.ubc.ca -p 22000
