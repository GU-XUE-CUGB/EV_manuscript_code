# Load necessary library
library(npreg)

# Apply the ss function to the dataset, specifying the 'Midpoint' and 'IQR' columns, and spar parameter
mod.ss <- with(dataset, ss(Midpoint, IQR, spar = 0.2))

# View the result
mod.ss

# Plot the result
plot(mod.ss)
