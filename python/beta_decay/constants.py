"""Physical constants for beta decay calculations."""

# Fine structure constant
ALPHA = 1.0 / 137.035999084

# Electron mass energy [MeV]
ME_MEV = 0.51099895

# Free nucleon axial coupling
G_A = 1.2756
G_V = 1.0

# Ratio g_A / g_V
G_A_OVER_V = G_A / G_V

# Weak decay constant D [s]
# t_{1/2} * f_0 = D / (g_A^2 * B_GT)  for pure GT allowed transition
D_CONST = 6289.0  # seconds

# kappa = D * g_A^2  [s]
KAPPA = D_CONST * G_A ** 2
