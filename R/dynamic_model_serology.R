# - - - - - - - - - - - - - - - - - - - - - - - 
# Define serological data
# Author: Adam Kucharski
# github.com/adamkucharski/fiji-denv3-2014
# - - - - - - - - - - - - - - - - - - - - - - - 


# - - - - - - - - 
# ELISA and MIA results - AGE GROUPS: 0<= x <20 and x>=20

n_ELISA_C_D = c(21,42,86)   # DENV-3 immunity ELISA
n_ELISA_A_D = c(133,154,176) # DENV-3 immunity ELISA

n_Luminex_C_D3 = c(6,36,85)   # DENV-3 MIA immunity
n_Luminex_A_D3 = c(81,104,176) # DENV-3 MIA immunity

# Choose whether to use
if(use.ELISA.data == T){n_Luminex_C_D3 = n_ELISA_C_D; n_Luminex_A_D3 = n_ELISA_A_D}
