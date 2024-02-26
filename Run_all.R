setwd('~/Documents/Taco/')
source('TACO.R')

folder_names <- basename(list.dirs())[-1]


# for (i in 14:202){
#   print(folder_names[i])
# 
#     Taco(folder_names[i] , '~/Documents/Taco/', th = 50 ,
#          burnin = 3.0e+3, s_size = 2.5e+3, newrun = T )
# 
# }



# 5
tac <- Taco(folder_names[12] , '~/Documents/Taco/', th = 50,
     burnin = 5.0e+3, s_size = 2.5e+3, newrun = T )
# 
# 
# tac$Bd_mat
# 
# tac <- Taco('Dead_Island_Swindles' , '~/Documents/Taco/', th = 50 ,
#                        burnin = 3.0e+3, s_size = 2.5e+3, newrun = T )


# tac <- Taco('Wylde Lake bog' , '~/Documents/Taco/', th = 50 ,
            # burnin = 3.0e+3, s_size = 2.5e+3, newrun = T )


# print(tac$Bd_mat)
 # tac <-Taco("Galvarne_Moraine_Ridge_Bjorck", '~/Documents/Taco/', th = 25 ,
            # burnin = 2e+3, s_size = 2.5e+3, newrun = T )

# 'Aero_M4' # Este da un erroir desde el inicio investigar

# '"JBL_2"' dows not iniciate 

# Nota: revisar unidades

# 186
# "Tasiusaq"

# 199
# "Wylde Lake bog"

# 140
# "PAT-HB-2010"