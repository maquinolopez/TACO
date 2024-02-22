setwd('~/Documents/Taco/')
source('TACO.R')

folder_names <- basename(list.dirs())[-1]


for (i in 58:202){
  print(folder_names[i])
  tac <- try(Taco(folder_names[i] , '~/Documents/Taco/', th = 25 ,
                  burnin = 2e+3, s_size = 2.5e+3, newrun = F ), silent = TRUE)
  if(inherits(tac, "try-error")) {
    Taco(folder_names[i] , '~/Documents/Taco/', th = 25 ,
         burnin = 2e+3, s_size = 2.5e+3, newrun = T )
  }
  
}


                                        


 # tac <-Taco("Galvarne_Moraine_Ridge_Bjorck", '~/Documents/Taco/', th = 25 ,
            # burnin = 2e+3, s_size = 2.5e+3, newrun = T )

# 'Aero_M4' # Este da un erroir desde el inicio investigar

# Nota: revisar unidades

# 186
# "Tasiusaq"

# 199
# "Wylde Lake bog"

# 140
# "PAT-HB-2010"