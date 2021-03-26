#This code is for automatically inputting chemical formulas into ode
#This code is highly based on Dr. Alex's code:https://github.com/ata27/CH4_ISA_Box_Model

#set constant we need:
assign('M',2.50e19)#air conc



require(deSolve)

#reading chemical formula data from excel file, automatically produce ODE output
require(readxl)
input<- data.frame(read_excel("input.xlsx",sheet="input"))
statedata<-data.frame(read_excel("input.xlsx",sheet="state"))


# functions
convert_emissions <- function(mass, molar_mass) {
  # function to convert the input emissions
  # constants
  C_fac = 2.5E19 # convert to molecules cm-3
  sec_per_year = 86400*360
  #Vol_Earth = 4*pi*(6378E3)^2 * (5000) # volume = area * height - 5km
  top <- (mass/molar_mass) # the mass of the species converted into moles
  bottom <- (5.2E21/28.8) # mass of air = 5.2E18 Kg and rmm = 28.8 g/mol = moles of air
  temp <- top/bottom
  return(temp*C_fac/sec_per_year) # return the result in molecules cm-3 s-1
}

# What we are going to do is solve some simple, but stiff, ODE's. 
# OH + CH4 = CO+H2 k1
# Cl + CH4 = CO k2
# OH + CO = CO2 (don't care about CO2 here) k3
# OH + X = LOSS (don't care how) k4
#H2+OH=H+H2O k5 (don't care about H and H2O)
#H2+Cl=H+HCl k6 (don't care about HCl)(ignore this path first)
#H2+X=soil ptake k7

# Let's write d[X]/dt as dX


# set times for the main loop - in seconds so will need to extend a lot for this slow chemistry! 
n_years <- 100
t_step  <- 86400*30
times <- seq(0, (n_years*360*86400), t_step)

# generate a data frame of ramp data
# at a lower time interval to the main loop
times2 <- seq(0, (n_years*360*86400), t_step*4)

# generate a linear sequence to increase the emissions up to a maximum factor (max_fac)
max_fac <- 10
lin.seq <- seq(1, max_fac, length.out = length(times2))
emis_frame  <- data.frame(day=times2, emis_scale=lin.seq )
# now generate a function to interpolate the emis factors
emis.func <- approxfun(emis_frame)

#produce state:initial concentration of gases
#order of elements in state must match the order of ode
state<-c()
state_name<-c()
ode<-c()
for (i in 1:length(row.names(statedata))){
  #  state[i]<-assign(statedata[i,1],statedata[i,2])#produce state
  state[i]<-statedata[i,2]
  state_name[i]<-statedata[i,1]#produce name of state
  ode[i]<-paste("d",statedata[i,1],sep="",collapse = "")
}
names(state)<-state_name#give the name to the states, this is necessary

#produce parameter:reaction rate constant
parameters<-c()
parameters_name<-c()
for (i in 1:length(row.names(input))){
  #parameters[i]<-assign(input[i,9],input[i,10])#produce parameter
  parameters[i]<-input[i,10]
  parameters_name[i]<-input[i,9]
}
names(parameters)<-parameters_name#give the name to the parameters



kinetics <- function(t, state, parameters){
  
  with(as.list(c(state, parameters)), {
    
    #find reaction that contains the ode, and put the reaction into ode expression
    expre<-function(ode_name){
      minus<-c(input[which(input$reactant1==ode_name),"expression"],
               input[which(input$reactant2==ode_name),"expression"],
               input[which(input$reactant3==ode_name),"expression"],
               input[which(input$reactant4==ode_name),"expression"])
      temp_minus=0
      if(!identical(minus,character(0))){
        for(i in 1:length(minus)){
          temp_minus=temp_minus-eval(parse(text=minus[i]))
        }
      }
      
      plus<-c(input[which(input$product1==ode_name),"expression"],
              input[which(input$product2==ode_name),"expression"],
              input[which(input$product3==ode_name),"expression"],
              input[which(input$product4==ode_name),"expression"])
      temp_plus=0
      if(!identical(plus,character(0))){
        for(j in 1:length(plus)){
          temp_plus=temp_plus+eval(parse(text=plus[j]))
        }
      }
      
      return(temp_minus+temp_plus)
    }
    
    emis_factor <- emis.func(t)
    anth_emi<-c()
    anth_name<-c()
    for (i in 1:length(row.names(statedata))){
      anth_name[i]<-statedata[i,3]
      #emission increases, here we set emission of h2 increases in the next 100 years
      if(statedata[i,3]=="SH2"){
        anth_emi[i]=assign(statedata[i,3],statedata[i,6]*emis_factor)
        next
      }
      anth_emi[i]=assign(statedata[i,3],statedata[i,6])#create variables of anthropogenic emission
      

    }
    names(anth_emi)<-anth_name
    
    
    # define emission factor to be time dependent 
    # and set emission rates
    # rate of change
    real_ode<-c()
    for(i in 1:length(state)){
      real_ode[i]<-assign(ode[i],anth_emi[i]+expre(state_name[i]))#produce ode
    }
    names(real_ode)<-ode
    # return list of output
    list(real_ode, 
         emis_factor=emis_factor, 
         anth_emi)
  })
}
# solve!
out <- ode(y=state, times=times, func=kinetics, parms=parameters,method= "radau")
print('ODE is solved')
###########print(state)print(state)################################## now analyse the data

# # quick plot
print('quick plot')
plot(out, type="l")
#plot(out, type="l", lwd#=2, 
#     xlab="Time (s)")
#grid()

# the output of the ODE solver is a bit odd so will convert to a data frame for analysis
print('put output data into out_data')
out_data <- as.data.frame(out)