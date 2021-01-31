using DifferentialEquations
using Plots; plotly()
using DifferentialEquations.EnsembleAnalysis
using CSV
using DataFrames
#using DiffEqGPU
#using CuArrays

beta_p = 0.8
r_1 = 0.6
r_2 = 0.6
b = 0.5
beta_v = 0.8
theta = 0.818469652
mu = 0.058596601
gamma = 0.06
sigma_L = 0.378173863
sigma_I = 0.030015876
sigma_v = 0.254146903
N_v = mu/gamma
r= max(r_1,r_2)
u_0 = [97.0,1.0,2.0,3.0,4.0,100,3.0,4.0]
T = 1000.0
time = (0.0,T)
N_p = u_0[1]+u_0[2]+u_0[3]
dt=0.001
#=
p = Dict("beta_p"=>beta_p,"r_1"=>r,"r_2"=>r_2,"b"=>b,"beta_v"=>beta_v,"theta"=>theta,"mu"=>mu,
"gamma"=>gamma,"sigma_L"=>sigma_L,"sigma_I"=>sigma_I,"sigma_v"=>sigma_v,"N_v"=>N_v,"N_p"=>N_p)
=#
#=p[1]= beta_p
p[2]= r
p[3]= b
p[4]= beta_v
p[5]= theta
p[6]= mu
p[7]= gamma
p[8]= sigma_p
p[9]= sigma_v
p[10]= N_v
p[11]= N_p
=#
Rs0 = beta_p*beta_v/(r*gamma)
print("Rs0=",Rs0)

function F_Drift(du,u,p,t)
  du[1] = -beta_p*u[1]*u[5]/N_v+r_1*u[2]+r_2*u[3]
  du[2] = beta_p*u[1]*u[5]/N_v-b*u[2]-r_1*u[2]
  du[3] = b*u[2]-r_2*u[3]
  du[4] = -beta_v*u[4]*u[3]/N_p-gamma*u[4]+(1-theta)*mu
  du[5] = beta_v*u[4]*u[3]/N_p-gamma*u[5]+theta*mu
  du[6] =  -beta_p*u[1]*u[5]/N_v+r_1*u[2]+r_2*u[3]+beta_p*u[1]*u[5]/N_v-b*u[2]-r_1*u[2]+b*u[2]-r_2*u[3]
  du[7] = -beta_v*u[4]*u[3]/N_p-gamma*u[4]+(1-theta)*mu+ beta_v*u[4]*u[3]/N_p-gamma*u[5]+theta*mu
end

function G_Diffusion(du,u,p,t)
  du[1] = sigma_L*u[2]*u[1]/N_p+sigma_I*u[1]*u[3]/N_p
  du[2] = -sigma_L*u[1]*u[2]/N_p
  du[3] = -sigma_I*u[1]*u[3]/N_p
  du[4] = -sigma_v*u[4]
  du[5] = -sigma_v*u[5]
  du[6] = sigma_L*u[2]*u[1]/N_p+sigma_I*u[1]*u[3]/N_p-sigma_L*u[1]*u[2]/N_p-sigma_I*u[1]*u[3]/N_p
  du[7] = -sigma_v*u[4]-sigma_v*u[5]
end

prob_sde_tomato_sys = SDEProblem(F_Drift,G_Diffusion,u_0,time)
sol = solve(prob_sde_tomato_sys)
title = plot(title = "R_s =$Rs0", grid = false, showaxis = false, bottom_margin = -50Plots.px)
p1=plot(sol,vars=(1),color="darkgreen",title="Susc. p.")
p2=plot(sol,vars=(2),color="darkorange", title ="Lat. p.")
p3=plot(sol,vars=(3),color="darkred",title = "Infec. p.")
p4=plot(sol,vars=(4),color="green", title = "Susc. v.")
p5=plot(sol,vars=(5),color="red",title ="Infec. v.")

p6=plot(sol,vars=(6),color="red",title ="Law C. P.")
p7=plot(sol,vars=(7),color="red",title ="Law C. v.")


plot(p1,p2,p3,p4,p5,title,layout = @layout([ [A B C]; [D E F]]), label="")
ensembleprob = EnsembleProblem(prob_sde_tomato_sys)
#Sol = solve(ensembleprob,EnsembleGPUArray(),trajectories=100)
Sol = solve(ensembleprob,p=p,EnsembleThreads(),trajectories=2000)

summ_1 = EnsembleSummary(Sol,0:0.1:T) #Media
plotly()
plot(summ_1,idxs=(3),labels="Middle 95%")
summ_2 = EnsembleSummary(Sol,0:0.1:T;quantiles=[0.25,0.75])
plot!(summ_2,idxs=(3),labels="Middle 50%",legend=true)
#summ_3 = EnsembleSummary(Sol,0:0.1:1000;quantiles=[0.75])#,0.75])
#plot!(summ_2,idxs=(1),labels="Middle 50%",legend=true)
#plot!(summ_3,idxs=(1),labels="Middle 50%",legend=true)
################################################################################
########################    data  Media    #####################################
################################################################################
time = summ_1.t
xu = summ_1.u
xu_glued = hcat(xu...)
Xu1 = xu_glued[1:5:end]
Xu2 = xu_glued[2:5:end]
Xu3 = xu_glued[3:5:end]
Xu4 = xu_glued[4:5:end]
Xu5 = xu_glued[5:5:end]


DF1 = DataFrame(t = time,S_p = Xu1,L_p =Xu2,I_p = Xu3, S_v = Xu4, I_v = Xu5)
DF1_red = DF1[1:10:end,1:end]
#CSV.write("/home/gabrielsalcedo/Dropbox/Artículos/JuliaPro_code/R0+1//DataSolutionR0greather1.csv",DF1_red)
CSV.write("C:/Users/adria/Dropbox/Artículos/JuliaPro_code/R0+1//DataSolutionR0greather1.csv",DF1_red)


################################################################################
########################   quartile data  Media    #############################
################################################################################
time = summ_1.t
xqlow = summ_1.qlow.u
X_glued = hcat(xqlow...)
X1 = X_glued[1:5:end]
X2 = X_glued[2:5:end]
X3 = X_glued[3:5:end]
X4 = X_glued[4:5:end]
X5 = X_glued[5:5:end]

yqhigh = summ_1.qhigh.u
Y_glued = hcat(yqhigh...)
Y1 = Y_glued[1:5:end]
Y2 = Y_glued[2:5:end]
Y3 = Y_glued[3:5:end]
Y4 = Y_glued[4:5:end]
Y5 = Y_glued[5:5:end]

DF3 = DataFrame(t = time,S_p = X1,L_p =X2,I_p = X3, S_v = X4, I_v = X5)
DF3_red = DF3[1:10:end,1:end]
DF4 = DataFrame(t = time,S_p = Y1,L_p =Y2,I_p = Y3, S_v = Y4, I_v = Y5)
DF4_red = DF4[1:10:end,1:end]
#CSV.write("/home/gabrielsalcedo/Dropbox/Artículos/JuliaPro_code/R0+1//DataQLowR0greather1.csv",DF3_red)
#CSV.write("/home/gabrielsalcedo/Dropbox/Artículos/JuliaPro_code/R0+1//DataQHighR0greather1.csv",DF4_red)
CSV.write("C:/Users/adria/Dropbox/Artículos/JuliaPro_code/R0+1//DataQLowR0greather1.csv",DF3_red)
CSV.write("C:/Users/adria/Dropbox/Artículos/JuliaPro_code/R0+1//DataQHighR0greather1.csv",DF4_red)



################################################################################
##########################   Data quantiles    #################################
################################################################################
time = summ_2.t
xxqlow = summ_2.qlow.u
XX_glued = hcat(xxqlow...)
XX1 = XX_glued[1:5:end]
XX2 = XX_glued[2:5:end]
XX3 = XX_glued[3:5:end]
XX4 = XX_glued[4:5:end]
XX5 = XX_glued[5:5:end]

yyqhigh = summ_2.qhigh.u
YY_glued = hcat(yyqhigh...)
YY1 = YY_glued[1:5:end]
YY2 = YY_glued[2:5:end]
YY3 = YY_glued[3:5:end]
YY4 = YY_glued[4:5:end]
YY5 = YY_glued[5:5:end]

DF_qantile_low = DataFrame(t = time,S_p = XX1,L_p =XX2,I_p = XX3, S_v = XX4,
 I_v = XX5)
DF_qantile_low_red = DF_qantile_low[1:10:end,1:end]

DF_qantile_high = DataFrame(t = time,S_p = YY1,L_p =YY2,I_p = YY3, S_v = YY4,
 I_v = YY5)
DF_qantile_high_red = DF_qantile_high[1:10:end,1:end]

#CSV.write("/home/gabrielsalcedo/Dropbox/Artículos/JuliaPro_code/R0+1//Data_qlowR0greather1.csv",
#DF_qantile_low_red)
#CSV.write("/home/gabrielsalcedo/Dropbox/Artículos/JuliaPro_code/R0+1//Data_qhighR0greather1.csv",
#DF_qantile_high_red)
CSV.write("C:/Users/adria/Dropbox/Artículos/JuliaPro_code/R0+1//Data_qlowR0greather1.csv",DF_qantile_low_red)
CSV.write("C:/Users/adria/Dropbox/Artículos/JuliaPro_code/R0+1//Data_qhighR0greather1.csv",DF_qantile_high_red)
###############################################################################
DF1 = DataFrame(time = summ_1.t,x = summ_1.u[1,:],y =summ_1.u[2,:], z =
summ_1.u[3,:], u = summ_1.u[4,:], w = summ_1.u[5,:])

DF2 = DataFrame(time = summ_1.t,x = summ_1.v[1,:],y =summ_1.v[2,:], z =
summ_1.v[3,:], u = summ_1.v[4,:], w = summ_1.v[5,:])

DF3 = DataFrame(t = time,S_p = X1,L_p =X2,I_p = X3, S_v = X4, I_v = X5)
DF3_red = DF3[1:10:end,1:end]
DF4 = DataFrame(t = time,S_p = Y1,L_p =Y2,I_p = Y3, S_v = Y4, I_v = Y5)
DF4_red = DF4[1:10:end,1:end]
DFQL = DataFrame(ql)
CSV.write("/home/gabrielsalcedo/Dropbox/Artículos/JuliaPro_code//DataU.csv",DF1)
CSV.write("/home/gabrielsalcedo/Dropbox/Artículos/JuliaPro_code//DataV.csv",DF2)


#CSV.write("/home/gabrielsalcedo/Dropbox/Artículos/JuliaPro_code//QlData.csv",DFQL)
