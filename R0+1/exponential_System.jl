using DifferentialEquations
using Plots; plotly()
using DifferentialEquations.EnsembleAnalysis
using CSV
using DataFrames
using StochasticDiffEq
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
u_0 = [log(97.0),log(1.0),log(2.0),log(3.0),log(4.0)]
T = 1000.0
time = (0.0,T)
N_p = exp(u_0[1])+exp(u_0[2])+exp(u_0[3])
dt=0.001
u0= [97.0/N_p,1.0/N_p,2.0/N_p,3.0/N_v,4.0/N_v]

Rs0 = beta_p*beta_v/(r*gamma)
print("Rs0=",Rs0)

function F_Drift(du,u,p,t)
  du[1] = -(beta_p/N_v)*exp(u[5])+r_1*(exp(u[2])/exp(u[1]))+r_2*(exp(u[3])/exp(u[1]))
  du[2] = (beta_p/N_v)*exp(u[1])*exp(u[5])/exp(u[2])-(b+r_1)
  du[3] = b*exp(u[2])/exp(u[3])-r_2
  du[4] = -(beta_v/N_p)*exp(u[3])-gamma+(1-theta)*mu/exp(u[4])
  du[5] = (beta_v/N_p)*exp(u[4])*exp(u[3])/exp(u[5])-gamma+theta*mu/exp(u[5])
end

function G_Diffusion(du,u,p,t)
  du[1] = (sigma_L/N_p)*exp(u[2])+(sigma_I/N_p)*exp(u[3])
  du[2] = -(sigma_L/N_p)*exp(u[1])
  du[3] = -(sigma_I/N_p)*exp(u[1])
  du[4] = -sigma_v
  du[5] = -sigma_v
end

prob_sde_tomato_sys = SDEProblem(F_Drift,G_Diffusion,u_0,time)
sol = solve(prob_sde_tomato_sys,SROCK1(),dt=0.00001)

title = plot(title = "R_s =$Rs0", grid = false, showaxis = false, bottom_margin = -50Plots.px)
p1=plot(exp.(sol),vars=(1),color="darkgreen",title="Susc. p.")
p2=plot(sol,vars=(2),color="darkorange", title ="Lat. p.")
p3=plot(sol,vars=(3),color="darkred",title = "Infec. p.")
p4=plot(sol,vars=(4),color="green", title = "Susc. v.")
p5=plot(sol,vars=(5),color="red",title ="Infec. v.")

plot(p1,p2,p3,p4,p5,title,layout = @layout([ [A B C]; [D E F]]), label="")
