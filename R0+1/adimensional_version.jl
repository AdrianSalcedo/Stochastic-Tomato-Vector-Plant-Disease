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
mu = 0.1#0.058596601
gamma = 0.6
sigma_L = 0.378173863
sigma_I = 0.030015876
sigma_v = 0.254146903
N_v = mu/gamma
r= max(r_1,r_2)
u_0 = [97.0,1.0,2.0,0.03,0.05]
T = 1000.0
time = (0.0,T)
dt=0.01
u0= [97.0/N_p,1.0/N_p,2.0/N_p,0.04,0.06]
N_p = u0[1]+u0[2]+u0[3]


Rs0 = beta_p*beta_v/(r*gamma)
print("Rs0=",Rs0)

function F_Drift(du,u,p,t)
 @inbounds begin
  du[1] = -beta_p*u[1]*u[5]+r_1*u[2]+r_2*u[3]
  du[2] = beta_p*u[1]*u[5]-b*u[2]-r_1*u[2]
  du[3] = b*u[2]-r_2*u[3]
  du[4] = -beta_v*u[4]*u[3]-gamma*u[4]+(1-theta)*mu
  du[5] = beta_v*u[4]*u[3]-gamma*u[5]+theta*mu
  end
  nothing
end

function G_Diffusion(du,u,p,t)
 @inbounds begin
  du[1,1] = sigma_L*u[2]*u[1]/N_p+sigma_I*u[1]*u[3]/N_p
  du[1,2] = 0
  du[2,1] = -sigma_L*u[1]*u[2]/N_p
  du[2,2] = 0
  du[3,1] = -sigma_I*u[1]*u[3]/N_p
  du[3,2] = 0
  du[4,1] = 0
  du[4,2] = -sigma_v*u[4]
  du[5,1] = 0
  du[5,2] = -sigma_v*u[5]
  end
 nothing
end



prob_sde_tomato_sys = SDEProblem(F_Drift,G_Diffusion,u0,time,noise_rate_prototype=zeros(5,2))

sol = solve(prob_sde_tomato_sys,Tsit5(),dt=dt)

title = plot(title = "R_s =$Rs0", grid = false, showaxis = false, bottom_margin = -50Plots.px)
p1=plot(sol,vars=(1),color="darkgreen",title="Susc. p.")
p2=plot(sol,vars=(2),color="darkorange", title ="Lat. p.")
p3=plot(sol,vars=(3),color="darkred",title = "Infec. p.")
p4=plot(sol,vars=(4),color="green", title = "Susc. v.")
p5=plot(sol,vars=(5),color="red",title ="Infec. v.")

plot(p1,p2,p3,p4,p5,title,layout = @layout([ [A B C]; [D E F]]), label="")


ensembleprob = EnsembleProblem(prob_sde_tomato_sys)
Sol = solve(ensembleprob,EnsembleThreads(),trajectories=500)

summ_1 = EnsembleSummary(Sol,0:0.01:T) #Media
plotly()
plot(summ_1,idxs=(1),labels="Middle 95%")
summ_2 = EnsembleSummary(Sol,0:0.01:T;quantiles=[0.25,0.75])
plot!(summ_2,idxs=(1),labels="Middle 50%",legend=true)





















