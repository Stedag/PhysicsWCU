%% 
% Language: Matlab
% Author: Stephen Dages

% I was inspired by <https://www.youtube.com/watch?v=yf-PF1kkPNU https://www.youtube.com/watch?v=yf-PF1kkPNU> 
% to model our homework here


a=1;
nmax=30
spacial_res=50;
hbar=6.582e-16; %eVs
m=.511e6; %eV
Trev=(4*m*a^2)/(pi*hbar)
Trange=Trev;
temporal_res=150;
[X,T]=meshgrid(0:a/spacial_res:a,0:Trange/temporal_res:Trange);
syms psin(x,t)
syms c(n)
c(n) =-4*(sqrt(6)*(-1^((n-1)/2)))*sin(n*pi/2)/((n^2)*pi^2);
syms Psi_n(n,x,t)
Psi_n(n,x,t)= c(n)*(sqrt(2/a))*sin((n*pi/a)*x)*exp(-1i*pi^2*hbar*n^2*t/(2*m*a^2));
syms Psi(x,t)
Psi(x,t)=0
for nval=1:2:nmax
    Psi(x,t)= Psi(x,t)+ Psi_n(nval,x,t)
end

%%
syms P re
P(x,t)=conj(Psi(x,t))*Psi(x,t)
re(x,t)=real(Psi(x,t))
fun=@(xtest) double(P(xtest,0))
integral(fun,0,a)

fun=@(xtest) double(P(xtest,0))
integral(fun,0,a)
%%
surf(X,T,double(P(X,T)))
xlabel("x")
ylabel("t")
%%
surf(X,T,double(re(X,T)))
xlabel("x")
ylabel("t")
%%
x=linspace(0,a,spacial_res);
plot(x,Psi(x,0))
plot(x,P(x,0))