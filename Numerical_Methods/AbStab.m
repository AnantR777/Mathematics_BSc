
lambda = -20;
Tmax=1;
h=0.1;
N=Tmax/h;


% ODE: y' = lambda y, y(0)=1

%% Implicit Euler 

BE=zeros(N+1,1);
BE(1)=1;

for i=1:N

    BE(i+1)=BE(i)/(1-h*lambda);

end

%% Explicit Euler

FE=zeros(N+1,1);
FE(1)=1;

for i=1:N

    FE(i+1)=FE(i)*(1+h*lambda);

end


%% Plot the results

x=linspace(0,Tmax,N+1);
u_exact=exp(lambda*x);
plot(x,u_exact,'r-',...
    x,BE,'*-b',...
    x,FE,'o:m');