%%  Tutorial
%%%%%%%%%%%%%%%%%%%%%
%% Comments : 
%The values that can be modified by the user are indicated by 
% the comment  %% THIS VALUE CAN BE CHANGED BY THE USER 
% All other lines should be left as they are by generic users
% NOTHING SHOULD BE CHANGED BELOW LINE 88
% Experienced users are of course welcome to modify the code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
epsilon=0;
noise_init=0;

%%% size of the sample (initial number of particles)
% Set N=0 by default: the sample is infinite 
%N=200; 
N=0;      %% THIS VALUE CAN BE CHANGED BY THE USER


%%%%
%% 3 ways to use the code: scenario 1, 2 or 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% In all the scenarii
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Choose the time points  at which the distribution profile is calculated
%time=[0 0.1 0.5 1 2 5 8 13 18];
%time=[0 1 2 5 8 13 18];
time=[0 1 2];  %% THIS VALUE CAN BE CHANGED BY THE USER
%%% The time at which the data are gathered can be in any unit (for instance in seconds).
%The unit of the time points t affects the unit of the parameters alpha and gamma.




%% 1. To generate the distribution profiles corresponding to fixed parameters alpha, gamma, and kappa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Set scenario= 1 

%%% Choose here the parameters gamma and alpha
gamma=1.3;   %% THIS VALUE CAN BE CHANGED BY THE USER
alpha=1;     %% THIS VALUE CAN BE CHANGED BY THE USER


%%% Choose here your initial condition 
choice_init=2;   %% THIS VALUE CAN BE CHANGED BY THE USER
                 %%% To see how,  look at the auxilliary function intial
                 %%% condition_paper
                 
                 
%%% Choose here your division kernel kappa                    
choicek0=5;      %% THIS VALUE CAN BE CHANGED BY THE USER
R=0.01;          %% THIS VALUE CAN BE CHANGED BY THE USER
                 %%% To see how,  look at the auxilliary function k0_paper



                 
%% 2. To extract the values of alpha and gamma from your experimental distribution profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set scenario =2

%%% Enter here a vector where is stocked the set of siezs that was measured.
SET_OF_SIZES = [1 2 4 8 16];   %% THIS VALUE CAN BE CHANGED BY THE USER

%%% Enter here the proportion of particules (the frequencies)for the
%%% different time points. 
%%% Here, the length of vector time is 28
FREQUENCIES{1}= [0 0 0 0 12];
FREQUENCIES{2}= [0 0 0 6 6];
FREQUENCIES{3}= [0 0 3 6 3];
.......
FREQUENCIES{28}= [];
    


%% 3. To generate the distribution profiles corresponding to fixed parameters alpha, gamma, and kappa
%% and then extract the values of alpha and gamma from your simulated distribution profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set scenario = 3
%%% AND choose the values for alpha, gamma, kappa, and the initial condition
%%% directly in scenario 1. 
%% 



scenario=2;  %% THIS VALUE CAN BE CHANGED BY THE USER

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NOTHING SHOULD BE CHANGED BELOW THIS LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T=length(time);




if scenario==1 | scenario==3


%%  Simulate data
%%%%%%%%%%%%%%%%%%%%
Z=0:0.01:1; 
dt=0.1;
fminOpt = optimset('Display', 'off', ...
    'MaxFunEvals', 1000000, 'MaxIter', 5000, ...
    'TolX', 1e-20, 'TolFun', 1e-20);


[X,u,Y,t_store,mass]=DETX_simu_direct_paper(alpha,gamma,choicek0,R,choice_init,N,time,dt,0,0);

figure(1)
hold on;
for i=1:length(time)
h(i)=plot(X,u{i},'LineWidth',2);
legh{i}=sprintf('t= %d',time(i));
end
legend(h,legh);
box on
xlabel('size x','fontsize',16)
ylabel('distribution u(t,x)','fontsize',16)
xlim([0,max(X)])
 set(gca,'FontSize',16);

figure(2)
hold on;
for i=1:length(time)
hh(i)=plot(Y{i},time(i)^(-2/gamma)*u{i},'LineWidth',2);
leghh{i}=sprintf('t= %d',time(i));
end
box on
legend(hh,leghh);
xlabel('rescaled size xt^{1/\gamma}','fontsize',16)
ylabel('rescaled distribution t^{2/\gamma}u(t,xt^{1/\gamma})','fontsize',16)
xlim([0,max(Y{2}/2)]);
 set(gca,'FontSize',16);
%xlim([0,max(Y{2}/10)]);


   f=u;
for i=1:T
    Numb(i)=moments(X,u{i},0);  %% Numb= total number of particles
    f{i}= u{i}/Numb(i);
            NN(i)=floor(N*Numb(i)/Numb(1));
       for j=1:NN(i)
      Y{i}(j)=acc_rej_algo(X,f{i});
    end;
    M(i)= 1/NN(i)*sum(Y{i}.^1);
end

end;


if scenario==2
    X=SET_OF_SIZES;
   for j=1:T
   u{j}=FREQUENCIES{j};
end;
    
%        f=u;
% for i=1:T
%     Numb(i)=moments(X,u{i},0);  %% Numb= total number of particles
%     f{i}= u{i}/Numb(i);
%            NN(i)=floor(N*Numb(i)/Numb(1));
%        for j=1:NN(i)
%       Y{i}(j)=acc_rej_algo(X,f{i});
%     end;
%     M(i)= 1/NN(i)*sum(Y{i}.^1);
% end


end;
if scenario==2 | scenario==3

%% 2. Recover the parameters gamma and alpha from your data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% For this tutorial, we use simulated data (X,u) we just obtained. 
%%% They can be replaced by your data
%%% X= the tabular of the sizes of your particles
%%%% write in the code intial_condition_paper the intial frequency 
%%%% of hte sizes (the first that was measured)
%%% u{1}= the frequency of each size at time t(1)
%%% u{2}= the frequency of each size at time t(2)
%%%.....

 [gamma_e,alpha_e,Te,C]=DETX_simu_inverse_paper(X,u,time,N,1,0);
 CS=10^C;

 
 
 %%%% build the theoretical fit found by the code
 time2= linspace(0,Te,10);
 size2=length(time2);
 time3=linspace(Te,time(end),10);
 curve2=  CS*Te^(-1/gamma_e) ;
 curve3=  CS*time3.^(-1/gamma_e) ;

%    figure(3)
%   hold on;
%   box on
%   plot(time,M,'r+','LineWidth',2);
%   plot(time,M,'b-','LineWidth',2);
%   plot(time2,curve2*ones(1,size2) ,'g-','LineWidth',2);
%  plot(time3,curve3 ,'g-','LineWidth',2);
%   %plot(log10(Te),0,'+');
%   xlabel('time t','fontsize',16);
%   ylabel('M_1(t)','fontsize',16);
%  set(gca,'FontSize',16);

 
 name = 'Alice';   

TODISP = ['The estimated value for gamma is ',num2str(gamma_e),'.'];
disp(TODISP);
TODISP_A = ['The estimated value for alpha is ',num2str(alpha_e),'.'];
disp(TODISP_A);

end;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Auxilliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial condition 
%
%
%
% Usage:  [Z,U] = initial_condition_paper(choice_init)
%
% Arguments:
%       choice_init    - integer  that specifies the type of
%       initial condition

% Returns:
%       Z             -logarithmic grid on which the initial condition is
%                       evaluated
%       U             - f_0(X)=U vector of the initial condition with total
%       mass 1



function [Z,U] = initial_condition_paper(choice_init)

if choice_init == 1 %decreasing exponential
    Z=-10:0.01:2; 
    X= exp(Z);
    U=exp(-X);
    
elseif choice_init == 2 %one picked Gaussian
    S=0.02;
    Z=-10:0.01:1; X= exp(Z);
    U=2/(sqrt(pi)*sqrt(S))*exp(-(X-1/2).^2/S);
    U=2/(sqrt(pi)*sqrt(S))*exp(-(X-1).^2/S);
    
elseif choice_init == 3 % two picked gaussian
    S=0.01;
    Z=-10:0.01:1; X= exp(Z);
    U= 2/(2*sqrt(pi)*sqrt(S))*exp(-(X-1/3).^2/S)+ 2/(2*sqrt(pi)*sqrt(S))*exp(-(X-2/3).^2/S);

 elseif choice_init== 4
    load('aSyn_DETX.mat')
    [U,Z,h]=ksdensity(log(L{15}),'weights',w{15});
    X=exp(Z);
%%%% INCLUDE HERE INITIAL CONDITIONS CORRESPONDING TO EXPERIMENTAL DATA
    %[U,Z,h]=ksdensity(log(L{15}),'weights',w{15});
    % X=exp(Z);
    %% THIS VALUE CAN BE CHANGED BY THE USER
    %% Other intial conditions can be considered
end

% 
%% Normalization of the initial condition: integral= 1
I=length(Z);
Y=X; Y(1)=[];
V=U;
V(I)=[];
X(I)=[];
dX=Y-X;
U=U/sum(V.*X.*dX);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k0_paper  
%
% Fragmentation kernel
%
% Usage: [y] = k0_paper(x,choice_k0,R)
%
% Arguments:
%       x            - grid of [0,1]  on which we evaluate the kernel k0
%       choice_k0    - integer between 1 and 8 that specifies the type of
%                       kernel
%       R            - spreading parameter when the kernel k0 is a sum of gaussians
%
% Returns:
%       y            - y=k0(x)

%Remark: R needs only to be specified for the gaussian kernel and the boundary kernel, 
% where it specifies the sreading. R is not equal to zero.
%In all other cases, the value chosen for R has no impact.

function [y] = k0_paper(x,choice_k0,R)  

if choice_k0==2    %kernel charging the boundaries
    y=2*(1+cos(2*pi*x)).^(R)./integral(@(z) ((1+cos(2*pi*z)).^R),0,1);
elseif choice_k0==3   %uniform kernel
      y=2*ones(1,length(x)); 
elseif choice_k0==4    %one peak Gaussian kernel
          y=kernel_gaussian(x,R);         
  elseif choice_k0==5   %two peaks Gaussian kernel
       y=kernel_2gaussian(x,R); % 2 gaussians in 1/3 and 2/3 :
     %% THIS VALUE CAN BE CHANGED BY THE USER
    %% Other fragmentation kernels can be considered
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function y=kernel_gaussian(x,R)   
%C=2/(sqrt(pi)*sqrt(R));
y= 2*exp(-(x-1/2).^2/R)/ integral(@(z) (exp(-(z-1/2).^2/R)),0,1);
end
    
function y=kernel_2gaussian(x,R)  
%C=1/(sqrt(pi)*sqrt(R));
y= 2*(exp(-(x-1/3).^2/R)+exp(-(x-2/3).^2/R))/integral(@(z) (exp(-(z-1/3).^2/R)+ exp(-(z-2/3).^2/R)),0,1);
end

%DETX_simu_direct

%solver for the fragmentation equation .
%
%
% Usage:  [X,u,Y,t_store,mass]=DETX_simu_direct_paper(alpha,gamma,choicek0, R,choice_init, N,time,dt,doPlot, doVideo);


%
% Arguments:
 %        alpha       - positive real number
 %        gamma        - positive real number
%         choicek0      - integer between 1 and 5 that specifies the type of
%                      fragmentation kernel
%                      1= delta, 2= edge ,  3= uniform, 4= 1Gaussian, 5= 2
%                      Gaussians
 %                     If k0 is a Gaussian or a cosinus kernel (2 or 4), its spreading R should be specified.
 %                     if k0 is not a Gaussian, the parameter R is useless. 


%          choice_init     -   choice of initial condition:
                           %1 = decreasing exponential
                           %2 = 1 gaussian
                           %3 = 2 gaussians
                           %4 = uniform
                           %5 = adapted to the paper of Ziff/grady
                           %6 = decreasing exponential, grid refined around
                           %zero
                          
% N - the sampling noise    % If N>0, we pick a sample of size N distributed along the initial condition and 
                            %build a new initial condition using ksdensity. If N=0, the initial condition is left unchanged.
                           
%              time         vector of any size of ordered non negative real numbers
%                           It represents the times at which the solution is computed
%                           It is assumed that the initial condition is given at time t=0.
%                           If t=0 does not belong to the set of times, the initial condition will
%                           not be plotted
%
%               dt          % upper bound for the time step. It may be
%                           automatically
%                           diminished if dt does not meet the CFL condition
%
% doPlot and doVideo        %booleans (1 or 0) that specifies if the plots or a video should be diplayed


% Returns:
%                X        - grid where the solution f is defined
%                  u       - array where are stocked the distributions profiles f for different the
%                             times of the vector time
%                  Y       -rescaled grid, adapted to thee convergence to
%                             stationary state
%            t_store        - actual times where the solution is computed (in case where
%                        for some i, t(i) / Delta t is not an integer
%              mass        -vector of the mass of the system for the computed times (should be constant)
%                          useful if we want to plot the rescaled stationary profile



function [X,u,Y,t_store,mass]=DETX_simu_direct_paper(alpha,gamma,choicek0,R,choice_init,N,time,dt,doPlot, doVideo);

if doVideo==1
vilm=VideoWriter('movie.avi');
open(vilm);
GF=[];
figure(1);clf;
end;

%initial distribution on a log scale and its uniform grid
[Zin,Uin]= initial_condition_paper(choice_init); % in the logarithmic scale
Uin_old=Uin;

% Adding the sampling noise
if N>0
Xin=exp(Zin);
for i=1:N
 Y(i)= acc_rej_algo(Xin,Uin);
end;
[Uin,Xin]= ksdensity(Y,Xin);
Uin=Uin*moments(Xin,Uin_old,0); %rescaling since the intial distribution may be of integral not equal to 1
end;

%renormalised distribution at the given times
%on the grid Y=exp(Xin)

%1. Choice of delta_t: CFL condition + condition on log(2)
Z=Zin; %uniform grid
I=length(Z);
dz =Z(2)-Z(1); %space step
%dt=min(time(1)+ (time(1)==0)*time(2),1/alpha*exp(-gamma*max(Z))); %%dt satisfies the CFL condition, and we take it even smaller if we want to evaluate the solution for very small times.
if time(1)==0&length(time)>1
dt=min(time(2),dt);
elseif time(1)==0
    dt=dt;
else
 dt=min(time(1),dt);
end;
t=0;
tmin=time(1);
%time=time-tmin;
n_times=1; 
t_store(1)=tmin;
%n=n/sumn; %normalisation so that integral (n dy) =1 but uses the regular grid to be more precise
q= zeros(I,I);
%tmax=time(2);
u=cell(size(time));
v=cell(size(time));
Y=cell(size(time)); %adapted grid for the asymptotic rescaling

    kk=1;

 X=exp(Z);
n=Uin.*X.*X;  
mass(1)=sum(n)*dz;% n is the conservative variable
sumn=sum(n.*dz); 
Y{1} = exp(Z+(1/gamma)*log(tmin));
u{1}=n./X./X; 
%number of equally spaced points
      while n_times<=length(time)
           % Storage of the computed profiles
     if abs(t-time(n_times))<min(abs(t-dt-time(n_times)),abs(t+dt-time(n_times)))
         u{n_times}=n./X./X;
        mass(n_times)=sum(n)*dz;
        Y{n_times} = exp(Z+(1/gamma)*log(time(n_times)+tmin));
         t_store(n_times)=t;
        n_times=n_times+1;
        end
     for i=1:I
     q(i,1:I-i)= exp(-2*(1:I-i)*dz).*k0_paper(exp(-(1:I-i)*dz),choicek0,R).*n(i+1:I).*exp(gamma*(1:I-i)*dz)*dz;
     F(i)= sum(q(i,:));  
     end
    n(1:I)=(n(1:I)+dt*alpha.*exp(gamma*Z(1:I)).*F(1:I))./(1+alpha*dt*exp(gamma*Z(1:I)));
    n=sumn*n./(sum(n)*dz); %normalisation, though the scheme should be almost conservative
   % if (doVideo==1)&&(n_times/50==n_times/50)  
     if (doVideo==1)
     % plot(X,n./X./X)
      YY = exp(Z+(1/gamma)*log(t+tmin))
      plot(YY,t^(-2/gamma)*n./X./X )
      xlim([0,max(YY)+0.1]);
      GF=[GF getframe];
    end;
t=t+dt;
kk=kk+1;
    
    
end

%% Plot the figures

if doPlot==1 
    
figure(1)
subplot(2,3,1)
hold on;
for i=1:length(time)
h(i)=plot(X,u{i});
legh{i}=sprintf('t= %d',time(i));
end
legend(h,legh);
box on
xlabel(' x (size)')
ylabel('f(t,x)')
xlim([0,max(X)])

if length(time)>1|time(1)>0
subplot(2,3,2)
hold on;
for i=1:length(time)
hh(i)=plot(Y{i},time(i)^(-2/gamma)*u{i});
leghh{i}=sprintf('t= %d',time(i));
end
box on
legend(hh,leghh);
xlabel('xt^{1/\gamma}')
ylabel('t^{2/\gamma}f(t,xt^{1/\gamma})')
xlim([0,max(Y{2}/2)]);
%xlim([0,max(Y{2}/10)]);
end

%%% Mass evolution in log-log scale

for i=1:length(time)
    Numb(i)=moments(X,u{i},0);  %% Numb= total number of particles
    f{i}= u{i}/Numb(i);
     M(i)=moments(X,f{i},1);
end
 
if N>0
 subplot(2,3,3)
  box on
  hold on;
  plot(X,Uin_old,'r-');
  plot(X,Uin,'b-');
  xlabel('x');
  ylabel('f(0,x)');
  title('t=0');
  legend('initial distribution','noisy initial distribution')
elseif N==0
     subplot(2,3,3)
  box on
  hold on;
  plot(X,Uin,'b-');
  xlabel('x');
  ylabel('f(0,x)');
  title('t=0');
  legend('initial distribution')
end;

  subplot(2,3,4)
  box on
  plot(0:0.01:1,k0_paper(0:0.01:1,choicek0,R));
  ylabel('k_0');
  
  
  subplot(2,3,5)
  hold on;
  box on
  plot(log(time),log(M),'r+');
  plot(log(time),log(M),'b-');
  xlabel('log(t)');
  ylabel('log(\mu((t))');
  
  subplot(2,3,6)
  hold on;
  box on
  plot(log(1+time/dt),log(M),'r+');
  plot(log(1+time/dt),log(M),'b-');
  xlabel('log(1+t/dt)');
  ylabel('log(\mu((t))');
  
 
end;    

if doVideo==1
 writeVideo(vilm,GF);
 close(vilm);
end;

end

%DETX_simu_inverse

%solve the inverse problem of the fragmentation equation
%
%
% Usage:  [alpha,gamma,k0]=DETX_simu_inverse_paper(X,u,time,N,s,epsilon);


%

%% Arguments:
%                X        - grid where the solution f is defined
%                u        - array where are stocked the distributions profiles f for different the
%                             times of the vector time
%                
%              time       - vector of any size of ordered non negative real numbers
%                           It represents the times at which the solution is computed
%                           It is assumed that the initial condition is given at time t=0.
%                 s      - positive real number. Moment we use to compute gamma.    


%% Returns 
 %        alpha       - positive real number
 %        gamma        - positive real number
%         k0            -fragmentation kernel
%         Te          -   the time at which the system is assumed to be at
%         equilibrium          1
%               N         - (large) positive integer : size of the sample. By convention, if N=0, we assume that we have a perfect knowledge 
%                            of the distribution function. 
%                            If N<0, then we assume that the size of the
%                            sample depends on the total number of
%                            particles. 
%              epsilon    - (small) positive integer : quantifies the
%                            deterministic noise

function [gamma,alpha,Te,C]=DETX_simu_inverse_paper(X,u,time,N,s,epsilon);

T=length(time);
L= length(X);
%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=u;
for i=1:T
    Numb(i)=moments(X,u{i},0);  %% Numb= total number of particles
    f{i}= u{i}/Numb(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fminOpt = optimset('Display', 'off', ...
    'MaxFunEvals', 1000000, 'MaxIter', 5000, ...
    'TolX', 1e-20, 'TolFun', 1e-20);


% Setting up optimisation
% initial values for parameter estimation


% sdt: steady distribution time
% empirically decide when to count as steady distribution
% sdt0 = time(end);
  gamma0 = 1;
  sdt0=time(1);
  Te=-1;
  Te=time(end)+1;
  logCs0=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% When we  only access a N size sample of data
if (N>0)
    %Y=f;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Generate N measurements from the distribution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    for i=1:T
       for j=1:N
      Y{i}(j)=acc_rej_algo(X,f{i});
    end;
   end;
   
    %%%%%%%%%
    %% gamma
    %%%%%%%%%
    
  %%%% Calculating the moments directly from the data 
    for i=1:T
    M(i)= 1/N*sum(Y{i}.^s);
    end
    


while(Te<time(1)|Te>time(end))
    sdt0= time(end)*rand(1);
 p = fminsearch(@FminFct, [gamma0 logCs0 sdt0], fminOpt);
    % Calc final parameters
    [~, msCalc, ~] = FminFct(p);
     gamma=p(1);
     Te=p(3);
     C=p(2);
end;

    
    %%%%%%%%%
    %% alpha
    %%%%%%%%%

    alpha=1/(time(T)*gamma*1/N*sum(Y{T}.^gamma));


  
elseif (N==0) 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% gamma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:T
    M(i)=moments(X,f{i},s);
end




logCs0 = log10(max(M));
logCs0 = log10(M(1));


while(Te<time(1)|Te>time(end))
    sdt0= time(end)*rand(1);
 p = fminsearch(@FminFct, [gamma0 logCs0 sdt0], fminOpt);
    % Calc final parameters
    [~, msCalc, ~] = FminFct(p);
     gamma=p(1);
     Te=p(3);
     C=p(2);
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% alpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


alpha=1/(time(T)*gamma*moments(X,f{T},gamma));





elseif N<0
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Generate N measurements from the distribution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   N=-N; 
   
    for i=1:T
        NN(i)=floor(N*Numb(i)/Numb(1));
       for j=1:NN(i)
      Y{i}(j)=acc_rej_algo(X,f{i});
    end;
   end;
   
    %%%%%%%%%
    %% gamma
    %%%%%%%%%
    
  %%%% Calculating the moments directly from the data 
    for i=1:T
    M(i)= 1/NN(i)*sum(Y{i}.^s);
    end
    
   logCs0 = log10(max(M));

%%%%%
while(Te<time(1)|Te>time(end))
    sdt0= time(end)*rand(1);
 p = fminsearch(@FminFct, [gamma0 logCs0 sdt0], fminOpt);
    % Calc final parameters
    [~, msCalc, ~] = FminFct(p);
     gamma=p(1);
     Te=p(3);
     C=p(2);
end;

    
    %%%%%%%%%
    %% alpha
    %%%%%%%%%

 alpha=1/(time(T)*gamma*1/NN(T)*sum(Y{T}.^gamma));

    
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Target function for optimisation
    function [ssr, msCalc, residuals] = FminFct(p)
        
        % Parameters
        gamma = p(1);
        %cs = 2.72.^p(2);
         cs = 10.^(p(2));
        sdt = abs(p(3));
        
        % Target function
        msCalc = cs.*(time.^(-s/gamma));
        % for t <= sdt, target function is a constant, objective method
        % but tradeoff is one extra (empirical) parameter
        msCalc(time <= sdt) = cs.*(sdt.^(-s/gamma));
        
        % Sums of squared residuals
    
        ssr = sum(((msCalc-M)).^2);
        % residuals
        residuals = msCalc-M;
    end

end
