function [t,y] = tacs_YCCPlosBIOL(startt,stopt, k, relph, relf,osfrq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function from Krause et al Plos Biology 2022 https://doi.org/10.1371/journal.pbio.3001650
% Modified by Y. Cabral-Calderin to have the oscillator frequency "osfrq"
% as input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solve the Hopf-Andronov equation over the time interval from startt to stopt,
% for external input of magnitude k, and phase relative to the ongoing oscillation relph.
% relf is frequency relative to ongoing oscillation (multiple)
% external stimulation turns on at time STIMON
% 
% Returns time indices t and output y (a vector containing the outputs of
% the two variables). You can then plot(t,y(:,1)), for example.
% Parameter values from the paper by Doelling & Assaneo:
lambda = 0.2;
%omega = 2*pi*0.5; % 1 cycle per unit time (works well with default parameters). 
omega = 2*pi*osfrq; % 1 cycle per unit time (works well with default parameters). 

gamma = 1;

[t,y] = ode45(@(t,y) odefcn(t,y,lambda,omega,gamma,k,relph, relf),[startt stopt],[0; -1]);
end

function dydt = odefcn(t,y,lambda,omega,gamma,k,relph, relf)
%Hopf-Andronov equation
STIMON = 10; % time step when tACS turns on

if t<STIMON 
    s=0;
else
    s = sin(t*(omega*relf)+relph); % setting the initial conditions to [0,-1] produces sine phase for y(1)
%    s=k; %tdcs
end

h = (y(1).^2+y(2).^2);

dy1 = lambda*y(1) - omega*y(2)-gamma*h.*y(1);
dy2 = lambda*y(2) + omega*y(1)-gamma*h.*y(2);

dydt = [dy1+k*s; dy2];
end
