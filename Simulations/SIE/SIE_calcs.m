
clear
% define superinfection exclusion dynamics from https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001941
percentage_sucess = [0.5 0.5 0.5 0.3 0.25 0.25 0.1 0.025 0 0];
time = [0:length(percentage_sucess)-2 10000];
% define the recovery rate as an exponential decay with mean 4 days.
recovery_rate = 1/5;
% assume a constant rate of infection, and compute what proportion of infections would lead to a superinfection before the first infection is cleared.
infection_rate = 1;
% define the number of simulations
n_simulations = 1000000;
totProb = 0;
for r = 1:n_simulations
    % define the simulation time as an exponential random variable with rate recovery_rate 
    simulationTime = exprnd(1/recovery_rate)*24;
    % define what percentage of infections would have lead to a superinfection before the first infection is cleared.
    % by integrating percentage sucess from 0 to simulationTime
    totPercentage = 0;
    i=1;
    while time(i+1) < simulationTime
        totPercentage = totPercentage + percentage_sucess(i)*(time(i+1)-time(i));
        i=i+1;
    end
    totPercentage = totPercentage + percentage_sucess(i)*(simulationTime-time(i));
    totProb = totProb+totPercentage/simulationTime;
end

SIE_prob = 1/(totProb/n_simulations);





nsamples = logspace(0,5,1000);
prevalence = logspace(-5,-1,1000);

meanEvents=nsamples'*prevalence*150/14*7/365;


