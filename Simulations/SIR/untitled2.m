% Parameters
N = 10000;
beta = 0.2;      % Transmission rate
gamma = 0.1;     % Recovery rate
k = 0.1;         % Dispersion parameter for the negative binomial distribution

% Variables for storing the time series data
T = [];
S_pop = [];
I_pop = [];
R_pop = [];


% Counter for the total number of infections
total_infections = 0;

% Repeat the simulation until at least 100 infections in total
while total_infections < 100
    % Reset initial conditions
    S = N - 1;
    I = 1;
    R = 0;
    t = 0;
    total_infections = I; % Start counting infections from the initial infected individual

    % Store time and populations
    T_current = [t];
    S_current = [S];
    I_current = [I];
    R_current = [R];
    offsprings = [0];

    while I > 0 && S > 0
        % compute the current R
        R_val = (beta * S * I / N)/(gamma * I);
        % Calculate rates
        infection_rate = beta * S * I / N/ R_val;
        recovery_rate = gamma * I;
        total_rate = infection_rate + recovery_rate;

        % Time to next event
        tau = exprnd(1/total_rate);
        t = t + tau;

        % Determine the next event
        if rand < (infection_rate / total_rate)
            % Infection event: draw the number of new infections from the negative binomial distribution
            % sample the new number of infections from the negative binomial distribution, such
            % that the mean is R and the variance is k.
            R_nb = beta/gamma;

            num_new_infections = nbinrnd(R_nb*k, k/(k+R_nb));
            num_new_infections = min(num_new_infections, S); % Can't infect more than the number of susceptibles
            % for each num_new infections, pick if it hit a S or not, if it hit an S, then it's an infection
            tot_new_infections=0;
            if num_new_infections>0
                for i = 1:num_new_infections
                    if rand < S/N
                        tot_new_infections = tot_new_infections + 1;
                    end
                end
            end
            S = S - tot_new_infections;
            I = I + tot_new_infections;
            total_infections = total_infections + tot_new_infections; % Update total infections
            offsprings = [offsprings, tot_new_infections];
        else
            % Recovery event
            I = I - 1;
            R = R + 1;
        end

        % Update the arrays
        T_current = [T_current, t];
        S_current = [S_current, S];
        I_current = [I_current, I];
        R_current = [R_current, R];
    end

    % Check if the simulation had at least 100 infections, if not repeat
    if total_infections >= 100
        T = T_current;
        S_pop = S_current;
        I_pop = I_current;
        R_pop = R_current;
    end
end

% Plot the results
plot(T, S_pop, 'b', T, I_pop, 'r', T, R_pop, 'g');
legend('Susceptible', 'Infected', 'Recovered');
title('Stochastic SIR Model Simulation with Gillespie Algorithm and Negative Binomial Offspring');
xlabel('Time');
ylabel('Number of Individuals');
