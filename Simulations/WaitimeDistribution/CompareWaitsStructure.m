clear
close all
% read in sim_1.log of the form "Sample SIREvents Lineages" as a string seperated by \t and with one header line
vals=[];

fwait = fopen('structWaitTimes.tsv', 'w');

for iteration = 1 : 500
    % read in the file fopen(['../SIR/master/SIR_simulations_' i '.xml'], find the lines  transmissionRate="......" and populationSize="..." and set the 
    % parameters transmission_rate and popS_Size to the values found
    % in the file
    
    f = fopen(['master/structuredSIR_simulations_' num2str(iteration) '.xml']);
    while ~feof(f)
        line = fgetl(f);
        if contains(line, '<transmissionRate spec="RealParameter" value="')            
            tmp = strsplit(line, '"');
            % convert the vector of transmission rates in tmp{4} of the form 1.2 1.4 1.2 into a vector transmission rate
            tmp2 = strsplit(tmp{4});
            transmission_rate = str2double(tmp2);
        elseif contains(line, '<populationSize spec="IntegerParameter"')
            tmp = strsplit(line, '"');
            tmp2 = strsplit(tmp{4});
            pop_Size = str2double(tmp2);
        end
    end
    fclose(f);

    
    data = textscan(fopen(['master/structuredSIR_simulations_' num2str(iteration) '.log']), '%s', 'delimiter', '\t', 'headerlines', 1);
    
    % the SIR events are the 2:3:end entries in the data cell array
    SIR_events = data{1}(2:3:end);
    % the lineages through time are the 3:3:end entries in the data cell array
    lineages = data{1}(3:3:end);
    
        
    % loop over the SIR and lineages
    for i = 1:length(SIR_events)
        times_between_reassortment = [];
        avg_rate_reassortment = [];
        avg_lins = [];
        avg_I = [];
        co_inf_events = [];
        avg_rate = [];

        disp(i)
        % replace ( and ) then split the SIR events into a cell array of strings and 
        sir = strsplit(strrep(strrep(SIR_events{i}, '[', ''), ']', ''), ', ');
        % do the same for lineages
        lins = strsplit(strrep(strrep(lineages{i}, '[', ''), ']', ''), ', ');
        % compute the number of infected individual over time
        I = zeros(1,20);
        time_I = 0;
        % get the initial location
        tmp = strsplit(sir{1}, ':');
        tmp2 = strsplit(tmp{end}, '_');
        I(str2double(tmp2{1})+1) = 1;
    
        % keep also track of the co-infection event times
        co_inf_times = [];
    
        % 0 is a co-infection event (no change to I), 1 a transmission event, 2 a recovery event and 2 a loss of immunity event
        for j = 1:length(sir)
            tmp = strsplit(sir{j}, ':');
            if strcmp(tmp{1}, '0')
                I = [I; I(end,:)];
                time_I = [time_I, str2double(tmp{2})];
                co_inf_times = [co_inf_times, time_I(end)];
            elseif strcmp(tmp{1}, '1')
                I = [I; I(end,:)];
                I(end, str2double(tmp{end})+1 )=I(end, str2double(tmp{end})+1 )+1;
                time_I = [time_I, str2double(tmp{2})];
            elseif strcmp(tmp{1}, '2')
                I = [I; I(end,:)];
                I(end, str2double(tmp{end})+1 )=I(end, str2double(tmp{end})+1 )-1;
                time_I = [time_I, str2double(tmp{2})];
            elseif strcmp(tmp{1}, '3')
                I = [I; I(end,:)];
                tmp2 = strsplit(tmp{end}, '_');
                I(end, str2double(tmp2{2})+1 )=I(end, str2double(tmp2{2})+1 )+1;
                I(end, str2double(tmp2{1})+1 )=I(end, str2double(tmp2{1})+1 )-1;
                time_I = [time_I, str2double(tmp{2})];
            end
        end
    
        % for lineages, 0 is a sampling event (+1 no_lins), 1 a reassortment event (+1 no_lins) and 2 a coalescent event (-1 no_lins)
        no_lins=zeros(1,20);
        tmp = strsplit(lins{1}, ':');    
        time = str2double(tmp{2});
        no_lins(str2double(tmp{end})+1) = 1;
    
        % keep track of how long it takes between reassortment events and how many lineages there were between them, weighted by the time between them
        time_since_reassortment = 0;
        % lin_time_since_reassortment = 0;
        rate_since_last_interval = 0;
        rt=[];

        disp(length(rate_since_last_interval))
        
        % for each interval, this keeps track of the probability of
        % reassortment from the co-infection probability and the 
        tot_reassortment_time = 0;
        for j = 2:length(lins)
            tmp = strsplit(lins{j}, ':');
            time = [time, str2double(tmp{2})];
            % lin_time_since_reassortment = lin_time_since_reassortment + no_lins(end)*(time(end) - time(end-1));
            time_since_reassortment = time_since_reassortment + (time(end) - time(end-1));   
    
            % for the times between this and the last reassortment event,
            % compute the average number of infected individuals
            curr_time = time(end);
            last_time = time(end-1);
    
            % find the index of the last time that is smaller or equal than the current time           
            [~,idx1]=min(abs(time_I-last_time));
            idx2 = find(time_I==curr_time);
    
            % count the number of co-infection events in this interval
            transmission = (I(idx2:idx1-1,:)-1).*transmission_rate./pop_Size;
            ratelins = transmission.*no_lins(end,:);
            timediffs = diff(time_I(idx2:idx1));

            result = sum(ratelins'*timediffs');
            if length(result)==0
                result = 0;
            end

            rate_since_last_interval = rate_since_last_interval + result;

            % add the next lins
            no_lins = [no_lins;no_lins(end,:)];
            if strcmp(tmp{1}, '0')                
                no_lins(end, str2double(tmp{end})+1) = no_lins(end, str2double(tmp{end})+1)+1;
            elseif strcmp(tmp{1}, '1')
                no_lins(end, str2double(tmp{end})+1) = no_lins(end, str2double(tmp{end})+1)+1;


                times_between_reassortment = [times_between_reassortment, time_since_reassortment];
                % avg_lins = [avg_lins, lin_time_since_reassortment/time_since_reassortment];
                avg_rate = [avg_rate, str2double(tmp{3})*rate_since_last_interval/time_since_reassortment];

                if length(avg_rate)~=length(times_between_reassortment)
                    fd
                end

                rt=[rt, time(end)];
                % for the times between this and the last reassortment event,
                % compute the average number of infected individuals
                curr_time = time(end);
                last_time = time(end)-times_between_reassortment(end);
    
                % find the index of the last time that is smaller or equal than the current time           
                [~,idx1]=min(abs(time_I-last_time));
                idx2 = find(time_I==curr_time);
    
                % count the number of co-infection events in this interval
                co_inf_events = [co_inf_events str2double(tmp{3})*sum(co_inf_times<last_time & co_inf_times>=curr_time)];
                avg_I = [avg_I sum(I(idx2:idx1).*time_I(idx2:idx1))/sum(time_I(idx2:idx1))];
                time_since_reassortment=0;
                lin_time_since_reassortment=0;
                rate_since_last_interval=0;
            elseif strcmp(tmp{1}, '2')
                no_lins(end, str2double(tmp{end})+1) = no_lins(end, str2double(tmp{end})+1)-1;

            else
                tmp2 = strsplit(tmp{end}, '_');
                no_lins(end, str2double(tmp2{2})+1) = no_lins(end, str2double(tmp2{2})+1)+1;
                no_lins(end, str2double(tmp2{1})+1) = no_lins(end, str2double(tmp2{1})+1)-1;


            end       
        end       
    
        % times_between_reassortment = [times_between_reassortment, time_since_reassortment];
        % avg_lins = [avg_lins, lin_time_since_reassortment/time_since_reassortment];
        % rt=[rt, time(end)];
        % % for the times between this and the last reassortment event,
        % % compute the average number of infected individuals
        % curr_time = time(end);
        % last_time = time(end)-times_between_reassortment(end);
        % 
        % % find the index of the last time that is smaller or equal than the current time           
        % [~,idx1]=min(abs(time_I-last_time));
        % idx2 = find(time_I==curr_time);
        % 
        % % count the number of co-infection events in this interval
        % co_inf_events = [co_inf_events sum(co_inf_times<last_time & co_inf_times>=curr_time)];
        % avg_I = [avg_I sum(I(idx2:idx1).*time_I(idx2:idx1))/sum(time_I(idx2:idx1))];
    end
    % remove the first interval, which is not properly accounted for in time
    % above
    % times_between_reassortment(first_element) = [];
    % avg_lins(first_element) = [];
    % avg_I(first_element) = [];
    % co_inf_events(first_element) = [];
    
    clear I time_I no_lins time lin_time_since_reassortment time_since_reassortment sir SIR_events
    
    %% Do the summary and relate the rates
    
    % the probability of a reassortment event being observed in an interval is
    % equal to the probability of there being a co-infection event divided by
    % the number of infected indviduals during that time % the probability of observing a coI
    % p_obs_coI = co_inf_events.*avg_lins./avg_I;
    
    
    % scatter(-times_between_reassortment, avg_rate)
    % set(gca, 'XScale', 'log');
    % set(gca, 'YScale', 'log');
    
    
    % times_between_reassortment
    for j =1:length(times_between_reassortment)    
        fprintf(fwait, '%d\t%f\n', iteration, times_between_reassortment(j)*avg_rate(j));
    end

    vals = [vals, times_between_reassortment.*avg_rate];
    fprintf('mean = %.3f std = %.3f\n', mean(vals),std(vals));
    fprintf('total reassortment events in network = %.3f\n', length(times_between_reassortment))
end


ksdensity(vals); hold on
ksdensity(exprnd(1,1000000,1)); 
legend('wait time distribution', 'exponential distribution')

fclose(fwait);


% compare the total number of reassortment event and the probability of
% observing those
% fprintf('total events = %d\nexpected observations = %.2f\n',...
%     length(times_between_reassortment),sum(p_obs_coI))
% c = 1;
% for i = floor(length(vals)/100):floor(length(vals)/100):length(vals)
%     [h(c), p(c), ksstat(c)] = kstest2(vals(1:i),vals2);
%     c=c+1;
% end
% figure()
% plot(ksstat)


% das
% reassortment_times = abs(reassortment_times);
% % 
% close all
% 
% % Calculate x and y values
% x = co_inf_rates;
% 
% y = 1 ./ (reassortment_times .* reassortment_lins);
% 
% scatter(co_inf_rates, reassortment_lins./reassortment_times, '.')
% 
% % Create a scatter plot with log scales
% % scatter(x, y);
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');
% 
% % Perform linear regression
% coefficients = polyfit(log(x), log(y), 1); % Fit a linear model to log-transformed data
% 
% % Extract slope and intercept from the coefficients
% slope = coefficients(1);
% intercept = coefficients(2);
% 
% % Create the regression line using the fitted coefficients
% x_fit = logspace(log10(min(x)), log10(max(x)), 100); % Generate x-values for the line
% y_fit = exp(intercept) * x_fit.^slope; % Calculate corresponding y-values
% 
% % Plot the regression line
% hold on; % Keep the current plot
% plot(x_fit, y_fit, 'r', 'LineWidth', 2); % Plot the regression line in red

% 
% 
% % hold on; plot(reassortment_times, 1./reassortment_lins, 'b.');
% % ksdensity((reassortment_lins.*co_inf_rates)./reassortment_times)
% % plot the same as a histogram plot
% % hist(-1*reassortment_times.*reassortment_lins.*co_inf_rates, 100)
% % add a plot with exponential random numbers with mean 1
% 
% close all
% ksdensity(reassortment_times.*reassortment_lins.*co_inf_rates);
% hold on; 
% exp_vals = exprnd(0.813,10000,1);
% ksdensity(exp_vals)
% 
% 
% 
% scatter(co_inf_rates.*reassortment_lins, 1./reassortment_times);
% 
