clear
close all
% for all the log files in ../SIR/master, get the log files and compute the prevalence
log = dir('../SIR/master/*.log');
for i = 1:length(log)
    data = textscan(fopen(strcat('../SIR/master/', log(i).name)), '%s', 'delimiter', '\t', 'headerlines', 1);
    SIR_events = data{1}(2:3:end);
    for i = 1:length(SIR_events)
        disp(i)
        % replace ( and ) then split the SIR events into a cell array of strings and 
        sir = strsplit(strrep(strrep(SIR_events{i}, '[', ''), ']', ''), ', ');
        % do the same for lineages
        lins = strsplit(strrep(strrep(lineages{i}, '[', ''), ']', ''), ', ');
        % compute the number of infected individual over time
        I = 1;
        time_I = 0;
    
        % keep also track of the co-infection event times
        co_inf_times = [];
    
        % 0 is a co-infection event (no change to I), 1 a transmission event, 2 a recovery event and 2 a loss of immunity event
        for j = 1:length(sir)
            tmp = strsplit(sir{j}, ':');
            if strcmp(tmp{1}, '0')
                I = [I, I(end)];
                time_I = [time_I, str2double(tmp{2})];
                co_inf_times = [co_inf_times, time_I(end)];
            elseif strcmp(tmp{1}, '1')
                I = [I, I(end) + 1];
                time_I = [time_I, str2double(tmp{2})];
            elseif strcmp(tmp{1}, '2')
                I = [I, I(end) - 1];
                time_I = [time_I, str2double(tmp{2})];
            elseif strcmp(tmp{1}, '3')
                I = [I, I(end)];
                time_I = [time_I, str2double(tmp{2})];
            end
        end                    
    end
end        


