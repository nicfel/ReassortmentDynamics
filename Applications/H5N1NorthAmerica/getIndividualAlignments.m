% get the individual alignments for the HA and NA segments from the gisaid
% data
clear
Segments = {'HA', 'NA', 'PB1', 'PB2', 'PA', 'NP', 'MP', 'NS'};
% read in the fasta file in data/H5N1northAmerica-07052024.fasta
fastaFile = fastaread('data/H5N1northAmerica-07052024.fasta');
% the headers are of form A/turkey/North_Dakota/22-012099-001-original/2022|EPI_ISL_16171250|2022-04-19|NA|original|A_/_H5N1
% split the alignment into the all the segments HA, NA, PB1, PB2, PA, NP, Mp, NS
c = ones(8,1); % counter
% check if data/HA_aligned.fasta exists, if so, skip to the end
if ~exist('data/HA_aligned.fasta', 'file')
    for i = 1:length(fastaFile)
        header = strsplit(fastaFile(i).Header, '|', 'CollapseDelimiters', false);
        if length(header) == 6
            if strcmp(header{4}, 'HA')
                HA(c(1)).Header = [header{1} '|' header{2} '|' header{3}];
                HA(c(1)).Sequence = fastaFile(i).Sequence;
                c(1) = c(1) + 1;
            elseif strcmp(header{4}, 'NA')
                NA(c(2)).Header = [header{1} '|' header{2} '|' header{3}];
                NA(c(2)).Sequence = fastaFile(i).Sequence;
                c(2) = c(2) + 1;
            elseif strcmp(header{4}, 'PB1')
                PB1(c(3)).Header = [header{1} '|' header{2} '|' header{3}];
                PB1(c(3)).Sequence = fastaFile(i).Sequence;
                c(3) = c(3) + 1;
            elseif strcmp(header{4}, 'PB2')
                PB2(c(4)).Header = [header{1} '|' header{2} '|' header{3}];
                PB2(c(4)).Sequence = fastaFile(i).Sequence;
                c(4) = c(4) + 1;
            elseif strcmp(header{4}, 'PA')
                PA(c(5)).Header = [header{1} '|' header{2} '|' header{3}];
                PA(c(5)).Sequence = fastaFile(i).Sequence;
                c(5) = c(5) + 1;
            elseif strcmp(header{4}, 'NP')
                NP(c(6)).Header = [header{1} '|' header{2} '|' header{3}];
                NP(c(6)).Sequence = fastaFile(i).Sequence;
                c(6) = c(6) + 1;
            elseif strcmp(header{4}, 'MP')
                MP(c(7)).Header = [header{1} '|' header{2} '|' header{3}];
                MP(c(7)).Sequence = fastaFile(i).Sequence;
                c(7) = c(7) + 1;
            elseif strcmp(header{4}, 'NS')
                NS(c(8)).Header = [header{1} '|' header{2} '|' header{3}];
                NS(c(8)).Sequence = fastaFile(i).Sequence;
                c(8) = c(8) + 1;
            else
                error('Unknown segment');
            end
        else
            error('Header does not have 6 fields');
        end
    end

    % write the HA and NA sequences to fasta files after deleting the existing ones
    if exist('data/HA_unaligned.fasta', 'file')
        delete('data/HA_unaligned.fasta');
        delete('data/NA_unaligned.fasta');
        delete('data/PB1_unaligned.fasta');
        delete('data/PB2_unaligned.fasta');
        delete('data/PA_unaligned.fasta');
        delete('data/NP_unaligned.fasta');
        delete('data/MP_unaligned.fasta');
        delete('data/NS_unaligned.fasta');        
    end

    fastawrite('data/HA_unaligned.fasta', HA);
    fastawrite('data/NA_unaligned.fasta', NA);
    fastawrite('data/PB1_unaligned.fasta', PB1);
    fastawrite('data/PB2_unaligned.fasta', PB2);
    fastawrite('data/PA_unaligned.fasta', PA);
    fastawrite('data/NP_unaligned.fasta', NP);
    fastawrite('data/MP_unaligned.fasta', MP);
    fastawrite('data/NS_unaligned.fasta', NS);


    % align the HA and NA sequences using mafft using the fastest possible settings in mafft
    system('/opt/homebrew/bin/mafft --auto data/HA_unaligned.fasta > data/HA_aligned.fasta');
    system('/opt/homebrew/bin/mafft --auto data/NA_unaligned.fasta > data/NA_aligned.fasta');
    system('/opt/homebrew/bin/mafft --auto data/PB1_unaligned.fasta > data/PB1_aligned.fasta');
    system('/opt/homebrew/bin/mafft --auto data/PB2_unaligned.fasta > data/PB2_aligned.fasta');
    system('/opt/homebrew/bin/mafft --auto data/PA_unaligned.fasta > data/PA_aligned.fasta');
    system('/opt/homebrew/bin/mafft --auto data/NP_unaligned.fasta > data/NP_aligned.fasta');
    system('/opt/homebrew/bin/mafft --auto data/MP_unaligned.fasta > data/MP_aligned.fasta');
    system('/opt/homebrew/bin/mafft --auto data/NS_unaligned.fasta > data/NS_aligned.fasta');    
else
    disp('skip generating alignments')
end

% in each dataset, check for outliers using raxml and tree time
% system('rm -r raxml');
% system('mkdir raxml');
if ~exist(['raxml/NS_timetree.nexus'], 'file') 
    % start a new file raxml/outliers_ Viruses{v} '_' num2str(y) '.csv'
    for s = 1:length(Segments)          
        % build a new fasta file for raxml that contains the fasta for for this season, virus and segment+the consensus sequences from the previous 5 seasons if available
        % read in the aligned sequences
        Sequences = fastaread(['data/' Segments{s} '_aligned.fasta']);
        % clean heade names by remoinv '
        for i =1:length(Sequences)
            Sequences(i).Header = strrep(Sequences(i).Header, '''','');
            Sequences(i).Header = strrep(Sequences(i).Header, '(','');
            Sequences(i).Header = strrep(Sequences(i).Header, ')','');
        end
        % check for dublicates, and only keep the first
        [~, idx] = unique({Sequences.Header});
        Sequences = Sequences(idx);

        % write the sequences to a new filed
        delete(['raxml/' Segments{s} '_raxml.fasta'])
        fastawrite(['raxml/' Segments{s} '_raxml.fasta'], Sequences);
        system(['/opt/homebrew/bin/iqtree2-nt 11 -s raxml/'  Segments{s} '_raxml.fasta -m GTR --prefix raxml/'  Segments{s} '']);

        % make a new file with the names and dates
        delete(['timetree/'  Segments{s} '_dates.csv'])
        fid = fopen(['timetree/' Segments{s} '_dates.csv'], 'w');
        fprintf(fid, 'name,date\n');
        for i=1:length(Sequences)
            header = strsplit(Sequences(i).Header, '|', 'CollapseDelimiters', false);
            fprintf(fid, '%s,%s\n', Sequences(i).Header, header{3});
        end

    end
end

foutlier = fopen(['raxml/H5N1_outliers.csv'], 'w');
for s = 1:length(Segments)      
    if ~exist(['raxml/NS_timetree.nexus'], 'file') 
        delete('raxml/outliers.tsv')
        system(['/opt/homebrew/Caskroom/miniconda/base/envs/plasmids/bin/python -m treetime --tree raxml/'  Segments{s}...
            '.treefile --aln raxml/' Segments{s} '_raxml.fasta'...
            ' --dates ' 'timetree/' Segments{s} '_dates.csv'...
            ' --outdir timetree/ ']);    
        % read in the raxml/outliers/tsv file, and attach it to the file raxml/outliers_ Viruses{v} '_' num2str(y) '.csv' files after skipping the first line
        if exist('raxml/outliers.tsv', 'file')
            fid = fopen('raxml/outliers.tsv');
            fgetl(fid);
            while ~feof(fid)
                line = fgetl(fid);
                fprintf(foutlier, '%s\t%s\n', line, Segments{s});
            end
        end
        % rename the file divergence_tree.nexus to the name of the alignment.nexus
        movefile('raxml/timetree.nexus', ['raxml/' Segments{s} '_timetree.nexus']);
        movefile('raxml/molecular_clock.txt', ['raxml/' Segments{s} '_molecular_clock.txt']);
        movefile('raxml/root_to_tip_regression.pdf', ['raxml/' Segments{s} '_root_to_tip_regression.pdf']);

    end
end

% read in the aligned HA and NA sequences
HA_aligned = fastaread('raxml/HA_raxml.fasta');
NA_aligned = fastaread('raxml/NA_raxml.fasta');
PB1_aligned = fastaread('raxml/PB1_raxml.fasta');
PB2_aligned = fastaread('raxml/PB2_raxml.fasta');
PA_aligned = fastaread('raxml/PA_raxml.fasta');
NP_aligned = fastaread('raxml/NP_raxml.fasta');
MP_aligned = fastaread('raxml/MP_raxml.fasta');
NS_aligned = fastaread('raxml/NS_raxml.fasta');


% from the HA segments, get the isolation times of each isolate
for i = 1:length(HA_aligned)
    header = strsplit(HA_aligned(i).Header, '|', 'CollapseDelimiters', false);
    HA_isolation_times(i) = datetime(header{3}, 'InputFormat', 'yyyy-MM-dd');
end

% randomly select 200 isolates from the HA segment, to ensure even sampling over time
rng(1);
% get the first and last isolation times, use only the month
first_isolation_time = datetime('2021-01-01', 'InputFormat', 'yyyy-MM-dd');
last_isolation_time = dateshift(max(HA_isolation_times), 'end', 'month');
num_months = calmonths(between(first_isolation_time, last_isolation_time));
dairy_samples = 0;

% read in the outliers from H5N1_outliers.csv
f = fopen('raxml/H5N1_outliers.csv');
outliers = cell(0,0);
while ~feof(f)
    tmp = strsplit(fgets(f));
    outliers{end+1, 1} = tmp{1};
end
system('rm -r xmls')
system('mkdir xmls')

% keep track of the sample IDs
samples = cell(0, 0);
for timeAveraged = {'averaged', 'proportional'}
    while length(samples) < 600
        if strcmp(timeAveraged, 'averaged')
            % Ramdomly pick a month between the first and last isolation times
            random_month_offset = randi([0, num_months]);
            random_isolation_time = dateshift(first_isolation_time, 'start', 'month', random_month_offset);    
            % Get all isolates sampled in this month
            isolates = find(HA_isolation_times >= random_isolation_time & HA_isolation_times < dateshift(random_isolation_time, 'start', 'month', 1));
        else
            isolates = find(HA_isolation_times >= first_isolation_time & HA_isolation_times <= last_isolation_time);
        end
    
        % randomly select one isolate from this month
        if ~isempty(isolates)
            random_isolate = isolates(randi(length(isolates)));
            % check that the random isolate was not already sampled
            if ~any(strcmp(samples, HA_aligned(random_isolate).Header))
                % ensure that the isolate is present in all other alignments as well
                if any(strcmp({NA_aligned.Header}, HA_aligned(random_isolate).Header)) && ...
                   any(strcmp({PB1_aligned.Header}, HA_aligned(random_isolate).Header)) && ...
                   any(strcmp({PB2_aligned.Header}, HA_aligned(random_isolate).Header)) && ...
                   any(strcmp({PA_aligned.Header}, HA_aligned(random_isolate).Header)) && ...
                   any(strcmp({NP_aligned.Header}, HA_aligned(random_isolate).Header)) && ...
                   any(strcmp({MP_aligned.Header}, HA_aligned(random_isolate).Header)) && ...
                   any(strcmp({NS_aligned.Header}, HA_aligned(random_isolate).Header))
                    if ~any(strcmp(outliers, HA_aligned(random_isolate).Header))
                        if contains(HA_aligned(random_isolate).Header, 'dairy')
                            if dairy_samples <= 3
                                disp(dairy_samples)
                                samples{end + 1} = HA_aligned(random_isolate).Header;
                            end
                            dairy_samples = dairy_samples + 1;
                        else
                            samples{end + 1} = HA_aligned(random_isolate).Header;
                        end
                    end
                end
            end
        end
    end
    
    % get the sequences of the selected isolates for both HA and NA and teh save them to HA.fasta and NA.fasta after deleting the existing ones
    for i = 1:length(samples)
        HA_index = find(strcmp({HA_aligned.Header}, samples{i}));
        NA_index = find(strcmp({NA_aligned.Header}, samples{i}));
        PB1_index = find(strcmp({PB1_aligned.Header}, samples{i}));
        PB2_index = find(strcmp({PB2_aligned.Header}, samples{i}));
        PA_index = find(strcmp({PA_aligned.Header}, samples{i}));
        NP_index = find(strcmp({NP_aligned.Header}, samples{i}));
        MP_index = find(strcmp({MP_aligned.Header}, samples{i}));
        NS_index = find(strcmp({NS_aligned.Header}, samples{i}));
        fastawrite(['xmls/' timeAveraged{1} '_HA.fasta'], HA_aligned(HA_index));
        fastawrite(['xmls/' timeAveraged{1} '_NA.fasta'], NA_aligned(NA_index));
        fastawrite(['xmls/' timeAveraged{1} '_PB1.fasta'], PB1_aligned(PB1_index));
        fastawrite(['xmls/' timeAveraged{1} '_PB2.fasta'], PB2_aligned(PB2_index));
        fastawrite(['xmls/' timeAveraged{1} '_PA.fasta'], PA_aligned(PA_index));
        fastawrite(['xmls/' timeAveraged{1} '_NP.fasta'], NP_aligned(NP_index));
        fastawrite(['xmls/' timeAveraged{1} '_MP.fasta'], MP_aligned(MP_index));
        fastawrite(['xmls/' timeAveraged{1} '_NS.fasta'], NS_aligned(NS_index));        
    end
end




