% get the individual alignments for the HA and NA segments from the gisaid
% data
clear
rng(1);
% split the alignment first into H3N2 and H1N1, theninto the all the segments HA, NA, PB1, PB2, PA, NP, MP, NS
c = ones(8,2); % counter
% vetor with the segment names
Segments = {'HA', 'NA', 'PB1', 'PB2', 'PA', 'NP', 'MP', 'NS'};
% define the two viruses
Viruses = {'H3N2', 'H1N1'};

% define which sequences to use as reference
reference = {'A/Alaska/232/2015', 'A/Arizona/33/2017'};

%define the max samples
target_sample_number = 600;

%% initialize data structure to store the alignments
Alignment = cell(size(c));
% check if data/H1N1_HA_aligned.fasta exists, if so, skip to the end
if ~exist('data/unaligned_H1N1_HA.fasta', 'file')
    % read in the fasta file in InfANorthAmerica-23052024.fasta
    fastaFile1 = fastaread('gisaid/InfANorthAmerica-23052024.fasta');
    fastaFile2 = fastaread('gisaid/InfANorthAmerica-202406051.fasta');
    fastaFile3 = fastaread('gisaid/InfANorthAmerica-202406052.fasta');
    fastaFile4 = fastaread('gisaid/InfANorthAmerica-23052029.fasta');
    % concatenate the two fasta files
    fastaFile = [fastaFile1; fastaFile2; fastaFile3; fastaFile4];
    
    % the headers are of form A/turkey/North_Dakota/22-012099-001-original/2022|EPI_ISL_16171250|2022-04-19|NA|original|A_/_H3N2 or H1N1

    for i = 1:length(fastaFile)
        header = strsplit(fastaFile(i).Header, '|', 'CollapseDelimiters', false);
        if length(header) == 6            
            % get the index of the segment
            segmentIdx = find(strcmp(Segments, header{4}));
            % do the same for virus, but only take the last 4 characters of header{6}
            if contains(header{6}, Viruses)
                virIdx = find(strcmp(Viruses, header{6}(end-3:end)));
                if ~isempty(virIdx)
                    % add the sequence to the alignment, then increase the counter for that alignment by 1
                    Alignment{segmentIdx, virIdx}(c(segmentIdx, virIdx)).Header = [header{1} '|' header{2} '|' header{3}];
                    Alignment{segmentIdx, virIdx}(c(segmentIdx, virIdx)).Sequence = fastaFile(i).Sequence;
                    c(segmentIdx, virIdx) = c(segmentIdx, virIdx) + 1;
                end
                if strcmp(header{1}, reference{virIdx})
                    % delete existing reference file
                    delete(['data/reference_' Viruses{virIdx} '_' Segments{segmentIdx} '.fasta']);
                    % write the reference sequence to a new file
                    ref = fastaFile(i);
                    ref(1).Header = 'Reference|Reference|1900-09-09'; % naming such that it drops out later by dates
                    fastawrite(['data/reference_' Viruses{virIdx} '_' Segments{segmentIdx} '.fasta'], ref);
                end
            else
                disp(header{6})
            end
        else
            error('Header does not have 6 fields');
        end
    end

    % write the alignments to the corresponding fasta files
    for i = 1:length(Viruses)
        for j = 1:length(Segments)
            delete(['data/unaligned_' Viruses{i} '_' Segments{j} '.fasta']);
            fastawrite(['data/unaligned_' Viruses{i} '_' Segments{j} '.fasta'], Alignment{j, i});
        end
    end
else
    disp('skip splitting fasta files')
end

if ~exist('data/aligned_H1N1_HA.fasta', 'file')
    % align the HA and NA sequences using mafft using the fastest possible settings in mafft
    % system('/opt/homebrew/bin/mafft --auto data/HA_unaligned.fasta > data/HA_aligned.fasta');
    for i = 1:length(Viruses)
        for j = 1:length(Segments)
            refFile = ['data/reference_' Viruses{i} '_' Segments{j} '.fasta'];
            unalignedFile = ['data/unaligned_' Viruses{i} '_' Segments{j} '.fasta'];
            alignedFile = ['data/aligned_' Viruses{i} '_' Segments{j} '.fasta'];
            if exist(refFile, 'file')
                system(['/opt/homebrew/bin/mafft --6merpair --thread -8 --keeplength --addfragments ' unalignedFile ' ' refFile ' > ' alignedFile]);                
            else
                error('reference not found')
            end
        end
    end
else
    disp('skip generating alignments')
end

%% clean the alignments by remove start and end, and identifying outlier sequences
% check if consensus exists
if ~exist('data/consensus_2010_H3N2_HA.fasta', 'file')
    for i = 1: length(Viruses)
        % Keep track of all the contextual sequences
        potential_sequences= cell(length(Segments),length(2010:2023));

        for s = 1:length(Segments)
            % read in the aligned sequences
            Sequences = fastaread(['data/aligned_' Viruses{i} '_' Segments{s} '.fasta']);
            delete(['data/alignedcleaned_' Viruses{i} '_' Segments{s} '.fasta'])

            % convert the sequences to a matrix
            seqMatrix = cell2mat({Sequences.Sequence}');
            % get the number of gaps per site in the alignment 
            numGaps = sum(seqMatrix == '-', 1);
            % define the first position as the first position where the number of gaps is less than 0.1 * the number of sequences and same for the last positions
            firstPos = find(numGaps < 0.05 * length(Sequences), 1, 'first')+50;
            lastPos = find(numGaps < 0.05 * length(Sequences), 1, 'last')-50;
            % get the sampling time for each sequence
            clear samplingTime
            for j = 1:length(Sequences)
                header = strsplit(Sequences(j).Header, '|', 'CollapseDelimiters', false);
                if length(header) == 3
                    samplingTime(j) = datetime(header{3}, 'InputFormat', 'yyyy-MM-dd');
                else
                    error('Header does not have 3 fields');
                end
            end
            % set up a vector that keeps track of which sequences will be discarded
            discard = true(length(Sequences), 1);

            % for every season, get a consensus sequence
            c=1;
            for year = 2010:2023
                first_isolation_time = datetime([num2str(year) '-09-01'], 'InputFormat', 'yyyy-MM-dd');
                last_isolation_time = datetime([num2str(year+1) '-05-01'], 'InputFormat', 'yyyy-MM-dd');
                isolates = find(samplingTime >= first_isolation_time & samplingTime < last_isolation_time);
                if length(isolates)>10
                    % get the consensus sequences over all sequences
                    consensus = mode(seqMatrix(isolates,firstPos:lastPos), 1);
                    % compute the distance of each sequence to the consensus
                    distances = sum(seqMatrix(isolates,firstPos:lastPos) ~= consensus, 2);
                    % get the sequence closest to the consensus
                    minDist = min(distances);
                    % save the one sequence closest to the consensus to a file consensus_year_segment.fasta
                    % remove all sequences that are more than 0.003*15 away on average from the consensus
                    totalLength = lastPos-firstPos;
                    maxDist = 0.003*15*totalLength;
                    % randomly pick 5 sequences with less than maxDist
                    use_isolates = isolates(distances < maxDist);
                    potential_sequences{s,c} = [potential_sequences{s,c}; Sequences(use_isolates)];
                    if year>=2015
                        % write to file all seqs that have distances <= maxDist
                        fastawrite(['data/alignedcleaned_' Viruses{i} '_' Segments{s} '.fasta'], Sequences(isolates(distances <= maxDist)));
                    end
                end
                c=c+1;
            end
        end
        c=1;
        for year = 2010:2023
            if length(potential_sequences{1,c})>0
                % for each year, pick 5 potential sequences that are in all segments
                potential_sequences_all = {potential_sequences{1,c}.Header};
                for s = 2:length(Segments)
                    potential_sequences_all = intersect(potential_sequences_all, {potential_sequences{s,c}.Header});
                end
                % randomly pick 5 sequences from potential_sequences_all
                if length(potential_sequences_all) > 5
                    seqs = randsample(potential_sequences_all, 5);
                    % for each Segment, print the seqs to file
                    for s = 1:length(Segments)
                        idx = find(ismember({potential_sequences{s,c}.Header}, seqs));
                        delete(['data/consensus_' num2str(year) '_' Viruses{i} '_' Segments{s} '.fasta'])
                        fastawrite(['data/consensus_' num2str(year) '_' Viruses{i} '_' Segments{s} '.fasta'], potential_sequences{s,c}(idx));
                    end
                end
            end
            c=c+1;
        end
    end
else
    disp('skip generating consensus sequences')
end


%% check if  H3N2_2016_HA.fasta exists
if ~exist('data/H3N2_2016_HA.fasta')

    % read in InfA_Cases.csv as a comma delimited file with header
    fid = fopen('InfA_Cases.csv');
    % loop over all lines
    c = 1;
    clear cases
    header = strsplit(fgets(fid), ',');
    while ~feof(fid)
        tline = fgets(fid);
        % split the line on ","
        split = strsplit(tline, ',');
        % get the date
        cases(c).date = datetime(split{11}, 'InputFormat', 'yyyy-MM-dd');
        % get the number of cases
        cases(c).h3 = str2double(split{4});
        cases(c).h1 = str2double(split{5});
        c=c+1;
    end

    % Make a vector with the population of each US state in the lower 48
    % based on 2020 census data
    states = {
        'Alabama', 5024279;
        'Arizona', 7151502;
        'Arkansas', 3011524;
        'California', 39538223;
        'Colorado', 5773714;
        'Connecticut', 3605944;
        'Delaware', 989948;
        'Florida', 21538187;
        'Georgia', 10711908;
        'Idaho', 1839106;
        'Illinois', 12812508;
        'Indiana', 6785528;
        'Iowa', 3190369;
        'Kansas', 2937880;
        'Kentucky', 4505836;
        'Louisiana', 4657757;
        'Maine', 1362359;
        'Maryland', 6177224;
        'Massachusetts', 7029917;
        'Michigan', 10077331;
        'Minnesota', 5706494;
        'Mississippi', 2961279;
        'Missouri', 6154913;
        'Montana', 1084225;
        'Nebraska', 1961504;
        'Nevada', 3104614;
        'New_Hampshire', 1377529;
        'New_Jersey', 9288994;
        'New_Mexico', 2117522;
        'New_York', 20201249;
        'North_Carolina', 10439388;
        'North_Dakota', 779094;
        'Ohio', 11799448;
        'Oklahoma', 3959353;
        'Oregon', 4237256;
        'Pennsylvania', 13002700;
        'Rhode_Island', 1097379;
        'South_Carolina', 5118425;
        'South_Dakota', 886667;
        'Tennessee', 6916897;
        'Texas', 29145505;
        'Utah', 3271616;
        'Vermont', 643077;
        'Virginia', 8631393;
        'Washington', 7693612;
        'West Virginia', 1793716;
        'Wisconsin', 5893718;
        'Wyoming', 576851;
        'District_Of_Columbia', 712816;
    };

    % compute the frequency of each state in the US
    state_freq = zeros(length(states), 1);
    for i = 1:length(states)
        state_freq(i) = states{i, 2};
    end
    state_freq = state_freq/sum(state_freq);

    % loop over the viruses
    for v = 1:length(Viruses)
        HA_aligned = fastaread(['data/alignedcleaned_' Viruses{v} '_HA.fasta']);
        NA_aligned = fastaread(['data/alignedcleaned_' Viruses{v} '_NA.fasta']);
        PB1_aligned = fastaread(['data/alignedcleaned_' Viruses{v} '_PB1.fasta']);
        PB2_aligned = fastaread(['data/alignedcleaned_' Viruses{v} '_PB2.fasta']);
        PA_aligned = fastaread(['data/alignedcleaned_' Viruses{v} '_PA.fasta']);
        NP_aligned = fastaread(['data/alignedcleaned_' Viruses{v} '_NP.fasta']);
        MP_aligned = fastaread(['data/alignedcleaned_' Viruses{v} '_MP.fasta']);
        NS_aligned = fastaread(['data/alignedcleaned_' Viruses{v} '_NS.fasta']);

        % from the HA segments, get the isolation times of each isolate
        clear HA_isolation_times names isolation_location names_part1
        for i = 1:length(HA_aligned)
            header = strsplit(HA_aligned(i).Header, '|', 'CollapseDelimiters', false);
            HA_isolation_times(i) = datetime(header{3}, 'InputFormat', 'yyyy-MM-dd');
            names{i} = HA_aligned(i).Header;
            names_part1{i} = header{1};
            isolation_location{i} = {};
            % get the state of each isolate
            for j = 1:length(states)
                if contains(header{1}, states{j, 1})
                    isolation_location{i} = states{j, 1};
                    break;
                end
            end
            % if isolation location for i not found, print the entire ehader
            if length(isolation_location)<i
                disp([header{1} ' ' header{2}   ' ' header{3}])
            end        
        end   

        % loop over the seasons 2010 until 2023
        for s = 2015:2023        
            % randomly select target_sample_number isolates from the HA segment, to ensure even sampling over time
            rng(1);
            % get the first and last isolation times, use only the month
            first_isolation_time = datetime([num2str(s) '-10-01'], 'InputFormat', 'yyyy-MM-dd');
            last_isolation_time = datetime([num2str(s+1) '-05-01'], 'InputFormat', 'yyyy-MM-dd');
            % num_months = calmonths(between(first_isolation_time, last_isolation_time));
            % get all isolates from that season
            isolates = find(HA_isolation_times >= first_isolation_time & HA_isolation_times < last_isolation_time);
            % check if the names_part1 of the isolates is unique, if there is one that is not unique, only keep the first one in isolates
            if length(unique(names_part1(isolates))) < length(isolates)
                [~, idx] = unique(names_part1(isolates));
                isolates = isolates(idx);
            end
            % get all the indices for cases between first and last isolation time
            cases_idx = find([cases.date] >= first_isolation_time & [cases.date] < last_isolation_time);
            % get the number of sequences for each state in isolates
            sample_number = zeros(length(states), 1);
            for i = 1:length(isolates)
                for j = 1:length(states)
                    if strcmp(isolation_location{isolates(i)}, states{j, 1})
                        sample_number(j) = sample_number(j) + 1;
                        break;
                    end
                end
            end
            % get the maximum sample number for each state allowed, but mulitplying
            % the state frequency by the number of min(400,length(isolates) by a factor of 5 and add 1
            max_samples = ceil(state_freq * min(target_sample_number,length(isolates)) * 5 + 1);            
            % while any entries in sample_number are larger than max_samples
            while any(sample_number > max_samples)
                sample_number = min(sample_number,max_samples);
                max_samples = ceil(state_freq * min(target_sample_number,sum(sample_number)) * 2 + 1);
            end

            if strcmp(Viruses{v}, 'H3N2')
                weekly_cases = [cases(cases_idx).h3];
            else
                weekly_cases = [cases(cases_idx).h1];
            end

            % define the most samples allowed per week based on the 
            max_samples_week = round(weekly_cases/sum(weekly_cases)*target_sample_number*1.25);
            cases_per_week = zeros(size(max_samples_week));

            % keep track of the sample IDs
            samples = cell(0, 0);
            last_samples = length(samples);
            same_sample_count = 0;
            samples_in_state = zeros(length(states), 1);
            while length(samples) < min(target_sample_number,sum(sample_number))
                % randomly add isolates to the samples
                if ~isempty(isolates)
                    % pick a week from where to sample based on the number of h3 or h1 cases
                    week_idx = randsample(length(cases_idx), 1, true, weekly_cases);

                    if (cases_per_week(week_idx)>= max_samples_week(week_idx))
                        continue;
                    end
                    
                    % pick a random isolate with sampling date between cases.date(week_idx) and cases.date(week_idx) + 7 days
                    potential_isolates = find(HA_isolation_times >= cases(week_idx).date & HA_isolation_times < cases(week_idx).date + 7);

                    if length(potential_isolates)==0
                        % skip skip the next steap
                        continue;
                    elseif length(potential_isolates) > 1
                        random_isolate = isolates(randi(length(isolates)));                
                    else
                        random_isolate = potential_isolates;
                    end
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
                            % get the sample location
                            for j = 1:length(states)
                                if strcmp(isolation_location{random_isolate}, states{j, 1})
                                    loc=i;
                                    if samples_in_state(j) < max_samples(j)
                                        samples_in_state(j) = samples_in_state(j) + 1;
                                        cases_per_week(week_idx) = cases_per_week(week_idx)+1;
                                        samples{end + 1} = HA_aligned(random_isolate).Header;
                                    end
                                    break;
                                end
                            end
                        end
                    end
                end
                if length(samples)==last_samples
                    same_sample_count=same_sample_count+1;
                else
                    same_sample_count=0;
                end
                last_samples = length(samples);
                if same_sample_count > 500
                    break;
                end
            end
            % add the first sample of the last season to the samples as outgroup
            % if s>2010
            %     lastseasonseq = fastaread(['data/' Viruses{v} '_' num2str(s-1) '_HA.fasta']);
            %     % get the first sampled sequence in time from the last season.
            %     clear times
            %     for i = 1:length(lastseasonseq)
            %         header = strsplit(lastseasonseq(i).Header, '|', 'CollapseDelimiters', false);
            %         times(i) = datetime(header{3}, 'InputFormat', 'yyyy-MM-dd');
            %     end
            %     [~, idx] = min(times);
            %     samples{end + 1} = lastseasonseq(idx).Header;
            % end

            % delet the existing fasta files
            if exist(['data/' Viruses{v} '_' num2str(s) '_HA.fasta'], 'file')
                delete(['data/' Viruses{v} '_' num2str(s) '_HA.fasta']);
                delete(['data/' Viruses{v} '_' num2str(s) '_NA.fasta']);
                delete(['data/' Viruses{v} '_' num2str(s) '_PB1.fasta']);
                delete(['data/' Viruses{v} '_' num2str(s) '_PB2.fasta']);
                delete(['data/' Viruses{v} '_' num2str(s) '_PA.fasta']);
                delete(['data/' Viruses{v} '_' num2str(s) '_NP.fasta']);
                delete(['data/' Viruses{v} '_' num2str(s) '_MP.fasta']);
                delete(['data/' Viruses{v} '_' num2str(s) '_NS.fasta']);
            end

            for i = 1:length(samples)
                HA_index = find(strcmp({HA_aligned.Header}, samples{i}));
                NA_index = find(strcmp({NA_aligned.Header}, samples{i}));
                PB1_index = find(strcmp({PB1_aligned.Header}, samples{i}));
                PB2_index = find(strcmp({PB2_aligned.Header}, samples{i}));
                PA_index = find(strcmp({PA_aligned.Header}, samples{i}));
                NP_index = find(strcmp({NP_aligned.Header}, samples{i}));
                MP_index = find(strcmp({MP_aligned.Header}, samples{i}));
                NS_index = find(strcmp({NS_aligned.Header}, samples{i}));
                fastawrite(['data/' Viruses{v} '_' num2str(s) '_HA.fasta'], HA_aligned(HA_index(1)));
                fastawrite(['data/' Viruses{v} '_' num2str(s) '_NA.fasta'], NA_aligned(NA_index(1)));
                fastawrite(['data/' Viruses{v} '_' num2str(s) '_PB1.fasta'], PB1_aligned(PB1_index(1)));
                fastawrite(['data/' Viruses{v} '_' num2str(s) '_PB2.fasta'], PB2_aligned(PB2_index(1)));
                fastawrite(['data/' Viruses{v} '_' num2str(s) '_PA.fasta'], PA_aligned(PA_index(1)));
                fastawrite(['data/' Viruses{v} '_' num2str(s) '_NP.fasta'], NP_aligned(NP_index(1)));
                fastawrite(['data/' Viruses{v} '_' num2str(s) '_MP.fasta'], MP_aligned(MP_index(1)));
                fastawrite(['data/' Viruses{v} '_' num2str(s) '_NS.fasta'], NS_aligned(NS_index(1)));
            end
        end
    end
else
    disp('skip sampling sequences')
end

% in each dataset, check for outliers using raxml and tree time
% system('rm -r raxml');
% system('mkdir raxml');
for v = 1:length(Viruses)
    for y = 2015:2023
        if ~exist(['raxml/' Viruses{v} '_' num2str(y) '_NS.treefile'], 'file') 
            for s = 1:length(Segments)          
                % build a new fasta file for raxml that contains the fasta for for this season, virus and segment+the consensus sequences from the previous 5 seasons if available
                % read in the aligned sequences
                if ~exist(['data/' Viruses{v} '_' num2str(y) '_' Segments{s} '.fasta'], 'file')
                    continue;
                end
                Sequences = fastaread(['data/' Viruses{v} '_' num2str(y) '_' Segments{s} '.fasta']);
                % get the sequences fromt he prev consensus
                for year = 2010:2023
                    if exist(['data/consensus_' num2str(year) '_' Viruses{v} '_' Segments{s} '.fasta'], 'file') && year~=y
                        consensus = fastaread(['data/consensus_' num2str(year) '_' Viruses{v} '_' Segments{s} '.fasta']);
                        % check that the consensus sequences are unique, if not, take the first one
                        [~, idx] = unique({consensus.Header});
                        Sequences = [Sequences; consensus(idx)];
                    end
                end

                % write the sequences to a new filed
                delete(['raxml/' Viruses{v} '_' num2str(y) '_' Segments{s} '_raxml.fasta'])
                fastawrite(['raxml/' Viruses{v} '_' num2str(y) '_' Segments{s} '_raxml.fasta'], Sequences);
                system(['/opt/homebrew/bin/iqtree2 -nt 11 -s  raxml/'  Viruses{v} '_' num2str(y) '_' Segments{s} '_raxml.fasta -m GTR --prefix raxml/'  Viruses{v} '_' num2str(y) '_' Segments{s} '']);
            end
        end
    end
end


for v = 1:length(Viruses)
    for y = 2015:2023
        if ~exist(['timetree/' Viruses{v} '_' num2str(y) '_NS_timetree.nexus'], 'file') 
            for s = 1:length(Segments)        
                if ~exist(['data/' Viruses{v} '_' num2str(y) '_' Segments{s} '.fasta'], 'file')
                    continue;
                end

                Sequences = fastaread(['raxml/' Viruses{v} '_' num2str(y) '_' Segments{s} '_raxml.fasta']);
                % make a new file with the names and dates
                delete(['raxml/' Viruses{v} '_' num2str(y) '_' Segments{s} '_dates.csv'])
                fid = fopen(['raxml/' Viruses{v} '_' num2str(y) '_' Segments{s} '_dates.csv'], 'w');
                fprintf(fid, 'name,date\n');
                for i=1:length(Sequences)
                    header = strsplit(Sequences(i).Header, '|', 'CollapseDelimiters', false);
                    fprintf(fid, '%s,%s\n', Sequences(i).Header, header{3});
                end
    
                delete('timetree/outliers.tsv')
                delete('timetree/timetree.nexus');
                delete('timetree/molecular_clock.txt');
                delete('timetree/root_to_tip_regression.pdf');
                delete('timetree/ancestral_sequences.fasta');
                delete('timetree/auspice_tree.json')
                delete('timetree/divergence_tree.nexus')
                delete('timetree/branch_mutations.txt')
                delete('timetree/sequence_evolution_model.txt')
                delete('timetree/timetree.pdf')
                delete('timetree/trace_run.log')

                system(['/opt/homebrew/Caskroom/miniconda/base/envs/plasmids/bin/python -m treetime --tree raxml/'  Viruses{v} '_' num2str(y) '_' Segments{s}...
                    '.treefile --aln raxml/' Viruses{v} '_' num2str(y) '_' Segments{s} '_raxml.fasta'...
                    ' --dates ' 'raxml/' Viruses{v} '_' num2str(y) '_' Segments{s} '_dates.csv'...
                    ' --outdir timetree/ ']);    

                % rename the file divergence_tree.nexus to the name of the alignment.nexus
                movefile('timetree/timetree.nexus', ['timetree/' Viruses{v} '_' num2str(y) '_' Segments{s} '_timetree.nexus']);
                movefile('timetree/molecular_clock.txt', ['timetree/' Viruses{v} '_' num2str(y) '_' Segments{s} '_molecular_clock.txt']);
                movefile('timetree/root_to_tip_regression.pdf', ['timetree/' Viruses{v} '_' num2str(y) '_' Segments{s} '_root_to_tip_regression.pdf']);
                movefile('timetree/ancestral_sequences.fasta', ['timetree/' Viruses{v} '_' num2str(y) '_' Segments{s} '_ancestral.fasta']);
                if exist('timetree/outliers.tsv', 'file')
                    movefile('timetree/outliers.tsv', ['timetree/' Viruses{v} '_' num2str(y) '_' Segments{s} '_outliers.tsv']);
                end
            end
        end
    end
end

useTreeTimeOutliers = true;

if useTreeTimeOutliers
    for v = 1:length(Viruses)
        for y = 2015:2023
            if ~exist(['timetree/' Viruses{v} '_' num2str(y) '_' Segments{1} '_ancestral.fasta'], 'file')
                continue;
            end
            % get all files in the folder timetree that start with Viruses{v}_year{v}_*_outliers.tsv
            files = dir(['timetree/' Viruses{v} '_' num2str(y) '_*_outliers.tsv']);
            % start a new file raxml/outliers_ Viruses{v} '_' num2str(y) '.csv'
            foutlier = fopen(['timetree/outliers_' Viruses{v} '_' num2str(y) '.csv'], 'w');
            for i = 1:length(files)
                % read in the outliers file
                outliers = readtable(['timetree/' files(i).name], 'FileType', 'text', 'ReadVariableNames', false);
                for j = 1:height(outliers)
                    fprintf(foutlier, '%s\n', outliers{j, 1}{1});
                end
            end
            fclose("all");

        end
    end
else
    for v = 1:length(Viruses)
        for y = 2015:2023
            if ~exist(['timetree/' Viruses{v} '_' num2str(y) '_' Segments{1} '_ancestral.fasta'], 'file')
                continue;
            % elseif exist(['timetree/outliers' Viruses{v} '_' num2str(y) '_outliers.png'], 'file')
            %     continue;
            end
    
            % start a new file raxml/outliers_ Viruses{v} '_' num2str(y) '.csv'
            foutlier = fopen(['timetree/outliers_' Viruses{v} '_' num2str(y) '.csv'], 'w');
    
            % make a 4 by 2 subplot
            figure
    
    
            for s = 1:length(Segments)        
                % read in ancestral_sequences
                Sequences = fastaread(['timetree/' Viruses{v} '_' num2str(y) '_' Segments{s} '_ancestral.fasta']);
                % compute the distance of every node that doesn't start ith NODE, relative to NODE_0000001
                % get the sequence of NODE_0000001
                refSeq='';
                % for i = 1 :length(Sequences)
                %     if startsWith(Sequences(i).Header, 'NODE_0000001')
                        refSeq = Sequences(1).Sequence;
                %     end
                % end
                
                clear time dist
                % compute the distance of every node to the reference sequence
                for i = 1:length(Sequences)
                    if ~startsWith(Sequences(i).Header, 'NODE')
                        % compute the distance for every base not N
                        dist(i) = sum(Sequences(i).Sequence ~= refSeq & Sequences(i).Sequence ~= 'N')/length(Sequences(i).Sequence);
                        tmp = strsplit(Sequences(i).Header, '|');                        
                        time(i) = datetime(tmp{3}, 'InputFormat', 'yyyy-MM-dd');
                    end
                end
                time_d = datenum(time(:));
                isnanvals = isnan(time_d); 
                % plot the distance over time
                subplot(4,2,s)
                plot(time_d, dist, 'o');hold on
                time_n = time_d(~isnanvals);
                dist_n = dist(~isnanvals);
                header_n = {Sequences(~isnanvals).Header}; 
    
                % sort the time and get the index
                [time_n, idx] = sort(time_n);
                dist_n = dist_n(idx);
                header_n = header_n(idx);
                % Fit a smooth curve using 'fit' function (e.g., using a smoothing spline)
                fitObj = fit(time_n, dist_n', 'poly1');                
                % Plot the smooth curve
                smoothCurve = feval(fitObj, time_n);
                plot(time_n, smoothCurve, '-r', 'LineWidth', 2);
                
                % Calculate residuals (difference between actual data and smooth curve)
                residuals = dist_n' - smoothCurve;                
                % Detect massive outliers
                threshold = 7 * std(residuals);
                outliers = abs(residuals) > threshold;
                
                no_outlier_idx = find(~outliers);

                % Fit a smooth curve using 'fit' function (e.g., using a smoothing spline)
                fitObj = fit(time_n(no_outlier_idx), dist_n(no_outlier_idx)', 'poly1');                
                % Plot the smooth curve
                smoothCurve = feval(fitObj, time_n(no_outlier_idx));
                plot(time_n(no_outlier_idx), smoothCurve, '-b', 'LineWidth', 2);
                
                % Calculate residuals (difference between actual data and smooth curve)
                residuals = dist_n(no_outlier_idx)' - smoothCurve;                
                % Detect massive outliers
                threshold = 5 * std(residuals);
                norm_outliers = abs(residuals) > threshold;
                outliers(no_outlier_idx(norm_outliers)) = true;

                % Plot outliers
                plot(time_n(outliers), dist_n(outliers), 'ro', 'MarkerFaceColor', 'r');


                
                % Add labels and legend
                xlabel('Time');
                ylabel('Distance');
                legend('Data', 'Smooth Curve', 'Outliers');
                title([Segments{s}]);
                % remove legend
                legend('off');
    
                % print the outliers to file
                for i = 1:length(outliers)
                    if outliers(i)>0
                        fprintf(foutlier, '%s\n', header_n{i});
                    end
                end
            end
            % save to file
            saveas(gcf, ['timetree/outliers' Viruses{v} '_' num2str(y) '_outliers.png']);
            close(gcf);

            if y==2022
                dsf
            end
                
        end
    end 
end

system('rm -r xmls');
system('mkdir xmls');
for v = 1:length(Viruses)
    for y = 2015:2023
        if exist(['timetree/outliers_' Viruses{v} '_' num2str(y) '.csv'], 'file')
            % get all sequences in the outliers fike
            f = fopen(['timetree/outliers_' Viruses{v} '_' num2str(y) '.csv']);
            outliers = cell(0,0);
            while ~feof(f)
                line = fgets(f);
                if line~=-1
                    tmp = strsplit(line);
                    outliers{end+1} = tmp{1};
                end
            end
            fclose(f);
            for s = 1:length(Segments)          
                % build a new fasta file for raxml that contains the fasta for for this season, virus and segment+the consensus sequences from the previous 5 seasons if available
                % read in the aligned sequences
                if ~exist(['data/' Viruses{v} '_' num2str(y) '_' Segments{s} '.fasta'], 'file')
                    continue;
                end
                Sequences = fastaread(['data/' Viruses{v} '_' num2str(y) '_' Segments{s} '.fasta']);

                % for yy = 1:5;
                %     year=y-yy
                %     if exist(['data/consensus_' num2str(year) '_' Viruses{v} '_' Segments{s} '.fasta'], 'file') && year~=y
                %         consensus = fastaread(['data/consensus_' num2str(year) '_' Viruses{v} '_' Segments{s} '.fasta']);
                %         % check that the consensus sequences are unique, if not, take the first one
                %         [~, idx] = unique({consensus.Header});
                %         Sequences = [Sequences; consensus(idx)];
                %     end
                % end


                % remove all sequences in outliers
                for i = length(Sequences):-1:1
                    if any(strcmp(Sequences(i).Header, outliers))
                        Sequences(i) = [];
                    end
                end
                
                % if s==1
                %     use_headers = {Sequences(randsample(length(Sequences), min(length(Sequences),400))).Header};
                % end
                % for i = length(Sequences):-1:1
                %     if ~any(strcmp(Sequences(i).Header, use_headers))
                %         Sequences(i) = [];
                %     end
                % end
                % 
                % % get the genomic distance of every sequence to the first sampled sequence, then plot it vs. time sampling time
                % % get the sampling times
                % times = zeros(length(Sequences), 1);
                % for i = 1 :length(Sequences)
                %     header = strsplit(Sequences(i).Header, '|', 'CollapseDelimiters', false);
                %     times(i) = datenum(datetime(header{3}, 'InputFormat', 'yyyy-MM-dd'));
                % end
                % % get the distance of every sequence to the first sequence
                % miseq = min(times);
                % for i = 1:length(Sequences)
                %     if times(i) == miseq
                %         refSeq = Sequences(i).Sequence;
                %     end
                % end
                % 
                % % get the distance of every sequence to the reference sequence
                % for i = 1:length(Sequences)
                %     dist(i) = sum(Sequences(i).Sequence ~= refSeq & Sequences(i).Sequence ~= 'N')/length(Sequences(i).Sequence);
                % end
                % % plot the distance over time
                % close all
                % plot(times, dist, 'o');hold on
                % % add title
                % title([Viruses{v} ' ' num2str(y) ' ' Segments{s}]);
                % pause(1)
                % 

                
                % write the sequences to a new filed
                delete(['xmls/' Viruses{v} '_' num2str(y) '_' Segments{s} '.fasta'])
                fastawrite(['xmls/' Viruses{v} '_' num2str(y) '_' Segments{s} '.fasta'], Sequences);
            end
        end
    end
end




