%Data Preprocessing: for spike-phase modeling in python
%Goal: grab interneuron spikes and rhythmic phases during odor sampling


%this script will organize the data into .csv files 
%columns:
%rat_id
%session_id
%odorblock/condition
%trial_labels
%trial_segment (pre/during/post)
%time (msec)
%lfp id #1
%...
%lfp id #2
%filtered/phase lfp #1 - theta
%filtered/phase lfp #1 - beta
%filtered/phase lfp #1 - lowgamma
%filtered/phase lfp #1 - highgamma
%...
%filtered/phase lfp #2 - theta
%filtered/phase lfp #2 - beta
%filtered/phase lfp #2 - lowgamma
%filtered/phase lfp #2 - highgamma
%neuron id #1
%...
%neuron id #2

clear 
close all
dbstop if error

addpath data\
addpath lfp_data\
addpath python_spkphase_odorsamp\
addpath python_spkphase_odorsamp_dg\
addpath python_spkphase_odorsamp_pyrs\
addpath python_spkphase_odorsamp_nel\
addpath python_spkphase_odorsamp_pyrs_nel\

% savefileto = 'python_spkphase_odorsamp_pyrs\';
% savefileto = 'python_spkphase_odorsamp\';
savefileto = 'python_spkphase_odorsamp_nel\';
% load('Interneuron_Identities_Automaze.mat')
load('tmp_identities.mat') %contains list that I believe are unique ints
% tab = readtable('Int_Identities_DG.csv'); %this doesn't work yet, but just 
                                          %run select_neurons.m and go
                                          %straight to making the table a
                                          %cell array
% load('CA1_Pyr_Identities.mat')


%decide ahead of time the time interval edges surrounding the event ts
time_surr_event_pre = [-0.25 0]; %add negative sign for pre event
time_surr_event_dur = [0 1.5];
time_surr_event_post = [1.5 1.75];

Int_Identities = Unique_Int_Identities; %only if you load 'tmp_identities.mat', 
                                        % which Pam created to avoid possible double-counts
% Int_Identities = Pyr_Identities_CA1;


for i = 1:length(Int_Identities)
    sprintf('on %d out of %d',i,length(Int_Identities))

    %initialize the variables that will become table columns
    rat_names = [];
    session_labels = [];
    odor_block_labels = [];
    trial_labels = [];
    trial_segment = [];
    accuracy = [];
    quarter_labels = [];

    odor_labels = [];
    pos_labels = [];

    time = [];
    lfp_odorsamp = [];

    filtered_theta = [];
    filtered_beta = [];
    filtered_lowgamma = [];
    filtered_highgamma = [];

    unit_spikes = [];
    
    
    gather_len_trials = []; %this serves as a nice check
    
    %decompose the first column of Int_Identities
    rat_expression = 'LH\d{1,2}';
    [rat_id] = regexp(Int_Identities(i,1),rat_expression,'match');
    sess_expression = 'D\d{1,2}';
    session_id = regexp(Int_Identities(i,1),sess_expression,'match');

    %get the struct containing this rat/session's metadata
    filestruct = find_files(Rat,rat_id,session_id);

    %to make string comparisons easier later (e.g. to find existing files
    %with these names, etc) make all the names the same length. this will
    %involve adding 0's to rats whose numbers are less than 10, and similar
    %for session days
    rat_id = char(rat_id{1});
    if numel(rat_id)<4
        %insert a zero between LH and the single digit
        exp = 'LH';
        insert = '0';
        rat_id = [exp,insert,rat_id(end)];
    end
    
    session_id = char(session_id{1});
    if numel(session_id)<3
        %somehow insert a zero between D and the single digit
        exp = 'D';
        insert = '0';
        session_id = [exp,insert,session_id(end)];
    end

    %============================================================
    %                          grab event_ts
    %============================================================
    
    [Intervals,Intervals_by_odorblock] = find_event_ts_automaze(filestruct,'data\');

    %now grab the "accuracy_byquarter"; older rats have var
    %"accuracy_byodor" but looks like it is the same thing
    if isfield(filestruct, 'accuracy_byquarter')
        
        acc_list = filestruct.accuracy_byquarter;

    elseif isfield(filestruct,'accuracy_byodor')

        acc_list = filestruct.accuracy_byquarter;

    end

    %TODO: if either accuracy var above is empty, you should really skip this 
    %session--Lara must have prevented the variable from filling in
    %for a reason!!

    %all block 1
    event_ts_block01 = Intervals_by_odorblock.odorblock1.correct(:,1);
    rm_ind = event_ts_block01 == 0;
    event_ts_block01(rm_ind)=[];
    
    odors_block01 = Intervals_by_odorblock.odorblock1.odor_id;
    odors_block01(rm_ind) = []; %to match the correct trials

    pos_block01 = Intervals_by_odorblock.odorblock1.pos_id;
    pos_block01(rm_ind) = [];

    %all block 2
    event_ts_block02 = Intervals_by_odorblock.odorblock2.correct(:,1);
    rm_ind = event_ts_block02==0;
    event_ts_block02(rm_ind)=[];

    odors_block02 = Intervals_by_odorblock.odorblock2.odor_id;
    odors_block02(rm_ind) = [];

    pos_block02 = Intervals_by_odorblock.odorblock2.pos_id;
    pos_block02(rm_ind) = [];
    
    %this will help you iterate by block later
    event_ts = {event_ts_block01,event_ts_block02};
    odorblock = {odors_block01,odors_block02};
    posblock = {pos_block01,pos_block02};


    %============================================================
    %                          grab lfp
    %============================================================
    %find the unit name first, so that you can figure out which wires you
    %can grab the lfp from
    unit_name = char(Int_Identities(i,2));
    wire_expression = '\d{2}';
    unit_wirenum = str2double(char(regexp(unit_name,wire_expression,'match')));

    %NOTE: must ensure that you are always selecting an lfp from the
    %same tetrode from which the spikes were drawn
    allowable_lfp_wires = unit_wirenum + [0,1,2,3];
    w1 = num2str(allowable_lfp_wires(1));
    w2 = num2str(allowable_lfp_wires(2));
    w3 = num2str(allowable_lfp_wires(3));
    w4 = num2str(allowable_lfp_wires(4));
    
    load(filestruct.LFP_file)
        
    vars = who;
    %the below formulation looks hacky, but if you try to consolidate it to
    %['TETFP(0)?(',w1,'|',w2,'|',w3,'|',w4,')\w?(ind|ts)?'], you tend to
    %grab more than just your intended variable
    if unit_wirenum < 9
        lfp_expression = ['TETFP0(',w1,'|',w2,'|',w3,'|',w4,')\w?(ind|ts)?'];
    elseif unit_wirenum == 9
        lfp_expression = ['TETFP(0',w1,'|',w2,'|',w3,'|',w4,')\w?(ind|ts)?'];
    elseif unit_wirenum > 9
        lfp_expression = ['TETFP(',w1,'|',w2,'|',w3,'|',w4,')\w?(ind|ts)?'];
    end
    lfpvars = unique(regexp(vars,lfp_expression,'match','once'));
    
    lfpvars{1}=[]; lfpvars(1)=[]; %removes the empty cell
    
    lfp = eval(lfpvars{1});
    lfp_ind = eval(lfpvars{2});
    lfp_ts = eval(lfpvars{3});
    full_time = ExpandTimestamps(lfp,lfp_ind,lfp_ts,1000);

    %============================================================
    %                       grab spikes
    %============================================================
    %to grab the spikes for this neuron, load the o.g. file and
    %look for the variable name in there
    load(filestruct.file)
    
    names = who;
    if unit_wirenum < 10
        neuron_expression = ['TETSPK0',num2str(unit_wirenum),'\w+'];
    else
        neuron_expression = ['TETSPK',num2str(unit_wirenum),'\w+'];
    end
    names = unique(regexp(names,neuron_expression,'match','once'));
    
    names{1}=[]; names(1)=[]; %removes the empty cell
    ind = find(strcmp(names, unit_name));

    if ~isempty(ind)
        spikes_ts = eval(names{ind});
    else 
    %if spikes_ts is still empty at this point, it's because the
    %spikes variables are listed under a different format
        tet_id = (unit_wirenum + 3)/4; %grab the tetrode (as opposed to wire) id
        unit_name_tet = ['T',num2str(tet_id),'_1',unit_name(end)]; %convert unit name 
        
        names = who;
        neuron_expression = ['T',num2str(tet_id),'_1\w+'];    
        names = unique(regexp(names,neuron_expression,'match','once'));
        
        names{1}=[]; names(1)=[]; %removes the empty cell
        ind = find(strcmp(names, unit_name_tet));

        spikes_ts = eval(names{ind}); %yay! got the spike times

    end

    spikes = ts2binary(spikes_ts,lfp,lfp_ind,lfp_ts,1000);

    %set up the name you will save the table under
    savename = [rat_id, '_', session_id, '_correct.csv'];

    if exist(fullfile(savefileto,savename)) ~=0 
        %if a file with this rat/session combination exists, you want to
        %add the current unit's data to it

        for block = 1:length(event_ts)

            pre_spikes = timeseries_by_event(spikes,event_ts{block},time_surr_event_pre,1);
            dur_spikes = timeseries_by_event(spikes,event_ts{block},time_surr_event_dur,1);
            post_spikes = timeseries_by_event(spikes,event_ts{block},time_surr_event_post,1);
    
            for trial = 1:size(dur_spikes,1)

                %set up spikes
                trial_spikes = [pre_spikes(trial,:)'; dur_spikes(trial,:)'; post_spikes(trial,:)'];
                unit_spikes = [unit_spikes; trial_spikes];

                %every other variable exists if you entered this 
                %if statement, so no need to recreate all the other 
                %column names

            end

        end
        
        %load the previously saved table; to achieve this:
        %list the files in the saved data directory
        filelist = dir(savefileto);

        %find the table file that matches the rat and session ids
        table_file = regexp([filelist.name],savename,'match');
        
        %load that .csv file (using the function "readtable")
        ratsession_df = readtable(table_file{1});
    
        %add the spikes column to the table you just pulled up; 
        ratsession_df.(unit_name) = unit_spikes;

        %TODO:
        %if the new lfp is on a different tetrode from the saved one,
        %rename the lfp var and the filtered versions of it, and save them
        %to the table
        if ~any(lfpvars{1} == string(ratsession_df.Properties.VariableNames))

            %filter the lfp according to preferred frequency bands
            [full_filt_theta, ~, ~] = Butterworth_Hilbert_LR(lfp, 1000, [4,12]);
            [full_filt_beta, ~, ~] = Butterworth_Hilbert_LR(lfp, 1000, [15,25]);
            [full_filt_lowgamma, ~, ~] = Butterworth_Hilbert_LR(lfp, 1000, [35,55]);
            [full_filt_highgamma, ~, ~] = Butterworth_Hilbert_LR(lfp, 1000, [65,90]);

            %you'll cut up the lfp & filtered versions according to the 
            %events
            for block = 1:length(event_ts)
                %now select the correct trial data for odor block of choice specifically
                pre_time = timeseries_by_event(full_time,event_ts{block},time_surr_event_pre,0)';
                dur_time = timeseries_by_event(full_time,event_ts{block},time_surr_event_dur,0)';
                post_time = timeseries_by_event(full_time,event_ts{block},time_surr_event_post,0)';
    
                pre_lfp = timeseries_by_event(lfp,event_ts{block},time_surr_event_pre,1);
                dur_lfp = timeseries_by_event(lfp,event_ts{block},time_surr_event_dur,1);
                post_lfp = timeseries_by_event(lfp,event_ts{block},time_surr_event_post,1);
                
                pre_theta = timeseries_by_event(full_filt_theta,event_ts{block},time_surr_event_pre,1);
                dur_theta = timeseries_by_event(full_filt_theta,event_ts{block},time_surr_event_dur,1);
                post_theta = timeseries_by_event(full_filt_theta,event_ts{block},time_surr_event_post,1);
        
                pre_beta = timeseries_by_event(full_filt_beta,event_ts{block},time_surr_event_pre,1);
                dur_beta = timeseries_by_event(full_filt_beta,event_ts{block},time_surr_event_dur,1);
                post_beta = timeseries_by_event(full_filt_beta,event_ts{block},time_surr_event_post,1);
        
                pre_lowgamma = timeseries_by_event(full_filt_lowgamma,event_ts{block},time_surr_event_pre,1);
                dur_lowgamma= timeseries_by_event(full_filt_lowgamma,event_ts{block},time_surr_event_dur,1);
                post_lowgamma = timeseries_by_event(full_filt_lowgamma,event_ts{block},time_surr_event_post,1);
        
                pre_highgamma = timeseries_by_event(full_filt_highgamma,event_ts{block},time_surr_event_pre,1);
                dur_highgamma = timeseries_by_event(full_filt_highgamma,event_ts{block},time_surr_event_dur,1);
                post_highgamma = timeseries_by_event(full_filt_highgamma,event_ts{block},time_surr_event_post,1);

                for trial = 1:size(dur_lfp,1)

                    %set up time
                    time = [time; pre_time(trial,:)'; dur_time(trial,:)'; post_time(trial,:)'];
    
                    %set up lfp 
                    lfp_odorsamp = [lfp_odorsamp; pre_lfp(trial,:)'; dur_lfp(trial,:)'; post_lfp(trial,:)'];
    
                    %set up filtered lfp
                    filtered_theta = [filtered_theta; pre_theta(trial,:)'; dur_theta(trial,:)'; post_theta(trial,:)'];
                    filtered_beta = [filtered_beta; pre_beta(trial,:)'; dur_beta(trial,:)'; post_beta(trial,:)'];
                    filtered_lowgamma = [filtered_lowgamma; pre_lowgamma(trial,:)'; dur_lowgamma(trial,:)'; post_lowgamma(trial,:)'];
                    filtered_highgamma = [filtered_highgamma; pre_highgamma(trial,:)'; dur_highgamma(trial,:)'; post_highgamma(trial,:)'];

                end

            end

            ratsession_df.(lfpvars{1}) = lfp_odorsamp;
            ratsession_df.([lfpvars{1},'filt_theta']) = filtered_theta;
            ratsession_df.([lfpvars{1},'filt_beta']) = filtered_beta;
            ratsession_df.([lfpvars{1},'filt_lowgamma']) = filtered_lowgamma;
            ratsession_df.([lfpvars{1},'filt_highgamma']) = filtered_highgamma;
        end
        
        %then save the table again, with the same filename as before    
        writetable(ratsession_df,fullfile(savefileto,savename))

    else
 
        %===========================================================
        %                      process lfp
        %===========================================================
        [full_filt_theta, ~, ~] = Butterworth_Hilbert_LR(lfp, 1000, [4,12]);
        [full_filt_beta, ~, ~] = Butterworth_Hilbert_LR(lfp, 1000, [15,25]);
        [full_filt_lowgamma, ~, ~] = Butterworth_Hilbert_LR(lfp, 1000, [35,55]);
        [full_filt_highgamma, ~, ~] = Butterworth_Hilbert_LR(lfp, 1000, [65,90]);
       
        %============================================================
        %             model spike-phase relationships
        %============================================================
        %Note: iterate through odor blocks: interneuron data is analyzed
        %separately
    
        for block = 1:length(event_ts)
            %now select the correct trial data for odor block of choice specifically
            pre_time = timeseries_by_event(full_time,event_ts{block},time_surr_event_pre,0)';
            dur_time = timeseries_by_event(full_time,event_ts{block},time_surr_event_dur,0)';
            post_time = timeseries_by_event(full_time,event_ts{block},time_surr_event_post,0)';

            pre_lfp = timeseries_by_event(lfp,event_ts{block},time_surr_event_pre,1);
            dur_lfp = timeseries_by_event(lfp,event_ts{block},time_surr_event_dur,1);
            post_lfp = timeseries_by_event(lfp,event_ts{block},time_surr_event_post,1);
            
            pre_theta = timeseries_by_event(full_filt_theta,event_ts{block},time_surr_event_pre,1);
            dur_theta = timeseries_by_event(full_filt_theta,event_ts{block},time_surr_event_dur,1);
            post_theta = timeseries_by_event(full_filt_theta,event_ts{block},time_surr_event_post,1);
    
            pre_beta = timeseries_by_event(full_filt_beta,event_ts{block},time_surr_event_pre,1);
            dur_beta = timeseries_by_event(full_filt_beta,event_ts{block},time_surr_event_dur,1);
            post_beta = timeseries_by_event(full_filt_beta,event_ts{block},time_surr_event_post,1);
    
            pre_lowgamma = timeseries_by_event(full_filt_lowgamma,event_ts{block},time_surr_event_pre,1);
            dur_lowgamma= timeseries_by_event(full_filt_lowgamma,event_ts{block},time_surr_event_dur,1);
            post_lowgamma = timeseries_by_event(full_filt_lowgamma,event_ts{block},time_surr_event_post,1);
    
            pre_highgamma = timeseries_by_event(full_filt_highgamma,event_ts{block},time_surr_event_pre,1);
            dur_highgamma = timeseries_by_event(full_filt_highgamma,event_ts{block},time_surr_event_dur,1);
            post_highgamma = timeseries_by_event(full_filt_highgamma,event_ts{block},time_surr_event_post,1);
    
            pre_spikes = timeseries_by_event(spikes,event_ts{block},time_surr_event_pre,1);
            dur_spikes = timeseries_by_event(spikes,event_ts{block},time_surr_event_dur,1);
            post_spikes = timeseries_by_event(spikes,event_ts{block},time_surr_event_post,1);

            %visual check
%             figure;hold on
%             for trial = 1:size(dur_spikes,1)
% 
%                plot([pre_lfp(trial,:), dur_lfp(trial,:), post_lfp(trial,:)] + trial,'k')
%                plot([pre_theta(trial,:), dur_theta(trial,:), post_theta(trial,:)] + trial,'r')
% 
%             end

            current_odors = odorblock{block};
            current_positions = posblock{block};

            %set up the starts of each quarter so you can save the accuracy
            %by quarter data along with the other table variables
            if block == 1

                q1_st = 1; %quarter 1 begins at trial 1 of block 1
                odorstruct = odorblock{1};
        
                q1_odors = odorstruct(1);
                tmp = find(odorstruct ~= q1_odors);
                q1_odors = [q1_odors; odorstruct(tmp(1))];
        
                tmp = find(odorstruct ~= q1_odors(1) & odorstruct ~= q1_odors(2));
                q2_st = tmp(1); %quarter 2 begins on the first trial with 
                                     %an odor deviation from odors in quarter 1
                q2_odors = unique(odorstruct(tmp));

            end

            if block == 2

                q3_st = 1; %quarter 3 begins at trial 1 of block 2
                odorstruct = odorblock{2};

                q3_odors = odorstruct(q3_st);
                tmp = find(odorstruct ~= q3_odors);
                q3_odors = [q3_odors; odorstruct(tmp(1))];
                
                tmp = find(odorstruct ~= q3_odors(1) & odorstruct ~= q3_odors(2));
                q4_st = tmp(1);
                q4_odors = unique(odorstruct(tmp));

            end

            %now fill in table variables trial by trial
            for trial = 1:size(dur_spikes,1)

                %set up time
                time = [time; pre_time(trial,:)'; dur_time(trial,:)'; post_time(trial,:)'];

                %set up lfp 
                lfp_odorsamp = [lfp_odorsamp; pre_lfp(trial,:)'; dur_lfp(trial,:)'; post_lfp(trial,:)'];

                %set up filtered lfp
                filtered_theta = [filtered_theta; pre_theta(trial,:)'; dur_theta(trial,:)'; post_theta(trial,:)'];
                filtered_beta = [filtered_beta; pre_beta(trial,:)'; dur_beta(trial,:)'; post_beta(trial,:)'];
                filtered_lowgamma = [filtered_lowgamma; pre_lowgamma(trial,:)'; dur_lowgamma(trial,:)'; post_lowgamma(trial,:)'];
                filtered_highgamma = [filtered_highgamma; pre_highgamma(trial,:)'; dur_highgamma(trial,:)'; post_highgamma(trial,:)'];

                %set up spikes
                trial_spikes = [pre_spikes(trial,:)'; dur_spikes(trial,:)'; post_spikes(trial,:)'];
                unit_spikes = [unit_spikes; trial_spikes];

                %set up trial_segments 
                pre = repmat({'pre'},length(pre_spikes(trial,:)),1);
                dur = repmat({'dur'},length(dur_spikes(trial,:)),1);
                post = repmat({'post'},length(post_spikes(trial,:)),1);
                segs = [pre;dur;post];
                trial_segment = [trial_segment; segs];

                %set up block labels
                odor_block_labels = [odor_block_labels; repmat(block,length(trial_spikes),1)];

                %set up trial labels 
                trial_labels = [trial_labels; repmat(trial,length(trial_spikes),1)];

                %set up odor id labels
                odor_labels = [odor_labels; repmat(current_odors(trial),length(trial_spikes),1)];

                %set up position id labels
                pos_labels = [pos_labels; repmat(current_positions(trial),length(trial_spikes),1)];

                len_trial = numel(trial_labels);
                gather_len_trials = [gather_len_trials; len_trial];

            end

        end

        %add rat name to table only after the table elements are
        %complete
        rat_names = repmat({rat_id},length(lfp_odorsamp),1);
    
        %save day/session labels
        session_labels = repmat({session_id},length(lfp_odorsamp),1);
    
        %set table up 
        ratsession_df = table(rat_names,time,...
            session_labels,odor_block_labels,trial_labels,trial_segment, ...
            odor_labels, pos_labels);

        %add lfp vars: you should always commit to naming them according to
        %the tetrode wirenum to spot potential errors or mismatches between
        %the saved lfp and the saved units
        ratsession_df.(lfpvars{1}) = lfp_odorsamp;
        ratsession_df.([lfpvars{1},'filt_theta']) = filtered_theta;
        ratsession_df.([lfpvars{1},'filt_beta']) = filtered_beta;
        ratsession_df.([lfpvars{1},'filt_lowgamma']) = filtered_lowgamma;
        ratsession_df.([lfpvars{1},'filt_highgamma']) = filtered_highgamma;

        %add unit var
        ratsession_df.(unit_name) = unit_spikes;

        %set up accuracy by quarter information for the first block
        odorstruct1 = Intervals_by_odorblock.odorblock1.odor_id;
        quarter_edges_list = [q1_st,q2_st,length(odorstruct1)+1];
        for q = 1:length(quarter_edges_list)-1
            
            idx = ratsession_df.trial_labels >= quarter_edges_list(q) ...
                & ratsession_df.trial_labels <= quarter_edges_list(q+1)-1;

            tmptab = ratsession_df(ratsession_df.odor_block_labels == 1 & idx,:);

            accuracy = [accuracy; repmat(acc_list(q),size(tmptab,1),1)]; 
            quarter_labels = [quarter_labels; repmat(q,size(tmptab,1),1)]; 

        end

        odorstruct2 = Intervals_by_odorblock.odorblock2.odor_id;
        quarter_edges_list = [q3_st,q4_st,length(odorstruct2)+1];
        for q = 1:length(quarter_edges_list)-1
            
            idx = ratsession_df.trial_labels >= quarter_edges_list(q) ...
                & ratsession_df.trial_labels <= quarter_edges_list(q+1)-1;

            tmptab = ratsession_df(ratsession_df.odor_block_labels == 2 & idx,:);
            
            %remember acc_list is len = 4, while quarter_edges_list just
            %has the edges for two quarters, this is why you need to add
            %"2" to the acc_list index--this will have you starting at
            %quarter 3 and then 4
            accuracy = [accuracy; repmat(acc_list(q+2),size(tmptab,1),1)]; 
            quarter_labels = [quarter_labels; repmat(q+2,size(tmptab,1),1)]; 

        end

        ratsession_df.accuracy = accuracy;
        ratsession_df.quarter_labels = quarter_labels;

        writetable(ratsession_df,fullfile(savefileto,savename))

        clearvars -except Rat Int_Identities i ...
            time_surr_event_pre time_surr_event_dur time_surr_event_post ...
            savefileto

    end

end



