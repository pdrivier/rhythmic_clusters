%Data Preprocessing: for spike-phase modeling in python
%Goal: grab interneuron spikes and rhythmic phases during approach,
%trial-matched to odor sampling epochs


%this script will organize the data into .csv files 
%columns:
%rat_id
%session_id
%odorblock/condition
%trial_labels
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

%make sure matlab knows where you'll be pulling the raw data from
addpath data\
addpath lfp_data\

%set up your saving directory
savefileto = 'python_spkphase_approach_ints\';

%load the list of interneurons you'll organize data from
load('unique_ints_by_daytetfrate.mat') 

%decide ahead of time the time interval edges surrounding the event ts
time_surr_event_dur = [0 1.5]; %this gives time for doors down relative to 
                               %the context signaling event ts

Int_Identities = unique_cells; 

%you'll generate two dataframes per rat-session combination: 
%one containing LFP-based data
%the other containing position/velocity-based data
%they can't be combined into the same dataframe because the rows of each
%dataframe correspond to a given time sample, and position and lfp
%variables are sampled at different rates, so you end up with a different
%number of rows per set of lfp/position data

for i = 1:length(Int_Identities)
    sprintf('on %d out of %d',i,length(Int_Identities))

    %initialize the variables for lfp-based dataframe
    rat_names = [];
    session_labels = [];
    odor_block_labels = [];
    trial_labels = [];
    trial_segment = [];
    accuracy = [];
    quarter_labels = [];

    odor_labels = []; %keep these so you know how many trials per 
                      %condition were available; odor released on trial
    pos_labels = [];  %odor port being sampled from

    time = [];
    lfp_approach = [];

    filtered_theta = [];
    filtered_beta = [];
    filtered_lowgamma = [];
    filtered_highgamma = [];

    unit_spikes = [];

    gather_len_trials = []; %this serves as a nice check

    %initialize variables for position-based dataframe
    vsm = [];%incorporate the actual velocity traces for the epoch
    vraw = [];
    vposts = [];
    posx = [];%and also include where the rat is located during this time
    posy = [];
    odor_block_labels_poslength = [];
    trial_labels_poslength = [];
    
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


    epoch_duration = 1.5; %duration of an odor sampling trial, which you 
                          %will want to match in the approach!
    %all block 1
    %NOTE: Intervals_by_odorblock context signaling variable already
    %accounts for the time it takes for the doors to fully go down (6s post
    %the event marker). On the other hand, the Intervals context signaling
    %variable contains the timestamp of the actual event marker, so it
    %always lags 6s behind the first column of timestamps in
    %Intervals_by_odorblock
    corrts_block01 = Intervals_by_odorblock.odorblock1.correct(:,1);
    corrinds_block01 = corrts_block01 ~= 0;
    corrts_block01 = corrts_block01(corrinds_block01);
    tmp_event_ts_block01 = Intervals_by_odorblock.odorblock1.context_signaling(corrinds_block01,1);
    tmp_event_ts_block01(:,2) = corrts_block01(:,1) - 0.750;
    
    %be a bit more generous in terms of allowed running time for trials
    %where the run was very short
    short_run_trials = find(tmp_event_ts_block01(:,2) - tmp_event_ts_block01(:,1) <= 0);
    tmp_event_ts_block01(short_run_trials,2) = corrts_block01(short_run_trials,1) - 0.250;
    
    %for trials that are still too short despite having relaxed proximity
    %to the odor port
    trial_tmpdur = tmp_event_ts_block01(:,2) - tmp_event_ts_block01(:,1);
    still_too_short =  find(trial_tmpdur < (epoch_duration + 0.100));
    deltas = epoch_duration - trial_tmpdur(still_too_short);
    tmp_event_ts_block01(still_too_short,1) = tmp_event_ts_block01(still_too_short,1) - (deltas+.750);
    tmp_event_ts_block01(still_too_short,2) = corrts_block01(still_too_short,1) - 0.750;

    %then check for trials where for some reason odor sampling start time
    %is earlier than the earliest context signaling (post door-down -- did
    %the rat somehow jump out the door before the door was fully down and
    %completed the trial super quickly? regardless, when you find a trial
    %like this, subtract as much from the context signaling start time as 
    %its delta from the odor sampling time, plus 2 seconds, up to a cap of
    %6 seconds(which is when the doors start coming down). TODO: implemnt
    %the cap!
    late_ctxt_trials = find(tmp_event_ts_block01(:,1) >= corrts_block01);
    delta = tmp_event_ts_block01(late_ctxt_trials,1) - corrts_block01(late_ctxt_trials);
    delta = delta + 2;
    tmp_event_ts_block01(late_ctxt_trials,1) =  tmp_event_ts_block01(late_ctxt_trials) - delta;
    tmp_event_ts_block01(late_ctxt_trials,2) =  corrts_block01(late_ctxt_trials,1) - 0.750;

    odors_block01 = Intervals_by_odorblock.odorblock1.odor_id(corrinds_block01);

    pos_block01 = Intervals_by_odorblock.odorblock1.pos_id(corrinds_block01);

    load(filestruct.file) %to get position vars to identify high velocity 
                          %epochs

    %get the highest velocity epochs when the rat is also facing the
    %direction of the port (to avoid high velocity jumps where the rat is
    %doing random frantic shit)
    event_ts_block01 = zeros(length(tmp_event_ts_block01),2);
    dt = unique(1/LED1_X_ts); dt = dt(dt~=0);
    for trial=1:length(tmp_event_ts_block01)

        pos_ts = LED1_X_ts(LED1_X_ts>=tmp_event_ts_block01(trial,1)&LED1_X_ts<=tmp_event_ts_block01(trial,2));
        posX = LED1_X(LED1_X_ts>=tmp_event_ts_block01(trial,1)&LED1_X_ts<=tmp_event_ts_block01(trial,2));
        posY = LED1_Y(LED1_X_ts>=tmp_event_ts_block01(trial,1)&LED1_X_ts<=tmp_event_ts_block01(trial,2));
        
        destination_port = pos_block01(trial);

        [newst,newed,vrawwin,vsmwin,postswin,posxwin,posywin] = find_high_velocity_epochs(pos_ts,posX,posY,epoch_duration,dt,destination_port);

        event_ts_block01(trial,1) = newst;
        event_ts_block01(trial,2) = newed;

        %TODO: the challenge with these new variables is their variable
        %length on each trial...need to save them differently
        vraw01{trial} = vrawwin;
        vsm01{trial} = vsmwin;
        postsveloc01{trial} = postswin;
        posx01{trial} = posxwin;
        posy01{trial} = posywin;
    end

    %all block 2
    corrts_block02 = Intervals_by_odorblock.odorblock2.correct(:,1);
    corrinds_block02 = corrts_block02 ~= 0;
    corrts_block02 = corrts_block02(corrinds_block02);
    tmp_event_ts_block02 = Intervals_by_odorblock.odorblock2.context_signaling(corrinds_block02,1);
    tmp_event_ts_block02(:,2) = corrts_block02 - 0.750;

    %be a bit more generous in terms of allowed running time for trials
    %where the run was very short
    short_run_trials = find(tmp_event_ts_block02(:,2) - tmp_event_ts_block02(:,1) <= 0);
    tmp_event_ts_block02(short_run_trials,2) = corrts_block02(short_run_trials,1) - 0.250;
    
    %for trials that are still too short despite having relaxed proximity
    %to the odor port
    trial_tmpdur = tmp_event_ts_block02(:,2) - tmp_event_ts_block02(:,1);
    still_too_short =  find(trial_tmpdur < (epoch_duration + 0.100));
    deltas = epoch_duration - trial_tmpdur(still_too_short);
    tmp_event_ts_block02(still_too_short,1) = tmp_event_ts_block02(still_too_short,1) - (deltas+.750);
    tmp_event_ts_block02(still_too_short,2) = corrts_block02(still_too_short,1) - 0.750;

    %then check for trials where for some reason odor sampling start time
    %is earlier than the earliest context signaling (post door-down -- did
    %the rat somehow jump out the door before the door was fully down and
    %completed the trial super quickly? regardless, when you find a trial
    %like this, subtract as much from the context signaling start time as 
    %its delta from the odor sampling time, plus 2 seconds, up to a cap of
    %6 seconds(which is when the doors start coming down). TODO: implemnt
    %the cap!
    late_ctxt_trials = find(tmp_event_ts_block02(:,1) >= corrts_block02);
    delta = tmp_event_ts_block02(late_ctxt_trials,1) - corrts_block02(late_ctxt_trials);
    delta = delta + 2;
    tmp_event_ts_block02(late_ctxt_trials,1) =  tmp_event_ts_block02(late_ctxt_trials) - delta;
    tmp_event_ts_block02(late_ctxt_trials,2) =  corrts_block02(late_ctxt_trials,1) - 0.750;

    odors_block02 = Intervals_by_odorblock.odorblock2.odor_id(corrinds_block02);

    pos_block02 = Intervals_by_odorblock.odorblock2.pos_id(corrinds_block02);
 
    %find the high velocity epochs
    event_ts_block02 = zeros(length(tmp_event_ts_block02),2);

    for trial=1:length(tmp_event_ts_block02)

        pos_ts = LED1_X_ts(LED1_X_ts>=tmp_event_ts_block02(trial,1)&LED1_X_ts<=tmp_event_ts_block02(trial,2));
        posX = LED1_X(LED1_X_ts>=tmp_event_ts_block02(trial,1)&LED1_X_ts<=tmp_event_ts_block02(trial,2));
        posY = LED1_Y(LED1_X_ts>=tmp_event_ts_block02(trial,1)&LED1_X_ts<=tmp_event_ts_block02(trial,2));
       
        destination_port = pos_block02(trial);

        [newst,newed,vrawwin,vsmwin,postswin,posxwin,posywin] = find_high_velocity_epochs(pos_ts,posX,posY,epoch_duration,dt,destination_port);

        event_ts_block02(trial,1) = newst;
        event_ts_block02(trial,2) = newed;
        vraw02{trial} = vrawwin;
        vsm02{trial} = vsmwin;
        postsveloc02{trial} = postswin;
        posx02{trial} = posxwin;
        posy02{trial} = posywin;

    end
    
    %this will help you iterate by block later
    event_ts = {event_ts_block01(:,1),event_ts_block02(:,1)};
    odorblock = {odors_block01,odors_block02};
    posblock = {pos_block01,pos_block02};
    vrawblock = {vraw01,vraw02};
    vsmblock = {vsm01,vsm02};
    postsvelocblock = {postsveloc01,postsveloc02};
    posxblock = {posx01,posx02};
    posyblock = {posy01,posy02};


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
%     load(filestruct.file)
    
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

    %set up the name you will save the lfp table under
    savename = [rat_id, '_', session_id, '_approach.csv'];
   

    if exist(fullfile(savefileto,savename)) ~=0 
        %if a file with this rat/session combination exists, you want to
        %add the current unit's data to it

        for block = 1:length(event_ts)

            dur_spikes = timeseries_by_event(spikes,event_ts{block},time_surr_event_dur,1);
    
            for trial = 1:size(dur_spikes,1)

                %set up spikes
                trial_spikes = dur_spikes(trial,:)';
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

            %filter the lfp
            [full_filt_theta, ~, ~] = Butterworth_Hilbert_LR(lfp, 1000, [4,12]);
            [full_filt_beta, ~, ~] = Butterworth_Hilbert_LR(lfp, 1000, [15,25]);
            [full_filt_lowgamma, ~, ~] = Butterworth_Hilbert_LR(lfp, 1000, [35,55]);
            [full_filt_highgamma, ~, ~] = Butterworth_Hilbert_LR(lfp, 1000, [65,90]);

            %you'll cut up the lfp & filtered versions according to the 
            %events
            for block = 1:length(event_ts)
                %now select the correct trial data for odor block of choice specifically
                dur_time = timeseries_by_event(full_time,event_ts{block},time_surr_event_dur,0)';
    
                dur_lfp = timeseries_by_event(lfp,event_ts{block},time_surr_event_dur,1);
                
                dur_theta = timeseries_by_event(full_filt_theta,event_ts{block},time_surr_event_dur,1);
        
                dur_beta = timeseries_by_event(full_filt_beta,event_ts{block},time_surr_event_dur,1);
        
                dur_lowgamma= timeseries_by_event(full_filt_lowgamma,event_ts{block},time_surr_event_dur,1);
        
                dur_highgamma = timeseries_by_event(full_filt_highgamma,event_ts{block},time_surr_event_dur,1);

                for trial = 1:size(dur_lfp,1)

                    %set up time
                    time = [time; dur_time(trial,:)'];
    
                    %set up lfp 
                    lfp_approach = [lfp_approach; dur_lfp(trial,:)'];
    
                    %set up filtered lfp
                    filtered_theta = [filtered_theta; dur_theta(trial,:)'];
                    filtered_beta = [filtered_beta; dur_beta(trial,:)'];
                    filtered_lowgamma = [filtered_lowgamma; dur_lowgamma(trial,:)'];
                    filtered_highgamma = [filtered_highgamma; dur_highgamma(trial,:)'];

                end

            end

            ratsession_df.(lfpvars{1}) = lfp_approach;
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
            dur_time = timeseries_by_event(full_time,event_ts{block},time_surr_event_dur,0)';

            dur_lfp = timeseries_by_event(lfp,event_ts{block},time_surr_event_dur,1);
            
            dur_theta = timeseries_by_event(full_filt_theta,event_ts{block},time_surr_event_dur,1);
    
            dur_beta = timeseries_by_event(full_filt_beta,event_ts{block},time_surr_event_dur,1);
    
            dur_lowgamma= timeseries_by_event(full_filt_lowgamma,event_ts{block},time_surr_event_dur,1);
    
            dur_highgamma = timeseries_by_event(full_filt_highgamma,event_ts{block},time_surr_event_dur,1);
    
            dur_spikes = timeseries_by_event(spikes,event_ts{block},time_surr_event_dur,1);

            %visual check
            figure;hold on
            for trial = 1:size(dur_spikes,1)

               plot(dur_lfp(trial,:) + trial,'k')
               plot(dur_theta(trial,:) + trial,'r')

            end

            current_odors = odorblock{block};
            current_positions = posblock{block};

            %grab this odor block's position & velocity vars
            current_vsm = vsmblock{block};
            current_vraw = vrawblock{block};
            current_posts = postsvelocblock{block};
            current_posx = posxblock{block};
            current_posy = posyblock{block};
            
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
                time = [time; dur_time(trial,:)'];

                %set up lfp 
                lfp_approach = [lfp_approach; dur_lfp(trial,:)'];

                %set up filtered lfp
                filtered_theta = [filtered_theta; dur_theta(trial,:)'];
                filtered_beta = [filtered_beta; dur_beta(trial,:)'];
                filtered_lowgamma = [filtered_lowgamma; dur_lowgamma(trial,:)'];
                filtered_highgamma = [filtered_highgamma; dur_highgamma(trial,:)'];

                %set up spikes
                trial_spikes = dur_spikes(trial,:)';
                unit_spikes = [unit_spikes; trial_spikes];

                %set up trial_segments 
                dur = repmat({'dur'},length(dur_spikes(trial,:)),1);
                
                trial_segment = [trial_segment; dur];


                %set up block labels
                odor_block_labels = [odor_block_labels; repmat(block,length(trial_spikes),1)];
                odor_block_labels_poslength = [odor_block_labels_poslength;...
                    repmat(block,length(current_vsm{trial}))];

                %set up trial labels 
                trial_labels = [trial_labels; repmat(trial,length(trial_spikes),1)];
                trial_labels_poslength = [trial_labels_poslength; ...
                    repmat(trial,length(current_vsm{trial}))];

                %set up odor id labels
                odor_labels = [odor_labels; repmat(current_odors(trial),length(trial_spikes),1)];

                %set up position id labels
                pos_labels = [pos_labels; repmat(current_positions(trial),length(trial_spikes),1)];

                %set up position/veloc vars
                vsm = [vsm; current_vsm{trial}];
                vraw = [vraw; current_vraw{trial}];
                vposts = [vposts; current_posts{trial}];
                posx = [posx; current_posx{trial}];
                posy = [posy; current_posy{trial}];

                len_trial = numel(trial_labels);
                gather_len_trials = [gather_len_trials; len_trial];

            end

        end

        %add rat name to table only after the table elements are
        %complete
        rat_names = repmat({rat_id},length(lfp_approach),1);
    
        %save day/session labels
        session_labels = repmat({session_id},length(lfp_approach),1);
    
        %set lfp table up 
        ratsession_df = table(rat_names,time,...
            session_labels,odor_block_labels,trial_labels,trial_segment, ...
            odor_labels, pos_labels, max_v, mean_v, med_v, min_v, std_v);

        %add lfp vars: you should always commit to naming them according to
        %the tetrode wirenum to spot potential errors or mismatches between
        %the saved lfp and the saved units
        ratsession_df.(lfpvars{1}) = lfp_approach;
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


        %set position/velocity table up 
        rat_names_poslength = repmat({rat_id},length(vsm),1);
        session_labels_poslength = repmat({session_id},length(vsm),1);


        ratpos_df = table(rat_names_poslength,vposts,...
            session_labels_poslength,odor_block_labels_poslength,...
            trial_labels_poslength,vraw, vsm, posx, posy);

        possavename = [rat_id, '_', session_id, '_posvars_approach.csv'];
        writetable(ratpos_df,fullfile(savefileto,pos_savename))


        clearvars -except Rat Int_Identities i ...
            time_surr_event_pre time_surr_event_dur time_surr_event_post ...
            savefileto

    end

end



