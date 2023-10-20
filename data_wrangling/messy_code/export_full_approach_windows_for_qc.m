%export full approach windows & lfps during these epochs to then quality
%control them and define the approach epochs in python

clear 
close all
dbstop if error

%make sure matlab knows where you'll be pulling the raw data from
addpath data\
addpath lfp_data\
addpath figures\quality_control_approach\

%set up your saving directory
savefileto = 'python_spkphase_approach_ints\';
if ~exist(savefileto,'dir') 
    mkdir(savefileto)
end
addpath(savefileto)

%load the list of interneurons you'll organize data from
load('unique_ints_by_daytetfrate.mat') 

Int_Identities = unique_cells; 

for i = 1:length(Int_Identities)
    sprintf('on %d out of %d',i,length(Int_Identities))
    
    %initialize position/veloc variable dataframe columns
    posx = [];
    posy = [];
    vposts = [];
    vraw = [];
    vsm = [];
    hdir = [];
    odor_port = [];

    rat_names_poslen = [];
    session_labels_poslen = [];
    odor_block_labels_poslen = [];
    trial_labels_poslen = [];



    %initialize lfp dataframe columns
    time = [];
    lfp_approach = [];

    rat_names_lfplen = [];
    session_labels_lfplen = [];
    odor_block_labels_lfplen = [];
    trial_labels_lfplen = [];
           


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

    pos_block01 = Intervals_by_odorblock.odorblock1.pos_id(corrinds_block01);

    load(filestruct.file) %to get position vars to compute velocity from 
  
 
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


    pos_block02 = Intervals_by_odorblock.odorblock2.pos_id(corrinds_block02);

 
    %put the events into the same cell array to iterate over later
    event_ts{1} = tmp_event_ts_block01;
    event_ts{2} = tmp_event_ts_block02;
    pos_block{1} = pos_block01;
    pos_block{2} = pos_block02;

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

    dt = unique(1/LED1_X_ts); dt = dt(dt~=0);
    for block = 1:length(event_ts)

        current_block_trials = event_ts{block};
        current_block_pos = pos_block{block};

        for trial=1:length(current_block_trials)
    
            startts = current_block_trials(trial,1);
            endts = current_block_trials(trial,2);
    
            pos_ts = LED1_X_ts(LED1_X_ts >= startts & LED1_X_ts <= endts);
            posX = LED1_X(LED1_X_ts >= startts & LED1_X_ts <= endts);
            posY = LED1_Y(LED1_X_ts >= startts & LED1_X_ts <= endts);

            %fix any issues in the tracking
            [posX,posY] = NearestNeighborforMaze_LR(posX,posY);
            
            %to assess heading direction, get the destination
            destination_port = current_block_pos(trial);
    
            %compute the velocity for the full trial
            delta_x = diff(posX)*2.5; %2.5cm for every unit change in x
            delta_y = diff(posY)*2.5; %2.5cm for every unit change in y
            dist = sqrt( (delta_x.^2) + (delta_y.^2) );
            vrawwin = dist./dt;
            
            %smooth the velocity vector
            vsmwin = smoothdata(vraw,"gaussian",40);
            
            %make the time vector match the velocity vector
            pos_ts = pos_ts(2:end);
            
            %get head direction angles
            dx = diff(posX)*2.5;
            dy = diff(posY)*2.5;
            hd_angles = atan2(dy,dx);
            hd_angles = hd_angles + pi;
            
            if destination_port == 3 || destination_port == 4
                hd_ind = find(hd_angles>=4);
                hd_ind = [hd_ind; find(hd_angles<=2)];
                hd_ind = sort(hd_ind);
            
            elseif destination_port == 1 || destination_port == 2
                hd_ind = find(hd_angles>= 2 & hd_angles<= 4);
            
            end

            %initialize position/veloc variable dataframe columns
            posx = [posx; posX(2:end)]; %needs same len as velocity
            posy = [posy; posY(2:end)];
            vposts = [vposts; pos_ts];
            vraw = [vraw; vrawwin];
            vsm = [vsm; vsmwin];
            hdir = [hdir; hd_angles];
            odor_port = [odor_port; repmat(destination_port,length(vposts),1)];

            rat_names_poslen = [rat_names_poslen; repmat({rat_id},length(vposts),1)];
            session_labels_poslen = [session_labels_poslen; repmat({session_id},length(vposts),1)];
            odor_block_labels_poslen = [odor_block_labels_poslen; repmat(block,length(vposts),1)];
            trial_labels_poslen = [trial_labels_poslen; repmat(trial_labels,length(vposts),1)];
        
            %grab the corresponding lfp and lfp timestamps
            lfp_snippet = lfp(full_time >= startts & full_time <= endts);
            lfp_time = full_time(full_time >= startts & full_time <= endts);
        
            time = [time; lfp_time];
            lfp_approach = [lfp_approach; lfp_snippet];

            rat_names_lfplen = [rat_names_lfplen; repmat({rat_id},length(time),1)];
            session_labels_lfplen = [session_labels_lfplen; repmat({session_id},length(time),1)];
            odor_block_labels_lfplen = [odor_block_labels_lfplen; repmat(block,length(time),1)];
            trial_labels_lfplen = [trial_labels_lfplen; repmat(trial_labels,length(time),1)];
           
            
            %visualize both velocity vectors and the corresponding head directions
            %during the run
            figure;
            subplot(211)
            plot(pos_ts,vrawwin,color=[.9 .9 .9]);axis tight
            hold on;plot(pos_ts,vsmwin,color='k')
            hold on;xline(pos_ts(hd_ind),'r')
            
            tmp_lfp = lfp(full_time>=startts & full_time<=endts);
            tmp_lfp_ts = full_time(full_time>=startts & full_time<=endts);
    
            subplot(212)
            plot(lfp_time,lfp_snippet,color='k');axis tight
            
            title([rat_id, session_id, ' ', lfpvars{1}, ' block1 trial ',num2str(trial)])
    
            figname = [rat_id, session_id, '_', lfpvars{1}, '_block1_trial',num2str(trial),'.pdf'];
            saveas(figure(1),fullfile('figures/quality_control_approach/',figname))
    
            close
        end
    end


end
