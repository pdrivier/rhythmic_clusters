function [epoch_start_ts,epoch_end_ts,maxv,meanv,medv,minv,stdv] = find_high_velocity_epochs(pos_ts,posX,posY,epoch_duration,dt,destination_port)
%pos_ts: vector of timestamps to be searched for high velocities
%         which you will search for the highest velocity epoch
%epoch_duration: scalar, dictates how long the epoch you find should be
%posX: vector, of tracking data, x-coordinate
%posY: vector, of tracking data, y-coordinate
%dt: scalar, of sampling time between successive tracking variable entries
%destination port: scalar, this trial's corresponding pos id of the sampled 
%                  odor port 


%fix any issues in the tracking
[posX,posY] = NearestNeighborforMaze_LR(posX,posY);

%compute the velocity for the full trial
delta_x = diff(posX)*2.5; %2.5cm for every unit change in x
delta_y = diff(posY)*2.5; %2.5cm for every unit change in y
dist = sqrt( (delta_x.^2) + (delta_y.^2) );
vraw = dist./dt;

%smooth the velocity vector
vsm = smoothdata(vraw,"gaussian",40);

%make the time vector match the velocity vector
pos_ts = pos_ts(2:end);

%get head direction angles
dx = diff(posX);
dy = diff(posY);
hd_angles = atan2(dy,dx);
hd_angles = hd_angles + pi;

if destination_port == 3 || destination_port == 4
    hd_ind = find(hd_angles>=4);
    hd_ind = [hd_ind; find(hd_angles<=2)];
    hd_ind = sort(hd_ind);

elseif destination_port == 1 || destination_port == 2
    hd_ind = find(hd_angles>= 2 & hd_angles<= 4);

end
   
%get the distribution of velocities and find the top ten percent
% figure;histogram(vsm)

%get the timestamps corresponding to the velocities in the top ten
%percentile of the distribution of smoothed velocities
sortvsm = sort(vsm(hd_ind), 'descend');
topvs = sortvsm(1:ceil(length(sortvsm)*0.1));
topinds = [];
for t=1:numel(topvs)
    topinds = [topinds; find(vsm==topvs(t));];
end


%visualize the spatial variables
% [AllPosX,AllPosY] = NearestNeighborforMaze_LR(LED1_X,LED1_Y);
% figure;plot(AllPosX,AllPosY,color=[.9 .9 .9])
% hold on;plot(posX,posY,LineWidth=5)
% 
% %visualize both velocity vectors and the corresponding head directions
% %during the run
% figure;plot(pos_ts,vraw);axis tight
% hold on;plot(pos_ts,vsm)
% hold on;xline(pos_ts(hd_ind),'k')
% hold on;xline(pos_ts(topinds),'r') 


%find out how close each of these top indices is to the end or start of the 
%full segment
dur_to_edges = []; %index, dur to start, dur to end
for t=1:length(topinds)
    dur_to_edges = [dur_to_edges; [topinds(t), ...
                                   pos_ts(topinds(t)) - pos_ts(1),...
                                   pos_ts(end) - pos_ts(topinds(t))]];
end

%grab multiple windows of size `epoch_duration`
w=1;
for t=1:size(dur_to_edges,1)
    
    if dur_to_edges(t,2) >= epoch_duration && dur_to_edges(t,3) >= epoch_duration
        
        %set up your index in the middle of the window
        potential_window_edges(w,:) = [pos_ts(dur_to_edges(t,1))-epoch_duration/2, ...
                                     pos_ts(dur_to_edges(t,1))+epoch_duration/2];

       
    elseif dur_to_edges(t,2) >= epoch_duration && dur_to_edges(t,3) < epoch_duration

        %set up your index to be the end of the window
        potential_window_edges(w,:) = [pos_ts(dur_to_edges(t,1))-epoch_duration, pos_ts(dur_to_edges(t,1))];

       
    elseif dur_to_edges(t,2) < epoch_duration && dur_to_edges(t,3) >= epoch_duration

        %set up your index to be the start of the window
        potential_window_edges(w,:) = [pos_ts(dur_to_edges(t,1)),pos_ts(dur_to_edges(t,1))+epoch_duration];

     

    elseif dur_to_edges(t,2) < epoch_duration && dur_to_edges(t,3) < epoch_duration && pos_ts(end) - pos_ts(1) >= epoch_duration

        %set up your index somewhere inside the window
        end_at = pos_ts(end-1);
        potential_window_edges(w,:) = [end_at-epoch_duration, end_at];

    elseif dur_to_edges(t,2) < epoch_duration && dur_to_edges(t,3) < epoch_duration && pos_ts(end) - pos_ts(1) < epoch_duration

        potential_window_edges(w,:) = [0,0];

     w = w+1;  

        
    end
    
end
if sum(potential_window_edges)==0

    %just call it a day and grab whatever epoch precedes the odor sampling
    %epoch and is the correct duration
    epoch_start_ts = pos_ts(1);
    epoch_end_ts = pos_ts(1) + epoch_duration;
    
    maxv = max(vsm(pos_ts>=epoch_start_ts & pos_ts<=epoch_end_ts));
    meanv = mean(vsm(pos_ts>=epoch_start_ts & pos_ts<=epoch_end_ts));
    medv = median(vsm(pos_ts>=epoch_start_ts & pos_ts<=epoch_end_ts));
    minv = min(vsm(pos_ts>=epoch_start_ts & pos_ts<=epoch_end_ts));
    stdv = std(vsm(pos_ts>=epoch_start_ts & pos_ts<=epoch_end_ts));
 

else
    %select the window with the highest mean smoothed velocity values 
    %throughout
    for w=1:size(potential_window_edges,1)
        st = potential_window_edges(w,1);
        ed = potential_window_edges(w,2);
        velocity_windows(w,:) = mean(vsm(pos_ts>=st & pos_ts<=ed));
    
    end
    
    
    epoch_start_ts = potential_window_edges(velocity_windows == max(velocity_windows),1);
    epoch_end_ts = potential_window_edges(velocity_windows == max(velocity_windows),2);

    maxv = max(vsm(pos_ts>=epoch_start_ts & pos_ts<=epoch_end_ts));
    meanv = mean(vsm(pos_ts>=epoch_start_ts & pos_ts<=epoch_end_ts));
    medv = median(vsm(pos_ts>=epoch_start_ts & pos_ts<=epoch_end_ts));
    minv = min(vsm(pos_ts>=epoch_start_ts & pos_ts<=epoch_end_ts));
    stdv = std(vsm(pos_ts>=epoch_start_ts & pos_ts<=epoch_end_ts));
    
    % hold on;xline([epoch_start_ts,epoch_end_ts],'b')
    % hold on;yline(velocity_windows(velocity_windows==max(velocity_windows)),'blue')
    % figure(10);plot(posX(pos_ts>=epoch_start_ts & pos_ts<=epoch_end_ts),posY(pos_ts>=epoch_start_ts & pos_ts<=epoch_end_ts),'k')

   
    %output the raw & smoothed velocities from the selected window, the
    %corresponding pos_ts for the window, and the corresponding posx and
    %posy -- values (hopefully) in cm/s due to distance computation at the
    %top of this function (using the 2.5 cm x 2.5 cm figure, check this
    %with Lara)
    vraw_window = vraw(pos_ts>=epoch_start_ts & pos_ts<=epoch_end_ts);
    vsm_window = vsm(pos_ts>=epoch_start_ts & pos_ts<=epoch_end_ts);
    posts_window = pos_ts(pos_ts>=epoch_start_ts & pos_ts<=epoch_end_ts);
    posx_window = posX(pos_ts>=epoch_start_ts & pos_ts<=epoch_end_ts);
    posy_window = posY(pos_ts>=epoch_start_ts & pos_ts<=epoch_end_ts);

end


end