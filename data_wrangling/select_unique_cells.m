function [unique_cells] = select_unique_cells(list_to_select_from,delimratsess,delimtetunit,min_accept_distance)
%SELECT_UNIQUE_CELLS, under the assumption that an interneuron (or other
%cell type) can be recorded multiple consecutive days, even after turning
%the tetrode, this function will help you narrow down the list to what are
%likely to be unique interneurons (or other cell type). 

%Given the sparsity of interneurons in the hippocampus, it is unlikely that
%you were fortunate enough to have recorded multiple *unique* interneurons 
%in the ***same rat***, on the ***same tetrode***, across consecutive days.
%So these are the parameters you will be on the lookout for. 

%INPUT
%list_to_select_from,  (n x 3) cell array, where n corresponds to the
%                       number of rat_session combinations
%delim,                 string, to indicate where to split the rat_session
%                       name at (e.g. delimratsess = 'D', delimtetunit =
%                       'TETSPK')
%min_accept_distance,   scalar, minimum number of days between recordings 
%                       of the same tetrode that we are willing to accept

%TODO - START HERE next work session
            %Another issue with the code is that
            %the last six cells in the dataset don't have a firing rate
            %attached to them--perhaps it's worth going back to take care
            %of this first in another snippet of code that searches out the
            %corresponding file in Rat struct? remember, these wouldn't be
            %saved in the 'frate' subfield of the Rat struct because I cut
            %them later!

            %Another issue is that sometimes, especially for LH16, the same
            %cell could have been spike sorted as a different unit letter,
            %just because he had so many pyramidal cells that the cutting
            %order for the interneuron would have changed, so need to do a
            %second pass to catch and pluck those out!

%first get the list of all unique rats
for n = 1:length(list_to_select_from)
    
    split_rat_sess = split(list_to_select_from(n,1),delimratsess);
    split_tet_unit = split(list_to_select_from(n,2),delimtetunit);
    
    wireset = split_tet_unit{2};
    

    
    rat = split_rat_sess{1};
    session = str2num(split_rat_sess{2});
    wirenum = str2num(wireset(1:2));
    wireletter = wireset(end);
    frate = list_to_select_from{n,3};


    reformatted_list{n,1} = rat;
    reformatted_list{n,2} = session;
    reformatted_list{n,3} = wirenum;
    reformatted_list{n,4} = wireletter;
    reformatted_list{n,5} = frate;


end
list_of_rats = unique([reformatted_list(:,1)]);

tmp_cell_array = [];
%with this list, you will now filter the newly reformatted cell array
for r =1:length(list_of_rats)
    
    rat = list_of_rats{r};
    keep_ind = regexp(reformatted_list(:,1),rat,'match');
    keep_ind = find(~cellfun(@isempty,keep_ind));

    %filter the cell array for just the rows for this rat
    subrat = reformatted_list(keep_ind,:);

    %now you'll have to determine if the tetrode number occurs 
    %accross multiple sessions
    unique_sessions = unique([subrat{:,2}]);

    %here is where you will need to determine if not just the wirenum, but 
    %the wirenum + wireletter combination is present across days! 
    wirenumletter = [];
    for i = 1:size(subrat,1)
        if subrat{i,3} < 10
            num = ['0',num2str(subrat{i,3})];
        else 
            num = num2str(subrat{i,3});
        end
        letter = subrat{i,4};
        wirenumletter = [wirenumletter; [num,letter]];
        subrat{i,6} = [num,letter];
    end
    unique_wires_letters = unique(wirenumletter,'rows');

    for wl = 1:length(unique_wires_letters)
        %count the number of times this tetrode appears in the rat's
        %dataset
        wirenumletter_expression = unique_wires_letters(wl,:);
        
        %find the number of perfect matches in subrat
        wirenum_matches = regexp(subrat(:,6),wirenumletter_expression,'match');
        n_appearances = numel(find(~cellfun(@isempty,wirenum_matches)));

        if n_appearances == 1

            subwire = subrat(find(~cellfun(@isempty,wirenum_matches)),:);
            tmp_cell_array = [tmp_cell_array; subwire];

        elseif n_appearances == 2 
            %if there are only two options, check for the number of days between 
            %recordings, and grab both if they exceed the minimum acceptable
            %distance, or grab only the highest frate one if fewer days than
            %the acceptable minimum distance
            subwire = subrat(find(~cellfun(@isempty,wirenum_matches)),:);
    
            %find the difference between the session days
            days_distance = diff([subwire{:,2}]);
        
            if days_distance > min_accept_distance

                tmp_cell_array = [tmp_cell_array;subwire];

            else

                frate_list = [subwire{:,5}];
                highest_fr_recording =  subwire(frate_list == max(frate_list),:);
                tmp_cell_array = [tmp_cell_array; highest_fr_recording];

            end

        elseif n_appearances > 2
            %if the tetrode has been included more than once, you'll have to
            %figure out if this happened over the course of consecutive days
            %as a rough ball park, more than 10 days distance will be considered
            %un-consecutive; the min_accept_distance can be varied depending on
            %priors over a specific region, or cell type, or known size of turns
            %(e.g. don't expect to be able to record from the same granule cell
            %after a 35 um turn across days, but yes an interneuron)
    
            %side note: the best way to do this would actually be to do it by 
            %tetrode turning distance (although even this is not always accurate
            %depending on the quality of the tops and the precision of the
            %turn, as well as fractional tissue movement (settling)
            
            %filter the cell array by wire letter coombination
            subwire = subrat(find(~cellfun(@isempty,wirenum_matches)),:);

            %determine if there are large fluctuations in firing rate over time
            %the idea is that the same tetrode could, in principle, have
            %recorded two different interneurons over a long turn span, and
            %this might show up as a gradual decrease in the firing rate of the
            %first interneuron followed by a gradual increase of the second
            %interneuron. It's always possible the tetrode somehow moved up
            %again and re-recorded the same interneuron, but it is unlikely
            %that the interneuron would have remained so healthy, is the
            %assumption
            frate_list = [subwire{:,5}];
            mean_frate = mean(frate_list);
            mean_sq_diffs = (frate_list - mean_frate).^2/length(frate_list);
            [pks,locs] = findpeaks(mean_sq_diffs);
            
            if length(pks) == 1 
               
                %there might be two interneurons so:
                %split the subwire data in two, at the locs index location
                put_nrn1_frates = [subwire{1:locs,5}];
                put_nrn2_frates = [subwire{locs+1:end,5}];

                if length(put_nrn1_frates) > 2 && length(put_nrn2_frates) > 2
    
                    highest_fr_recording1 = subwire(put_nrn1_frates == max(put_nrn1_frates),:);
                    tmp_cell_array = [tmp_cell_array; highest_fr_recording1];
        
        
                    highest_fr_recording2 = subwire(locs+find(put_nrn2_frates == max(put_nrn2_frates)),:);
                    tmp_cell_array = [tmp_cell_array; highest_fr_recording2];

                else
                    %if else is true, then the identified local peak is only 
                    %likely separating something like the following frate
                    %list = [30 31 29 10 13] -- where the smallest
                    %identified value (the output of local peak finding) is
                    %at 10, and it doesn't seem meaningful to have later
                    %increased to 13 Hz (e.g. still likely the same neuron
                    %in the process of being lost by the tetrode, with some
                    %variation in recording noise or true spiking)

                    %so just grab the highest firing rate recording
                    highest_fr_recording =  subwire(frate_list == max(frate_list),:);
                    tmp_cell_array = [tmp_cell_array; highest_fr_recording];


                end
    
            else
                %if there are a ton of local peaks, this might be an
                %indication of firing rate variations that correspond to a 
                %single neuron being recorded on super noisy days sometimes
                %and great, low-noise days other times, so just grab the
                %highest firing rate day.
                highest_fr_recording = subwire(frate_list == max(frate_list),:);
                tmp_cell_array = [tmp_cell_array; highest_fr_recording];
    
    
            end
        end

    end
    
end

unique_ish_cells = tmp_cell_array;

    %here, do a second pass to trim out units that might be repeats but
    %were, for cluster cutting reasons, assigned a slightly different unit
    %letter even though they were cut on the same tetrode.

    %since this is the second pass, we will reset the minimum acceptable
    %distance of days, assuming our first pass was quite rigorous--this
    %time around, we will be looking wire by wire (not wire-unit
    %combinations), so there may have been enough days between recordings
    %for a wire-unit combo, but a different wire-unit may have been
    %selected in between those other two units' days and would throw off
    %our estimate of days_distance between recordings. 
    % e.g. 
    %{'LH9'}    {[ 4]}    {[13]}    {'a'}    {[45.5565]}    {'13a'}
    %{'LH9'}    {[ 8]}    {[13]}    {'b'}    {[42.2210]}    {'13b'}
    %{'LH9'}    {[15]}    {[13]}    {'a'}    {[38.7408]}    {'13a'}

    %LH9's 13a and 13b often appear simultaneously in recordings across
    %multiple days, but we would have chosen only the highest firing rate
    %version of each of them, so it feels too conservative to drop 2/3 of
    %these, and we would have kept 13a on D4 and 13a on D15 because these
    %days are actually beyond the 10-day minimum acceptable distance. 

    min_accept_distance = 5;

    %first, check for cells recorded on the same tetrode, on the same
    %day--these are likely unique neurons (unless they were artificially
    %separated, but we can't fix that at this stage). 

cleaner_tmp_array = [];
for r = 1:length(list_of_rats)

    rat = list_of_rats{r};
    keep_ind = regexp(unique_ish_cells(:,1),rat,'match');
    keep_ind = find(~cellfun(@isempty,keep_ind));

    %filter the cell array for just the rows for this rat
    subrat = unique_ish_cells(keep_ind,:);

    unique_wires = unique([subrat{:,3}]);

    %filter the data frame for the unique wires 
   
    for w = 1:length(unique_wires)
        wirenum = unique_wires(:,w);
        
        %find the number of perfect matches in subrat
        subwire = subrat([subrat{:,3}] == wirenum, :);
        subwire = sortrows(subwire,2);
        n_appearances = size(subwire,1);

        days_distance = diff([subwire{:,2}]); %need to keep recordings 
                                              %of different units made
                                              %on the same day! (e.g. == 0)
        if n_appearances >= 3
        
            if all(days_distance < min_accept_distance) 
                %if all the day distances are too short, just grab the
                %highest firing rate recording for the tetrode, to be
                %maximally conservative, since the unit character change
                %might only reflect a different number of pyramidal cells
                %co-recorded on the tetrode, such that the interneuron was
                %cut as a different unit character on different days
                frate_list = [subwire{:,5}];
                highest_fr_recording =  subwire(frate_list == max(frate_list),:);
                cleaner_tmp_array = [cleaner_tmp_array;highest_fr_recording];

            elseif any(days_distance > min_accept_distance)
                %if any of the day distances is acceptable, find which
                %one(s) it(they) are and keep the earliest day's index for
                %sure, and then you'll have to see if the index+1 has any
                %other distances from other recording days to consider

                keep_ind = find(days_distance > min_accept_distance) ;
                
                if days_distance(end) > min_accept_distance
                    %if the last index in days_distance is greater than the
                    %minimum acceptable distance, also add the last index
                    %of subwire into keep_ind, because the first keep_ind
                    %specification will only keep the next-to-last index in
                    %subwire
                    keep_ind = [keep_ind, keep_ind(end)+1];
                end

                if any(days_distance == 0)
                    zero_inds = find(days_distance == 0) + 1;%grab the last
                                                             %session with
                                                             %units on the
                                                             %same wire,
                                                             %same day
                    keep_ind = [keep_ind, zero_inds];

                end
                
                cleaner_tmp_array = [cleaner_tmp_array;subwire(keep_ind,:)];

                other_inds = find(days_distance < min_accept_distance);%this
                                                            %also catches
                                                            %any zero-index
                                                            %ones

                new_subwire = subwire(other_inds,:);
                %and from this one, just grab the highest firing rate
                %recording!
                frate_list = [new_subwire{:,5}];
                highest_fr_recording = new_subwire(frate_list == max(frate_list),:);
                cleaner_tmp_array = [cleaner_tmp_array; highest_fr_recording];

            end

        elseif n_appearances == 2

            if (days_distance < min_accept_distance) & (days_distance > 0)
    
                frate_list = [subwire{:,5}];
                highest_fr_recording = subwire(frate_list == max(frate_list),:);
                cleaner_tmp_array = [cleaner_tmp_array; highest_fr_recording];
    
            else
    
                cleaner_tmp_array = [cleaner_tmp_array; subwire];
    
            end

        else
    
            cleaner_tmp_array = [cleaner_tmp_array;subwire];
    
            
    
        end

    end
end


%now, clean the cell array back up--second column needs to have the
%'D' added to whatever the session number is, and the third column needs to
%have the 'TETSPK' appended to the first **wire** number and re-append the
%letter to the end of it

cleaner_tmp_array = sortrows(cleaner_tmp_array,[1,2]);
for i = 1:length(cleaner_tmp_array)

    unique_cells{i,1} = [cleaner_tmp_array{i,1},'D',num2str(cleaner_tmp_array{i,2})];

    if cleaner_tmp_array{i,3} < 10
        unique_cells{i,2} = ['TETSPK0',num2str(cleaner_tmp_array{i,3}),cleaner_tmp_array{i,4}];
    else
        unique_cells{i,2} = ['TETSPK',num2str(cleaner_tmp_array{i,3}),cleaner_tmp_array{i,4}];
    end

    unique_cells{i,3} = cleaner_tmp_array{i,5};

end


