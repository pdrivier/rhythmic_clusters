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

        if n_appearances > 1
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
      
        %then find the difference between the session days
        days_distance = diff([subwire{:,2}]);

        %grab the indices corresponding to any pair of days that is under
        %the minimum threshold of days specified; remember these will be
        %indeces for subwire cell array
        danger_inds = find(days_distance < min_accept_distance);

        for d = 1:length(danger_inds)

            pair1 = danger_inds(d);
            pair2 = danger_inds(d) + 1;

            %START HERE next work session: it's not clear that this code
            %works also for more than 2 consecutive days (e.g. a scenario
            %may arise where the offending wire number combination has been
            %recorded four days in a row--the current code doesn't account
            %for this, but it needs to! Another issue with the code is that
            %the last six cells in the dataset don't have a firing rate
            %attached to them--perhaps it's worth going back to take care
            %of this first in another snippet of code that searches out the
            %corresponding file in Rat struct? remember, these wouldn't be
            %saved in the 'frate' subfield of the Rat struct because I cut
            %them later!

            %select the first and second days of the offending pair from
            %the cell array and select only the index corresponding to the
            %highest firing rate
            consider = subwire([pair1,pair2],:);
            ind_highest_fr = find([consider{:,5}]==max(consider{:,5}));
            tmp_cell_array = [tmp_cell_array; subwire(ind_highest_fr,:)]; 

        end

        end

    end

    
   
end
%remember to clean the cell array back up--second column needs to have the
%'D' added to whatever the session number is, and the third column needs to
%have the 'TETSPK' appended to the first **wire** number and re-append the
%letter to the end of it


end