%start_here_set_your_data_up

%First you'll need to make a list of interneurons you will analyze.
%Although challenging to establish without a shadow of a doubt,
%interneurons can often be recorded over multiple days, after having turned
%the tetrode down several micrometers. Whereas you would easily lose
%pyramidal cells with a tetrode turn, you would keep the interneuron in the
%vicinity (take a look at reconstructions by Lasztoczi & Klausberger 2014
%and others--the interneurons are massively layer-spanning). 

%Since we want to be as conservative as humanly possible, we will reject
%interneurons from the original list that were recorded on the same wire
%across days that are too close together--the user can set this tolerance
%if desired. Another heuristic we will use is that of firing rate
%variations across days. If there is a "neuron" whose firing rates steadily
%decrease, and then several days later suddenly increases, we might
%imagine these are potentially two different interneurons: one that we were
%turning away from and so its firing rate decreased gradually, and a new
%one that we are beginning to hit. This seems more likely than assuming the
%tetrode rose back to the original position were the highest firing rates
%were recorded for that neuron--we don't expect the tissue we've turned
%away from to be healthy like that, though it's of course possible. 

% -------------------------------------------------------------------------
%We'll first load the original list of interneuron identities--------------
% -------------------------------------------------------------------------
load('Int_Identities_Automaze.mat')

list_to_select_from = Int_Identities;
delimratsess='D';
delimtetunit = 'TETSPK';
min_accept_distance = 10;

unique_cells = select_unique_cells(list_to_select_from,delimratsess,delimtetunit,min_accept_distance);
clearvars -except unique_cells
save('unique_ints_by_daytetfrate.mat','unique_cells') %save your output!

% -------------------------------------------------------------------------
%now, run the script that makes your odor sampling intervals
% -------------------------------------------------------------------------
start_here_make_intervals_odorsamp

