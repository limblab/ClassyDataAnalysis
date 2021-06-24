function [] = loadRawMarkerDataDLC(cds,marker_data_path)
% Given a cds, marker data from color tracking, and a save location, will
% spatiotemporally align the markers to the data in the CDS, transform the coordinates for use in OpenSim and load into cds
%
% 
%% load marker data file
marker_table = readtable(marker_data_path);

% remove bad frames from maker_table
marker_table = marker_table(marker_table.is_good_frame==1,:);

%% get marker names
% get frame_num and remove from table
frame_nums = marker_table.fnum;
marker_table.fnum = [];
col_names = marker_table.Properties.VariableNames;
% remove suffixes
for i_col = 1:numel(col_names)
    underscore_idx = strfind(col_names{i_col},'_');
    col_names{i_col} = col_names{i_col}(1:underscore_idx-1);
end
marker_names = uniquecell(col_names);


%% get cds time for each frame
% find sync line. Likely called videosync
a_idx = 0;
sync_name = 'videosync';
for i_analog = 1:numel(cds.analog)
    if(any(strcmpi(cds.analog{i_analog}.Properties.VariableNames,sync_name)))
        a_idx = i_analog;
    end
end

% get frame_times based on videosync
frame_times = cds.analog{a_idx}.t(diff((cds.analog{a_idx}.(sync_name) - mean(cds.analog{a_idx}.(sync_name)) > 100)) > 0.5);

% remove first frame if it occurs > 1s before others
if(frame_times(2)-frame_times(1) > 1)
    frame_times = frame_times(2:end);
end

if(numel(frame_times) ~= numel(frame_nums))
    error('frame nums and frame times do not have the same length. Worth checking into...');
end


%% make marker table to put in cds
% currently, we are just throwing the whole table into the cds, figuring
% that the excess entries (score, n_cams) might be interesting to someone
% at some point...

marker_table.t = frame_times;

%% 1. PUT KINECT DATA INTO OPENSIM COORDINATES

% [marker_table,~] = transformForOpenSimDLC(marker_table,cds);


%% add to CDS (passed by reference)
%append new data into the analog cell array:
%stick the data in a new cell at the end of the cds.analog
%cell array:

cds.analog{end+1}=marker_table;

% set new data window
cds.setDataWindow()
cds.sanitizeTimeWindows();

logStruct=struct('fileName',marker_data_path);
evntData=loggingListenerEventData('loadRawMarkerData',logStruct);
notify(cds,'ranOperation',evntData)

