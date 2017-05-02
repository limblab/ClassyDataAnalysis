function loadOpenSimData(cds,folderPath)
    %this is a method of the cds class and should be stored in the
    %@commonDataStructure folder with the other class methods.
    %
    %attempts to load Open Sim data from the source directory of the cds.
    %This uses the meta field to try and find properly named files in the
    %folder specified. The prefix of the file myst match the name of the
    %source file of the cds
    
    
    if ~strcmp(folderPath(end),filesep)
        folderPath=[folderPath,filesep];
    end
    
    %interpolate onto times aligned with the existing kinematic data. If no
    %existing kinematic data, just use a 100hz signal aligned to zero.
    if ~isempty(cds.kin)
        dt=mode(diff(cds.kin.t));
    else
        dt=.01;
    end
    
    prefix=cds.meta.rawFileName;
    if ~iscell(prefix)
        prefix={prefix};
    end
    foundFiles={};
    for i=1:numel(prefix)
        %find and strip extensions if present
        extLoc=max(strfind(prefix{i},'.'));
        if ~isempty(extLoc)
            prefix{i}=prefix{i}(1:extLoc-1);
        end
        %look for joint kinematics *_Kinematics_q.sto files with matching
        %prefix:
        fileNameList={[folderPath,prefix{i},'_Kinematics_u.sto'];...
%             [folderPath,prefix{i},'_MuscleAnalysis_Length.sto'];...
            [folderPath,prefix{i},'_Dynamics_q.sto']};
        for j=1:numel(fileNameList)
            foundList=dir(fileNameList{j});
            if ~isempty(foundList)
                %load data from file into table 'kin':
                fid=fopen(fileNameList{j});
                %loop through the header till we find the first row of data:
                tmpLine=fgetl(fid);
                while ~strcmp(tmpLine,'endheader')
                    if ~isempty(strfind(tmpLine,'nRows'))
                        nRow=str2double(tmpLine(strfind(tmpLine,'=')+1:end));
                    elseif ~isempty(strfind(tmpLine,'nColumns'))
                        nCol=str2double(tmpLine(strfind(tmpLine,'=')+1:end));
                    elseif ~isempty(strfind(tmpLine,'inDegrees'))
                        if ~isempty(strfind(tmpLine,'yes'))
                            unitLabel='deg';
                        else
                            unitLabel='rad';
                        end
                    end
                    tmpLine=fgetl(fid);
                end
                header=strsplit(fgetl(fid));
                %convert 'time' to 't' to match cds format:
                idx=find(strcmp(header,'time'),1);
                if isempty(idx)
                    %look for a 't' column
                    idx=find(strcmp(header,'t'),1);
                    if isempty(idx)
                        error('loadOpenSimData:noTime',['could not find a time column in the file: ', kinFileName])
                    end
                else
                    %convert 'time into 't'
                    header{idx}='t';
                end
                a=reshape(fscanf(fid,repmat('%f',[1,nCol])),[nCol,nRow])';
                %sanity check time:
                SR=mode(diff(a(:,1)));
                if size(a,1)~=round((1+ (max(a(:,1))-min(a(:,1)))/SR))
                    warning('loadOpenSimData:badTimeSeries',['the timeseries in the detected opensim data is missing time points. expected ',num2str((1+ (max(a(:,1))-min(a(:,1)))/SR)),' points, found ',num2str(size(a,1)),' points'])
                    disp('data will be interpolated to reconstruct missing points')
                    cds.addProblem('kinect data has missing timepoints, data in the cds has been interpolated to reconstruct them')
                end
                %interpolate to desired time vector:
                desiredTime=roundTime(a(1,1):dt:a(end,1));%uniformly samples a with spacing dt, then shifts time bins to be zero aligned
                desiredTime=desiredTime(desiredTime>min(a(:,1)) & desiredTime<max(a(:,1)))';%clear out any points that fall outside the original time window due to the shift
                
                kin=array2table([desiredTime,interp1(a(:,1),a(:,2:end),desiredTime)],'VariableNames',header);
                unitsLabels=[{'s'},repmat({unitLabel},[1,nCol-1])];
                kin.Properties.VariableUnits=unitsLabels;
                %find sampling rate and look for matching rate in analog data:
                SR=round(1/mode(diff(kin.t)));
                cdsFrequencies=zeros(1,length(cds.analog));
                for k=1:length(cds.analog)
                    cdsFrequencies(k)=round(1/mode(diff(cds.analog{k}.t)));
                end
                match=find(cdsFrequencies==SR);
                %append new data into the analog cell array:
                if isempty(match)
                    %stick the data in a new cell at the end of the cds.analog
                    %cell array:
                    cds.analog{end+1}=kin;
                else
                    %append the new data to the table with the matching
                    %frequency:
                    cds.analog{match}=mergeTables(cds.analog{match},kin);
                end
                foundFiles=[foundFiles;fileNameList(j)];
            end
        end
           
    end
    
    
    cds.sanitizeTimeWindows
    logStruct=struct('folder',folderPath,'fileNames',foundFiles);
    evntData=loggingListenerEventData('loadOpenSimData',logStruct);
    notify(cds,'ranOperation',evntData)
end