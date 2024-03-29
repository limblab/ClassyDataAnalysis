function forceFromNSx(cds,opts)
    %finds analog channels in the NEVNSx with labels indicating they
    %contain force data and parses them based on the options in opts and
    %the filters in the cds. Because the cds is a member of the handle
    %superclass, this function does not return anything
    
    t=[];
    force=[];
    forces=[];
    handleforce=[];
    %forces for wf and other tasks that use force_ to denote force channels
    forceCols = find(~cellfun('isempty',strfind(lower(cds.NSxInfo.NSx_labels),'force_')));
    % Sometimes force columns cannot be found because the experimenters
    % used different names, so try to find again with a very common naming
    if isempty(forceCols)
        forceCols = find(contains(cds.NSxInfo.NSx_labels, 'F1')|...
                         contains(cds.NSxInfo.NSx_labels, 'F2')|...
                         contains(cds.NSxInfo.NSx_labels, 'Fx')|...
                         contains(cds.NSxInfo.NSx_labels, 'Force_X')|...
                         contains(cds.NSxInfo.NSx_labels, 'Force_Y')|...
                         contains(cds.NSxInfo.NSx_labels, 'Fy'));
    end
%     keyboard
    robotForceChannels = find(~cellfun('isempty',strfind(cds.NSxInfo.NSx_labels,'ForceHandle')));
    if isempty(forceCols)&&isempty(robotForceChannels)
        %if we didn't find any forces skip force processing
        return
    end
    if ~isempty(forceCols)
        [loadCellData,t]=cds.getResampledFromNSx(cds.kinFilterConfig.sampleRate,forceCols);
        %build our table of force data:
        labels=cell(1,length(forceCols));
        for i=1:length(forceCols)
            %if we have x or y force, give the field our special
            %label so that later processing can find it easily
            if strcmpi(cds.NSxInfo.NSx_labels(forceCols(i)),'force_x')||...
                    strcmpi(cds.NSxInfo.NSx_labels(forceCols(i)),'Fx')||...
                    strcmpi(cds.NSxInfo.NSx_labels(forceCols(i)), 'Force_X')||...
                    strcmpi(cds.NSxInfo.NSx_labels(forceCols(i)),'F1')
                labels(i)={'fx'};
            elseif strcmpi(cds.NSxInfo.NSx_labels(forceCols(i)),'force_y')||...
                    strcmpi(cds.NSxInfo.NSx_labels(forceCols(i)),'Fy')||...
                    strcmpi(cds.NSxInfo.NSx_labels(forceCols(i)), 'Force_Y')||...
                    strcmpi(cds.NSxInfo.NSx_labels(forceCols(i)),'F2')
                labels(i)={'fy'};
            else
                labels(i)=cds.NSxInfo.NSx_labels(forceCols(i));
            end
        end
        %truncate to deal with the fact that encoder data doesn't start
        %recording till 1 second into the file and store in a table
        t=roundTime(t,.00001);
        if opts.robot % non-robot devices don't have any encoder data
            force=array2table(loadCellData(t>=min(cds.enc.t) & t<=max(cds.enc.t),:),'VariableNames',labels);
        else
            force = array2table(loadCellData,'VariableNames',labels);
        end
    end
    %forces for robot:
    if opts.robot
        if length(robotForceChannels)==6
            if isempty(cds.enc)
                warning('forceFromNEVNSx:noEncoderAngles','Encoder data is required to compute handle forces from raw load cell inputs. 6 load cell inputs are present, but no encoder data was found. Load cell data not included in cds')
                cds.addProblem('missing encoder data: tried to load handle force data but had no encoder data to compute load direction from')
            else
                achan_index=-1*ones(1,6);
                for i=1:6
                    achan_index(i)=find(~cellfun('isempty',strfind(cds.NSxInfo.NSx_labels,['ForceHandle',num2str(i)])));
                end
                %pull filtered analog data for load cell:
                [loadCellData,t]=cds.getResampledFromNSx(cds.kinFilterConfig.sampleRate,achan_index);
                %truncate to handle the fact that encoder data doesn't start
                %recording until 1 second into the file and convert load cell 
                %voltage data into forces
                t=roundTime(t,.00001);
                handleforce=cds.handleForceFromRaw(loadCellData,t,opts);
                %write temp into the cds
                %sorry about the rounding time on the line below. there are some edge
                %cases where machine precision becomes an issue and rounding the time
                %takes care of that.
                timePrecision = 1/cds.kinFilterConfig.sampleRate;
                min_t = roundTime(max([min(t), min(cds.enc.t)]),timePrecision);
                max_t = roundTime(min([max(t), max(cds.enc.t)]),timePrecision);
                timeTable=table(roundTime(t(roundTime(t)>=min_t & roundTime(t)<=max_t)),'VariableNames',{'t'});
                forces=[timeTable,handleforce,force];
            end
        else
            handleforce=[];
            if isempty(cds.force)
                warning('forceFromNEVNSx:noForceSignal','No force handle signal found because calc_from_raw did not find 6 channels named ''ForceHandle*''');
            end
        end
    else
        forces = [array2table(t),force];
    end

    if ~isempty(forces)
        forces.Properties.VariableUnits=[{'s'} repmat({'N'},1,size(handleforce,2)+size(force,2))];
        forces.Properties.Description='a table containing force data. First column is time, all other columns will be forces. If possible forces in x and y are identified and labeled fx and fy';
    else
        forces=cell2table(cell(0,3),'VariableNames',{'t','fx','fy'});
        forces.Properties.VariableUnits=[{'s'} repmat({'N'},1,2)];
        forces.Properties.Description='an empty table. No force data was found in the data source';
    end
    
    if isempty(cds.force)
        set(cds,'force',forces);
    elseif ~isempty(force)
        set(cds,'force',mergeTables(cds.force,forces));
    end
    evntData=loggingListenerEventData('forceFromNSx',cds.kinFilterConfig);
    notify(cds,'ranOperation',evntData)
end