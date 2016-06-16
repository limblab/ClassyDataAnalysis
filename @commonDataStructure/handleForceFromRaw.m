function handleForce=handleForceFromRaw(cds,raw_force,opts)
    %this is a method function of the commonDataStructure class and must be
    %saved in the @commonDataStructure folder with the class definition
    %
    % computes the handle force using the data in the NEVNSx object and
    % puts it into the correct field of the cds. since the cds is a member
    % of the handle superclass no variable needs to be passed back. Assumes
    % EXACTLY 6 channels will be flagged as ForceHandle channels. If there
    % less than 6, there will be columns of zeros returned. If more than 6
    % columns exist, names with duplication will be ignored, and channels
    % after the first 6 will be ignored silently
    %
    %There is no particular reason this needs to be a function, but it 
    %cleans up the main code to move it to this sub-function, and local
    %variables will be automatically cleared, saving memory
        
    %calculate offsets for the load cell and remove them from the force:
    if sum(cds.kin.still) > 100  % Only use still data if there are more than 100 movement free samples                
        force_offsets = mean(raw_force(cds.kin.still,:));
    else
        %issue warning
        warning('NEVNSx2cds:noStillTime','Could not find 100 points of still time to compute load cell offsets. Defaulting to mean of force data')
        %make known problem entry
        cds.addProblem('No still data to use for computing load cell offsets. Offsets computed as mean of all load cell data')
        force_offsets = mean(raw_force);
    end

    % Get calibration parameters based on lab number            
    if isfield(opts,'labNum') && opts.labNum>0
        [fhcal,rotcal,Fy_invert]=getLabParams(opts.labNum,opts.dateTime,opts.rothandle);
    else
        error('handleForceFromRaw:LabNotSet','handleForceFromRaw needs the lab number in order to select the correct load cell calibration')
    end
    raw_force = (raw_force -  repmat(force_offsets, length(raw_force), 1)) * fhcal * rotcal;
    clear force_offsets;

    % fix left hand coords in some force data
    raw_force(:,2) = Fy_invert.*raw_force(:,2); 

    %rotate load cell data into room coordinates using robot arm
    %angle
    if size(raw_force,2)==2
        if isfield(opts,'labNum')&& opts.labNum==3 %If lab3 was used for data collection  
            handleForce=table( raw_force(:,1).*cos(-cds.enc.th2) - raw_force(:,2).*sin(cds.enc.th2),...
                raw_force(:,1).*sin(cds.enc.th2) + raw_force(:,2).*cos(cds.enc.th2),...
                'VariableNames',{'fx','fy'});
        elseif isfield(opts,'labNum')&& opts.labNum==6 %If lab6 was used for data collection         
            if datenum(opts.dateTime)<datenum('07-Mar-2016')
                %the old lab6 arm was on the opposite side of the robot
                %from the other labs and the new lab6 setup. This handles
                %that condition
                handleForce=table( raw_force(:,1).*cos(-cds.enc.th1) - raw_force(:,2).*sin(cds.enc.th1),...
                                    raw_force(:,1).*sin(cds.enc.th1) + raw_force(:,2).*cos(cds.enc.th1),...
                                    'VariableNames',{'fx','fy'});
            else
                handleForce=table( raw_force(:,1).*cos(-cds.enc.th2) - raw_force(:,2).*sin(cds.enc.th2),...
                                    raw_force(:,1).*sin(cds.enc.th2) + raw_force(:,2).*cos(cds.enc.th2),...
                                    'VariableNames',{'fx','fy'});
            end
        end
    elseif size(raw_force,2)==6
        if isfield(opts,'labNum')&& opts.labNum==3 %If lab3 was used for data collection  
            handleForce=table( raw_force(:,1).*cos(-cds.enc.th2) - raw_force(:,2).*sin(cds.enc.th2),...
                raw_force(:,1).*sin(cds.enc.th2) + raw_force(:,2).*cos(cds.enc.th2),...
                raw_force(:,3),...
                raw_force(:,4).*cos(-cds.enc.th2) - raw_force(:,5).*sin(cds.enc.th2),...
                raw_force(:,4).*sin(cds.enc.th2) + raw_force(:,5).*cos(cds.enc.th2),...
                raw_force(:,6),...
                'VariableNames',{'fx','fy','fz','mx','my','mz'});
        elseif isfield(opts,'labNum')&& opts.labNum==6 %If lab6 was used for data collection         
            if datenum(opts.dateTime)<datenum('07-Mar-2016')
                %the old lab6 arm was on the opposite side of the robot
                %from the other labs and the new lab6 setup. This handles
                %that condition
                handleForce=table( raw_force(:,1).*cos(-cds.enc.th1) - raw_force(:,2).*sin(cds.enc.th1),...
                                    raw_force(:,1).*sin(cds.enc.th1) + raw_force(:,2).*cos(cds.enc.th1),...
                                    raw_force(:,3),...
                                    raw_force(:,4).*cos(-cds.enc.th1) - raw_force(:,5).*sin(cds.enc.th1),...
                                    raw_force(:,4).*sin(cds.enc.th1) + raw_force(:,5).*cos(cds.enc.th1),...
                                    raw_force(:,6),...
                                    'VariableNames',{'fx','fy','fz','mx','my','mz'});
            else
                handleForce=table( raw_force(:,1).*cos(-cds.enc.th2) - raw_force(:,2).*sin(cds.enc.th2),...
                                    raw_force(:,1).*sin(cds.enc.th2) + raw_force(:,2).*cos(cds.enc.th2),...
                                    raw_force(:,3),...
                                    raw_force(:,4).*cos(-cds.enc.th2) - raw_force(:,5).*sin(cds.enc.th2),...
                                    raw_force(:,4).*sin(cds.enc.th2) + raw_force(:,5).*cos(cds.enc.th2),...
                                    raw_force(:,6),...
                                    'VariableNames',{'fx','fy','fz','mx','my','mz'});
            end
        end
    else
        error('handleForceFromRaw:BadConversion',['Expected either 2 or 6 channels in converted handle force. Instead got ',num2str(size(raw_force,2)),' channels'])
    end
    
end
