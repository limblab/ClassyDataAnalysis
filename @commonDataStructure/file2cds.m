function file2cds(cds,filePath,varargin)
    %this is a method function for the common_data_structure (cds) class, and
    %should be located in a folder '@common_data_structure' with the class
    %definition file and other method files
    %
    %file2cds(folderPath,fileName) loads the file(s) specified into
    %temporary fields in the cds, processes the raw data into the cds
    %format, and deletes the raw data
    %file2cds(filePath,options...)
    %Accepts the following options as additional input arguments:
    %'rothandle'    assumes the handle was flipped upside down, and
    %               converts force signals appropriately
    %'ignoreJumps' does not flag jumps in encoder output when converting
    %               encoder data to position. Useful for Lab1 where the
    %               encoder stream is used to process force, not angle
    %               steps. 
    %'ignoreFilecat'   Ignores the join time between two files. Normally
    %               the join time is flagged and included in the list of
    %               bad kinematic times so that artifacts associated with
    %               the kinematic discontinuity can be avoided in further
    %               processing
    %lab number:    an integer number designating the lab from the set 
    %               1,2,3,6
    %'taskTASKNAME' specifies the task performed during data collection.
    %               NEVNSx2cds looks for the first part of the argument to
    %               match the string 'task' and then takes the remainder of
    %               the string to be the task name. for example 'taskRW'
    %               would result in the task being set to 'RW'. Currently
    %               viable task strings are:
    %                   CO: center out
    %                   CObump: center out bump task
    %                   RW: random walk
    %                   BD: bump direction
    %                   WF: wrist flexion
    %                   multi_gadget: multigadget task
    %                   UNT: uncertainty
    %                   RP:resist perturbations
    %                   DCO: dynamic center out
    %'arrayARRAYNAME'   specifies the array used for data collection.
    %               file2cds looks for the first part of the argument to
    %               match the string 'array' and then takes the remainder of
    %               the string to be the array name. for example 'arrayM1'
    %               would result in the task being set to 'M1'
    %'monkeyMONKEYNAME' specifies the monkey that this data is from.
    %               file2cds looks for the first part of the argument to
    %               match the string 'monkey' and then takes the remainder of
    %               the string to be the monkey name. for example 'monkeyChips'
    %               would result in the task being set to 'Chips'
    %'ranByPERSONNAME'  specifies the person who handled the monkey for the
    %               experiment. file2cds looks for the first part of the
    %               argument to match the string 'ranBy' and then takes the
    %               remainder of the string as the person's name
    %'mapFileFILEPATH'  specifies the path to the correct mapfile for this
    %               array. file2cds looks for the first part fo the
    %               arguemtn to match the string 'mapFile' and then takes
    %               the remainder of the string as the full file path of
    %               the map file.
    %example: cds.file2cds('C:/datafolder/data.nev', 'rothandle', 3,'taskCO',
    %'arrayM1','monkeyChips','ranByTucker','mapFileC:/datafolder/map.cmp')
    %imports the data from data.nev and data.nsx into the fields of cds, 
    %assuming the robot handle was inverted, the data came from lab3, the 
    %task was the center out task, the array was in M1, the monkey was
    %Chips, Tucker conducted the experiment and the map file can be found
    %in map.cmp.
    
    
    %% construct opts structure:
    opts=struct('labNum',-1,'rothandle',false,'ignore_jumps',false,'ignore_filecat',false,'robot',false,'task','Unknown','hasChaoticLoad',false); 

        % Parse arguments
        if ~isempty(varargin)
            for i = 1:length(varargin)
                optStr = char(varargin{i});           
                if strcmp(optStr, 'rothandle')
                    opts.rothandle = true;
                elseif strcmp(optStr, 'chaoticLoad')
                    opts.hasChaoticLoad=true;
                elseif strcmp(optStr, 'ignoreJumps')
                    opts.ignore_jumps=true;
                elseif strcmp(optStr, 'ignoreFilecat')
                    opts.ignore_filecat=true;
                elseif ischar(optStr) && length(optStr)>4 && strcmp(optStr(1:4),'task')
                    opts.task=optStr(5:end);
                elseif ischar(optStr) && length(optStr)>5 && strcmp(optStr(1:5),'array')
                    opts.array=optStr(6:end);
                elseif ischar(optStr) && length(optStr)>5 && strcmp(optStr(1:6),'monkey')
                    opts.monkey=optStr(7:end);
                elseif ischar(optStr) && length(optStr)>5 && strcmp(optStr(1:5),'ranBy')
                    opts.ranBy=optStr(6:end);
                elseif ischar(optStr) && length(optStr)>5 && strcmp(optStr(1:5),'mapFile')
                    opts.mapFile=optStr(8:end);
                elseif isnumeric(varargin{i})
                    opts.labNum=varargin{i};    %Allow entering of the lab number               
                else 
                    error('Unrecognized option: %s', optStr);
                end
            end
        end
        %check the options and throw warnings if some things aren't set:
        flag=0;
        if strcmp(opts.task,'Unknown')
            flag=true;
            warning('NEVNSx2cds:taskNotSet','No task was passed as an input variable. Further processing can attempt to automatically identify the task, but success is not garaunteed')
        end
        if ~isfield(opts,'array')
            flag=true;
            warning('NEVNSx2cds:arrayNotSet','No array label was passed as an input variable.')
        end
        if opts.labNum==-1
            flag=true;
            warning('NEVNSx2cds:labNotSet','The lab number where this data was collected was not passed as an input variable')
        end
        if ~isfield(opts,'monkey')
            flag=true;
            warning('NEVNSx2cds:monkeyNotSet','The monkey from which this data was collected was not passed as an input variable')
        end
        if ~isfield(opts,'ranBy')
           flag=true;
           warning('NEVNSx2cds:monkeyNotSet','The person who collected this data was not passed as an input variable')
        end
        if ~isfield(opts,'mapFile')
            flag=true;
            warning('NEVNSx2cds:noMapFilePassed','The path to the array map file was not passed as an input variable')
        end
        if flag
            while 1
                s=input('Do you want to cancel and re-run this data load including the information missing above? (y/n)\n','s');
                if strcmpi(s,'y')
                    error('NEVNSx2cds:UserCancelled','User cancelled execution to re-run with additional input')
                elseif strcmpi(s,'n')
                    break
                else
                    disp([s,' is not a valid response'])
                end
            end
        end
        %set the robot flag if we are using one of the robot labs:
        if opts.labNum == 2 || opts.labNum == 3 || opts.labNum ==6
            opts.robot=true;
        end
    
    
%     %see if the 'noDB' flag was passed
%     noDB=0;
%     for i=1:length(varargin)
%         if strcmp(varargin{i},'noDB')
%             noDB=1;
%         end
%     end
%     %check the database to see if the file has already been processed and
%     %load it if it has:
%     dataFromDB=0;
%     if ~noDB
%         try
%             conn=database('LLTestingDB','LLMatlabScript','mvemjlht','Vendor','PostgreSQL','Server','vfsmmillerwiki.fsm.northwestern.edu');
%             data=fetch(conn,['select *** from session where sourceFile = ',fileName]
%             if ~isempty(data)
%                 %load from database
%                 cds.database2cds(conn,data,varargin)
%                 dataFromDB=1;
%             end
%             close(conn);
%         catch 
%             warning('file2cds:databaseError','Failed to connect or fetch data from the LLSessionsDB database')
%             
%         end
%     end
%     if ~dataFromDB
%         varargin=[varargin,{'dbEmpty'}];
        cds.nev2NEVNSx(filePath);
        cds.NEVNSx2cds(opts);
        cds.clearTempFields()
%        cds.writeSessionSummary()
        evntData=loggingListenerEventData('file2cds',[]);
        notify(cds,'ranOperation',evntData)
%     end
end