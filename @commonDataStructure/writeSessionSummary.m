function writeSessionSummary(cds) 
    %this is a method function for the common_data_structure (cds) class, and
    %should be located in a folder '@common_data_structure' with the class
    %definition file and other method files
    %
    %write_session_summary(cds)
    %this method parses the meta field of the cds and writes a flat text
    %file for parsing by humans and scripts. It will also write a line to
    %the limblabDataHistory.xlsx file containing all the same information
    %so that people can sort/organize the data in excel
    folderpath=[filesep, filesep,'fsmresfiles.fsm.northwestern.edu',filesep,'fsmresfiles',filesep,'Basic_sciences',filesep,'Phys','L_MillerLab',filesep,'data',filesep,'cds_archive',filesep];
    fnameDH=[folderpath,'limblabDataHistory.xlsx'];
    fnameSummary=[folderpath,'summary_files',filesep,cds.meta.cdsName,'.txt'];
    fnameBackup=[folderpath,'historyBackupMatlab.mat'];
    %build list of labels that we will want to write:
    dataLabels={'Source File',...
                    'dateTime',...
                    'cdsName',...
                    'processedTime',...
                    'cds version',...
                    'monkey name',...
                    'array location',...
                    'task',...
                    'collected by',...
                    'hasUnits',...
                    'hasKinematics',...
                    'hasForce',...
                    'hasEmg',...
                    'hasLfp',...
                    'hasAnalog',...
                    'hasTriggers',...
                    'hasBumps',...
                    'hasChaoticLoad',...
                    'hasSorting',...
                    'numSorted',...
                    'numWellSorted',...
                    'numDualUnits',...
                    'numTrials',...
                    'numReward',...
                    'numAbort',...
                    'numFail',...
                    'numIncomplete',...
                    'percentStill',...
                    };
    
    %open the limblabDataHistory.xlsx file and find either the line with a
    %prior entry for this data, or the first empyt line of the workbook:
    [~,~,dataHistory]=xlsread(fnameDH,'dataHistory');
    %load backupDH from mat file:
    if isempty(dir(fnameBackup))
        backupDH=[];
    else
        load(fnameBackup);
    end
    if ~isequal(backupDH,dataHistory)
        reWriteFlag=true;
    else
        reWriteFlag=false;
    end
    colNames=dataHistory(1,:);
    excelLineNum=find(strcmp(cds.meat.cdsName,dataHistory(:,strcmp('cdsName',colNames))));
    excelData=cell(1,length(colNames));
    if ~isempty(excelLineNum)
        excelLineNum=size(dataHistory,1)+1;
    else
        warning('writeSessionSummary:cdsHistoryExists','summary data for a cds with this name already exists. That summary will be overwritten')
    end
    %open the text file for writing:
    fhandle=fopen(fnameSummary,'w');
    %loop through the dataLabels and write data for each to file. Also put
    %into cell array for writing to excel file:
    for i=1:numel(dataLabels)
        itemName=dataLabels{i};
        itemData=cds.meta.(itemName);
        if ischar(itemData)
            fprintf(fhandle,'%s:\t%s\n\r',itemName,itemData);
        else
            fprintf(fhandle,'%s:\t%s\n\r',itemName,itemData);
        end
        excelData{strcmp(itemName,colNames)}=itemData;
    end
    
    %write line to limblabDataHistory.xlsx:
    if reWriteFlag
        status=xlswrite(fnameDH,[dataHistory;excelData],'dataHistory');
    else
        status=xlswrite(fnameDH,excelData,'dataHistory',excelLineNum);
    end
    if ~status
        warning('surrary data not written to limblabDataHistory.xlsx')
    end    
    
    evntData=loggingListenerEventData('writeSessionSummary',[]);
    notify(cds,'ranOperation',evntData)
end