function unitsFromNEV(cds,opts)
    %takes a cds handle an NEVNSx object and an options structure and
    %populates the units field of the cds
    unitList = unique([cds.NEV.Data.Spikes.Electrode;cds.NEV.Data.Spikes.Unit]','rows');
    %establish id table
    if isfield(opts,'array')
        array=opts.array;
    else
        warning('unitsFromNEV:noArrayName','The user did not specify an array name for this data. Please re-load the data and specify an array name')
        cds.addProblem('arrayNameUnknown: the user did not specify a name for the array this data comes from')
        array='?';
    end
    if isfield(opts,'monkey')
        monkey=opts.monkey;
    else
        warning('unitsFromNEV:noArrayName','The user did not specify an monkey name for this data. Please re-load the data and specify an array name')
        cds.addProblem('monkeyNameUnknown: the user did not specify a name for the monkey this data comes from')
        monkey='?';
    end
    %if we already have unit data, check that our new units come from a
        %different source so that we don't get duplicate entries
    if ~isempty(cds.units) && ~isempty(unitList) && ~isempty(find(strcmp({cds.units.array},opts.array),1,'first'))
        error('unitsFromNEV:sameArrayName','the cds and the current data have the same array name, which will result in duplicate entries in the units field. Re-load one of the data files using a different array name to avoid this problem')
    end
    %try loading a mapfile:
    noMap=true;
    if isfield(opts,'mapfile') && ~isempty(opts.mapFile)
        try
            arrayMap={};
            fid=fopen(opts.mapFile,'r');
            inheader=true;
            while ~feof(fid)
                tline=fgets(fid);
                %skip header lines:
                if inheader && ~isempty(strfind(tline,'Cerebus mapping'))
                    %this is our last header line, set the flag false and
                    %continue to the next line:
                    inheader=false;
                    continue
                elseif inheader
                    continue
                end
                if strcmp(tline(1:2),'//')
                    %this should be our column labels
                    colLabels=splitstr(tline(3:end));
                    
                    rowNumCol=strcmp(colLabels,'row');
                    colNumCol=trcmp(colLabels,'col');
                    bankCol=strcmp(colLabels,'bank');
                    pinCol=strcmp(colLabels,'elec');
                    labelCol=strcmp(colLabels,'label');
                    continue
                end
                %if we got to this point we are on an actual data line:
                tmp=textscan(tline,'%s');
                tmp=tmp{1};
                rowNum=num2str(tmp{rowNumCol});
                colNum=num2str(tmp{colNumCol});
                pin=num2str(tmp{pinCol});
                bank=char(tmp{bankCol});
                label=chan(tmp{labelCol});
                switch bank
                    case 'A'
                        chan=pin;
                    case 'B'
                        chan=pin+32;
                    case 'C'
                        chan=pin+64;
                    case 'D'
                        chan=pin+96;
                    otherwise
                        error('unitsFromNEV:badBankLabel',['unitsFromNEV is not configured to handle arrays with bank label: ',bank])
                end
                arrayMap=[arrayMap;{chan,pin,rowNum,colNum,bank,label}];
            end
            arrayMap=cell2table(arrayMap,'VariableNames',{'chan','pin','row','col','bank','label'});
            noMap=false;
        catch ME
            noMap=true;
        end
    end
    if noMap
        if exist('ME','var')
            problemData.description='tried to load mapfile and failed';
            problemData.error=ME;
        else
            problemData.description='no map file was passed';
        end
        cds.addProblem('No Map file. Electrode locations and bank ID are not available in the units structure',problemData);
    end
    
    %initialize struct array:
    units=struct('chan',cell(numel(unitList),0),...
                            'ID',cell(numel(unitList),0),...
                            'rowNum',cell(numel(unitList),0),...
                            'colNum',cell(numel(unitList),0),...
                            'pinNum',cell(numel(unitList),0),...
                            'bank',cell(numel(unitList),0),...
                            'label',cell(numel(unitList),0),...
                            'array',cell(numel(unitList),0),...
                            'wellSorted',cell(numel(unitList),0),...this is a stub as testSorting can't be run till the whole units field is populated
                            'monkey',cell(numel(unitList),0),...
                            'spikes',repmat( cell2table(cell({0,2}),'VariableNames',{'ts','wave'}),numel(unitList),1));
    %loop through and unit entries for each unit
    for i = 1:size(unitList,1)
        %we are avoiding using the set methor here in order to avoid
        %unnecessary duplication of data in memory.
%         cds.units(i)=struct('chan',unitList(i,1),...
%                             'ID',unitList(i,2),...
%                             'array',array,...
%                             'wellSorted',false,...this is a stub as testSorting can't be run till the whole units field is populated
%                             'monkey',monkey,...
%                             'spikes',table(...timestamps for current unit from the NEV:
%                                      [double(cds.NEV.Data.Spikes.TimeStamp(cds.NEV.Data.Spikes.Electrode==unitList(i,1) & ...
%                                         cds.NEV.Data.Spikes.Unit==unitList(i,2)))/30000]',... 
%                                     ...waves for the current unit from the NEV:    
%                                     double(cds.NEV.Data.Spikes.Waveform(:,cds.NEV.Data.Spikes.Electrode==unitList(i,1) ...
%                                     &  cds.NEV.Data.Spikes.Unit==unitList(i,2))'),...
%                                     'VariableNames',{'ts','wave'}));
        units(i).chan=unitList(i,1);
        if noMap
            units(i).rowNum=nan;
            units(i).colNum=nan;
            units(i).pinNum=nan;
            units(i).bank=nan;
            units(i).label=nan;
        else
            %find the correct row of our arrayMap:
            idx=find(arrayMap.chan==units(i).chan,1);
            %copy data from the map into the current unit entry:
            units(i).rowNum=arrayMap.row(idx);
            units(i).colNum=arrayMap.col(idx);
            units(i).pinNum=arrayMap.pin(idx);
            units(i).bank=arrayMap.bank(idx);
            units(i).label=arrayMap.label(idx);
        end
        units(i).ID=unitList(i,2);
        units(i).array=array;
        units(i).wellSorted=false;
        units(i).monkey=monkey;
        units(i).spikes=table(...timestamps for current unit from the NEV:
                                     [double(cds.NEV.Data.Spikes.TimeStamp(cds.NEV.Data.Spikes.Electrode==unitList(i,1) & ...
                                        cds.NEV.Data.Spikes.Unit==unitList(i,2)))/30000]',... 
                                    ...waves for the current unit from the NEV:    
                                    double(cds.NEV.Data.Spikes.Waveform(:,cds.NEV.Data.Spikes.Electrode==unitList(i,1) ...
                                    &  cds.NEV.Data.Spikes.Unit==unitList(i,2))'),...
                                    'VariableNames',{'ts','wave'});
        %check for resets in time vector
        idx=cds.skipResets(units(i).spikes.ts);
        if ~isempty(idx) && idx>1
            %if there were resets, remove everything before the resets
            units(i).spikes{1:idx,:}=[];
        end
        
    end
    cds.units=[cds.units;units];
    
%    unitscds.testSorting; %tests each sorted unit to see if it is well-separated from background and other units on the same channel
    opData.array=array;
    opData.numUnitsAdded=size(unitList,1);
    evntData=loggingListenerEventData('unitsFromNEV',opData);
    notify(cds,'ranOperation',evntData)
end