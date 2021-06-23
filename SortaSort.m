function varargout = SortaSort(varargin)
% SORTASORT MATLAB code for SortaSort.fig
%      SORTASORT, by itself, creates a new SORTASORT or raises the existing
%      singleton*.
%
%      H = SORTASORT returns the handle to a new SORTASORT or the handle to
%      the existing singleton*.
%
%      SORTASORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SORTASORT.M with the given input arguments.
%
%      SORTASORT('Property','Value',...) creates a new SORTASORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SortaSort_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SortaSort_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SortaSort

% Last Modified by GUIDE v2.5 18-Jun-2021 00:40:43

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @SortaSort_OpeningFcn, ...
                       'gui_OutputFcn',  @SortaSort_OutputFcn, ...
                       'gui_LayoutFcn',  [] , ...
                       'gui_Callback',   []);
    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end

    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end
    % End initialization code - DO NOT EDIT
end


% --- Executes just before SortaSort is made visible.
function SortaSort_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SortaSort (see VARARGIN)

    % Choose default command line output for SortaSort
    handles.output = hObject;
    handles.tblSort.Data = [];
    handles.TblSortInds = [];
    % Update handles structure
    guidata(hObject, handles);

    % UIWAIT makes SortaSort wait for user response (see UIRESUME)
    % uiwait(handles.figure1);
end


% --- Outputs from this function are returned to the command line.
function varargout = SortaSort_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
    if isfield(handles,'rez')
        varargout{1} = handles.rez;
    else
        varaargout{1} = [];
    end
end


% --- Executes on selection change in pmCol1.
function pmCol1_Callback(hObject, eventdata, handles)
% hObject    handle to pmCol1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pmCol1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pmCol1
    UpdateTable(hObject, handles)
end

% --- Executes during object creation, after setting all properties.
function pmCol1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pmCol1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end


% --- Executes on selection change in pmCol2.
function pmCol2_Callback(hObject, eventdata, handles)
% hObject    handle to pmCol2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pmCol2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pmCol2
    UpdateTable(hObject, handles)
end

% --- Executes during object creation, after setting all properties.
function pmCol2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pmCol2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end


% --- Executes on selection change in pmCol3.
function pmCol3_Callback(hObject, eventdata, handles)
% hObject    handle to pmCol3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pmCol3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pmCol3
    UpdateTable(hObject, handles)
end

% --- Executes during object creation, after setting all properties.
function pmCol3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pmCol3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end



% --- Executes on button press in pbMerge.
function pbMerge_Callback(hObject, eventdata, handles)
% hObject    handle to pbMerge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    selClusters = handles.SelClusters;
    rez = handles.rez;
    if length(selClusters)<2
        return;
    end
    [allSpikes, sInds] = sort(cell2mat(rez.SortaSort.GroupedSpikes(selClusters)'));
    %[~, uniqSpks] = unique(round(allSpikes/30)); 
    uniqSpks = [true; diff(allSpikes)>30]; % remove spikes that are adjacent by 1 ms
    allSpikes = allSpikes(uniqSpks);
    counts = length(allSpikes);
%     newGroupWave = sum(bsxfun(@times,rez.SortaSort.dWU(:,:,selClusters),...
%         permute(counts,[1 3 2]))/sum(counts),3);
%     rez.SortaSort.dWU(:,:,selClusters) = [];
%     rez.SortaSort.dWU = cat(3,rez.SortaSort.dWU, newGroupWave);
    rez.SortaSort.ClusterID(selClusters) = [];
    rez.SortaSort.ClusterID{end+1} = num2str(max(cellfun(@str2num,rez.SortaSort.ClusterID))+1);
    rez.SortaSort.ClusterType(selClusters) = [];
    rez.SortaSort.ClusterType(end+1) = categorical({'Unknown'},...
        {'Unknown' 'Noise' 'MU' 'SU'});
    
    rez.SortaSort.GroupedSpikes(selClusters) = [];
    rez.SortaSort.GroupedSpikes{end+1} = allSpikes;
    
    % not optimal, should find a good way to rescale amplitudes given the
    % new average spike waveforms.
    allAmps = cell2mat(rez.SortaSort.GroupedAmplitudes(selClusters)');
    allAmps = allAmps(sInds);
    allAmps = allAmps(uniqSpks);
    rez.SortaSort.GroupedAmplitudes(selClusters) = [];
    rez.SortaSort.GroupedAmplitudes{end+1} = allAmps;
    
    rez.SortaSort.SpkWaves(selClusters) = [];
    handles.rez = rez;
    rez.SortaSort.SpkWaves{end+1} = ExtractSpkWaveforms(length(rez.SortaSort.GroupedSpikes), 100, handles);

%     allPCs = cell2mat(rez.SortaSort.GroupedPCs(selClusters)');
%     allPCs = allPCs(sInds,:);
%     allPCs = allPCs(uniqSpks,:);
%     rez.SortaSort.GroupedPCs(selClusters) = [];
%     rez.SortaSort.GroupedPCs{end+1} = allPCs;
    
    rez.SortaSort.DistData = CalculateDistances(rez);
    handles.SelClusters = numel(rez.SortaSort.ClusterID);
    handles.rez = rez;
    
    guidata(hObject, handles);
    UpdateTable(hObject, handles);
    UpdatePlots(hObject, handles);
end

% --- Executes on button press in pbRenumberClusters.
function pbRenumberClusters_Callback(hObject, eventdata, handles)
% hObject    handle to pbRenumberClusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    newClusterIDs = arrayfun(@num2str, 1:numel(handles.rez.SortaSort.ClusterID), ...
        'uniformoutput',false);
    handles.rez.SortaSort.ClusterID = newClusterIDs;
    guidata(hObject, handles);
    UpdateTable(hObject, handles);
    UpdatePlots(hObject, handles);
end

% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

end

% --------------------------------------------------------------------
function menuOpen_Callback(hObject, eventdata, handles)
% hObject    handle to menuOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    if isfield(handles,'CurrFPath')
        [fName, fPath] = uigetfile([handles.CurrFPath '*.mat']);
    else
        [fName, fPath] = uigetfile('*.mat');
    end
    handles.CurrFPath = fPath;
    fullPath = fullfile(fPath,fName);
    load(fullPath);
    % create index into fil and dat data
    filProp = dir(rez.ops.fbinary);
    numTPts = filProp.bytes/(2*rez.ops.NchanTOT);
    handles.mmfFil = memmapfile(rez.ops.fbinary,'format',{'int16' [rez.ops.NchanTOT numTPts] 'file'});
    handles.mmfDat = memmapfile([rez.ops.fbinary(1:(end-3)) 'dat'],'format',{'int16' [rez.ops.NchanTOT numTPts] 'file'});
        
    if exist('rez')
        if ~isfield(rez,'SortaSort')
            rez.SortaSort.Settings = struct('uVScale', 2.343, ... % convert binary values to microvolts
                                            'uVSep', 100, ... % separation between voltage traces for trace plots
                                            'XCorrWindow', 0.1, ... % duration of the XCorr window in seconds
                                            'TraceWindow', 0.03, ... % separation between voltage traces for trace plots                                            
                                            'uVScaleBar', 25, ... % for scaling the waveforms plot
                                            'msScaleBar', 1); % for scaling the waveforms plot
            numClusters = size(rez.W,2);
%             for j = 1:numClusters
%                 rez.SortaSort.dWU(:,:,j) = squeeze(rez.W(:,j,:))*squeeze(rez.U(:,j,:))';
%             end
%             if isa(rez.SortaSort.dWU,'gpuArray')
%                 rez.SortaSort.dWU = gather(rez.SortaSort.dWU);
%             end
            % remove duplicate spikes
            [~,uniqSpks,~] = unique(rez.st3(:,[1 2]),'rows');
            rez.st3 = rez.st3(uniqSpks,:);
            %rez.cProj = rez.cProj(uniqSpks,:);
            %rez.cProjPC = rez.cProjPC(uniqSpks,:,:);
            rez.SortaSort.ClusterID = arrayfun(@num2str,1:numClusters,'uniformoutput',false);
            rez.SortaSort.ClusterType = categorical(repmat({'Unknown'},numClusters,1),...
                {'Unknown' 'Noise' 'MU' 'SU'});
            rez.SortaSort.GroupedSpikes = arrayfun(@(x)rez.st3(rez.st3(:,2)==x,1),...
                1:numClusters,'uniformoutput',false);
            rez.SortaSort.GroupedAmplitudes = arrayfun(@(x)rez.st3(rez.st3(:,2)==x,3),...
                1:numClusters,'uniformoutput',false);
            handles.rez = rez;
            for j = 1:numClusters
                rez.SortaSort.SpkWaves{j} = ExtractSpkWaveforms(j, 100, handles);
            end
            %rez.SortaSort.GroupedPCs = arrayfun(@(x)rez.cProjPC(rez.st3(:,2)==x,1:3,1),...
            %    1:numClusters,'uniformoutput',false);
            
            rez.SortaSort.DistData = CalculateDistances(rez);
        end
        handles.SelClusters = 1;
        handles.CurrPlotClusters = nan;
        handles.rez = rez;
        
        handles.figure1.Name = ['SortaSort' fPath]; 
        guidata(hObject, handles);
        UpdateTable(hObject, handles)
    else
        error('No rez data present');
    end

end

% --------------------------------------------------------------------
function menuSave_Callback(hObject, eventdata, handles)
% hObject    handle to menuSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    rez = handles.rez;
    uisave('rez',[handles.CurrFPath '*.mat']);
    disp('Done saving')
end

% --------------------------------------------------------------------
function menuRetWorkspace_Callback(hObject, eventdata, handles)
% hObject    handle to menuRetWorkspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    SortaSort_OutputFcn(hObject, eventdata, handles)
end

% --- Executes on button press in cbProfile.
function cbProfile_Callback(hObject, eventdata, handles)
% hObject    handle to cbProfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbProfile
    UpdatePlots(hObject, handles)
end

% --- Executes on button press in cbWaveform.
function cbWaveform_Callback(hObject, eventdata, handles)
% hObject    handle to cbWaveform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbWaveform
    UpdatePlots(hObject, handles)
end

% --- Executes on button press in cbXCorr.
function cbXCorr_Callback(hObject, eventdata, handles)
% hObject    handle to cbXCorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbXCorr
    UpdatePlots(hObject, handles)
end

% --- Executes on button press in cbTrace.
function cbTrace_Callback(hObject, eventdata, handles)
% hObject    handle to cbTrace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbTrace
    UpdatePlots(hObject, handles)
end

% --- Executes on button press in cbAmplitudes.
function cbAmplitudes_Callback(hObject, eventdata, handles)
% hObject    handle to cbAmplitudes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbAmplitudes
    UpdatePlots(hObject, handles)
end

% --- Executes on button press in pbSort.
function pbSort_Callback(hObject, eventdata, handles)
% hObject    handle to pbSort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    UpdateTable(hObject, handles)
end

% Calculates all the various distance and similarity metrics between units
function distData = CalculateDistances(rez)

    numClus = length(rez.SortaSort.SpkWaves);
    [numSamp, numElec, ~] = size(rez.SortaSort.SpkWaves{1});
    meanWaves = cell2mat(permute(cellfun(@(x)mean(x,3), rez.SortaSort.SpkWaves,...
                'uniformoutput',false),[3 1 2]));
    distData.MeanSpkWaves = meanWaves;
    % calculate profile map for each unit
    profMap = sqrt(squeeze(sum(meanWaves.^2,1)));
    profMap(isinf(profMap)) = 0;
    profMap = bsxfun(@rdivide,profMap,sum(profMap,1));
    distData.ProfileMap = profMap;
    
    % calculate distances between profile maps
    profMapDists = 1-corr(profMap);
    distData.ProfileMapDist = profMapDists;
    
    % FIX PEAK ELEC CHAN
    % SET XCORR BIN SIZE, update XCORR function
    % determine profile map peaks
    for j = 1:numClus
        [~, peakInd(j)] = max(abs(profMap(:,j)));
        peakX(j) = rez.xc(peakInd(j));
        peakY(j) = rez.yc(peakInd(j));
    end
    distData.Peaks = [peakX peakY];
    distData.PeakInds = peakInd;
    distData.PeakChans = rez.ops.chanMap(peakInd);
    
    % determine centroid distances
    peakDists = sqrt((bsxfun(@minus,peakX,peakX')).^2 + (bsxfun(@minus,peakY,peakY')).^2);
    distData.PeakDists = peakDists;
    
    % get waveform similarity across all electrodes
    waveAllDists = corr(reshape(meanWaves,[numElec*numSamp numClus]));
    distData.WaveAllDists = 1-waveAllDists;
    
    % get spectral similarity across all electrodes
    waveAllSpecDists = corr(reshape(abs(fft(meanWaves,[],2)),[numElec*numSamp numClus]));
    distData.WaveAllSpecDists = 1-waveAllSpecDists;
    
    % get peak waveform similarity
    clear peakWaves;
    for j = 1:numClus
        peakWaves(:,j) = meanWaves(:,peakInd(j),j);
    end
    distData.PeakWaves = peakWaves;
    
    wavePeakDists = corr(peakWaves);
    distData.WavePeakDists = 1-wavePeakDists;
    
    wavePeakSpecDists = corr(abs(fft(peakWaves,[],1)));
    distData.WavePeakSpecDists = 1-wavePeakSpecDists;
end


function UpdateTable(hObject, handles)
    sasData = handles.rez.SortaSort;
    selCols = {handles.pmCol1.String{handles.pmCol1.Value}; ...
               handles.pmCol2.String{handles.pmCol2.Value}; ...
               handles.pmCol3.String{handles.pmCol3.Value}};
    
    % to use when colorized plots are added, see https://www.mathworks.com/matlabcentral/answers/25038-how-to-change-each-cell-color-in-a-uitable     
    colorgen = @(color,text)['<table border=0 width=400 bgcolor=',color,'><TR><TD>',text,'</TD></TR> </table>'];
    % determine if a particular cluster is selected, otherwise orient to
    % first cluster
    if isempty(handles.SelClusters)
        selClusters = 1;
    else
        selClusters = handles.SelClusters;
    end
    
    sortCols = cell(4,1);
    sortCols{1} = arrayfun(@(x)char(x),sasData.ClusterType,'uniformoutput',false);
    colNames = {'AssignClusterType'};
    rowNames = sasData.ClusterID;
    colTypes{1} = {'Unknown' 'Noise' 'MU' 'SU'}; 
    colorYes = false;
    for j = 1:length(selCols)
        currCol = selCols{j};
        switch currCol
            case 'Cluster ID'
                sortCols{j+1} = str2num(strjoin(sasData.ClusterID,' '))';
                colNames{j+1} = 'ClusterID';
                colTypes{j+1} = 'char';
                colorYes(j+1) = false;
                colSortDir{j} = 'ascend';
            case 'Cluster Type'
                sortCols{j+1} = arrayfun(@(x)char(x),sasData.ClusterType,'uniformoutput',false);
                colNames{j+1} = 'ClusterType';
                colTypes{j+1} = 'char';
                colorYes(j+1) = false;
                colSortDir{j} = 'descend';
            case 'Spatial Profile'
                sortCols{j+1} = mean(sasData.DistData.ProfileMapDist(selClusters,:),1)';
                colNames{j+1} = 'SpatialProfile';
                colTypes{j+1} = 'numeric';
                colorYes(j+1) = true;
                colSortDir{j} = 'ascend';
            case 'Peak Location'
                sortCols{j+1} = mean(sasData.DistData.PeakDists(selClusters,:),1)';
                colNames{j+1} = 'PeakLocation';    
                colTypes{j+1} = 'numeric';
                colorYes(j+1) = true;
                colSortDir{j} = 'ascend';
            case 'Complete Waveform'
                sortCols{j+1} = mean(sasData.DistData.WaveAllDists(selClusters,:),1)';
                colNames{j+1} = 'CompleteWaveform';
                colTypes{j+1} = 'numeric';
                colorYes(j+1) = true;
                colSortDir{j} = 'ascend';
            case 'Complete Spectrum'
                sortCols{j+1} = mean(sasData.DistData.WaveAllSpecDists(selClusters,:),1)';
                colNames{j+1} = 'CompleteSpectrum';   
                colTypes{j+1} = 'numeric';
                colorYes(j+1) = true;
                colSortDir{j} = 'ascend';
            case 'Peak Waveform'
                sortCols{j+1} = mean(sasData.DistData.WavePeakDists(selClusters,:),1)';
                colNames{j+1} = 'PeakWaveform';
                colTypes{j+1} = 'numeric';
                colorYes(j+1) = true;
                colSortDir{j} = 'ascend';
            case 'Peak Spectrum'
                sortCols{j+1} = mean(sasData.DistData.WavePeakSpecDists(selClusters,:),1)';
                colNames{j+1} = 'PeakSpectrum';
                colTypes{j+1} = 'numeric';
                colorYes(j+1) = true;
                colSortDir{j} = 'ascend';
            case 'Peak Channel'
                sortCols{j+1} = sasData.DistData.PeakChans';
                colNames{j+1} = 'PeakChannel';
                colTypes{j+1} = 'numeric';
                colorYes(j+1) = true;
                colSortDir{j} = 'ascend';
        end      
    end
    
    distTable = table(sortCols{1},sortCols{2},sortCols{3},sortCols{4});
    [distTable, sortInds] = sortrows(distTable,[2 3 4],colSortDir);
    handles.tblSort.ColumnFormat = colTypes;
    handles.tblSort.Data = table2cell(distTable);
    handles.tblSort.ColumnEditable = [true false false false];
    handles.tblSort.ColumnName = colNames;
    handles.tblSort.RowName = rowNames(sortInds);
    handles.SortedInds = sortInds;
    % colorize numeric entries
    
    %handles.tblSort.RearrangeableColumns = false;
    guidata(hObject, handles);
end


% --- Executes when selected cell(s) is changed in tblSort.
function tblSort_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to tblSort (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
    handles.TblSortInds = eventdata.Indices;
    if ~isempty(eventdata.Indices) && ~any(eventdata.Indices(:,2)==1)
        selClusters = handles.tblSort.RowName(eventdata.Indices(:,1));
        for j = 1:length(selClusters)
            clusterInds(j) = find(strcmp(selClusters(j),handles.rez.SortaSort.ClusterID));
        end
        handles.SelClusters = clusterInds;
        guidata(hObject, handles); 
    end
    UpdatePlots(hObject, handles);
end


function UpdatePlots(hObject, handles)
    numPlots = 6;
    if ~strcmp(hObject.Tag,'pbReplotGraphs')
        if isempty(handles.SelClusters)
            return;
        elseif isequal(handles.SelClusters, handles.CurrPlotClusters)
            return;
        end
    end
    waitH = waitbar(0,'Starting plotting');
    if handles.cbProfile.Value
        if ~isfield(handles,'ProfileH')
            handles.ProfileH = figure('UserData',handles);
        elseif ~ishghandle(handles.ProfileH)
            handles.ProfileH = figure('UserData',handles);
        end
        figure(handles.ProfileH);
        
        clf('reset'); % FIX
        set(gcf,'name','Spatial Profile');
        waitbar(1/numPlots,waitH,'Plotting profile');
        PlotProfile(handles.SelClusters,handles.rez)
    end
    if handles.cbWaveform.Value
        if ~isfield(handles,'WaveformsH')
            handles.WaveformsH = figure('UserData',handles);
            firstTime = true;
        elseif ~ishghandle(handles.WaveformsH)
            handles.WaveformsH = figure('UserData',handles);
            firstTime = true;
        else
            firstTime = false;
            
        end
        figure(handles.WaveformsH);
        origXLim = xlim();
        origYLim = ylim();
        clf('reset'); % FIX
        set(gcf,'name','Waveforms');
        waitbar(2/numPlots,waitH,'Plotting waveforms');
        PlotWaveforms(handles.SelClusters,handles.rez, handles.WaveformsH)
        if ~firstTime
            xlim(origXLim);
            ylim(origYLim);
        end
    end
    if handles.cbTrace.Value
        if ~isfield(handles,'TraceH')
            handles.TraceH = figure('name','Traces','UserData',handles);
        elseif ~ishghandle(handles.TraceH)
            handles.TraceH = figure('name','Traces','UserData',handles);
        end
        figure(handles.TraceH);
        waitbar(3/numPlots,waitH,'Plotting traces');
        PlotTraces(handles.SelClusters,handles.rez, handles.TraceH, ...
            handles.mmfFil, handles.mmfDat)
    end
    if handles.cbXCorr.Value
        if ~isfield(handles,'XCorrH')
            handles.XCorrH = figure('UserData',handles);
        elseif ~ishghandle(handles.XCorrH)
            handles.XCorrH = figure('UserData',handles);
        end
        figure(handles.XCorrH);
        clf('reset'); % FIX
        set(gcf,'name','Cross-correlation');
        waitbar(4/numPlots,waitH,'Plotting cross-correlations');
        PlotXCorr(handles.SelClusters,handles.rez, handles.XCorrH, handles.mmfFil)
    end
    if handles.cbAmplitudes.Value
        if ~isfield(handles,'AmpsH')
            handles.AmpsH = figure('UserData',handles);
        elseif ~ishghandle(handles.AmpsH)
            handles.AmpsH = figure('UserData',handles);
        end
        figure(handles.AmpsH);
        clf('reset'); % FIX
        set(gcf,'name','Amplitudes');
        waitbar(5/numPlots,waitH,'Plotting amplitudes');
        PlotAmplitudes(handles.SelClusters,handles.rez, handles.AmpsH, ...
            handles.mmfFil, handles.mmfDat)
    end
    if handles.cbFeatures.Value
        if ~isfield(handles,'FetsH')
            handles.FetsH = figure('UserData',handles);
        elseif ~ishghandle(handles.FetsH)
            handles.FetsH = figure('UserData',handles);
        end
        figure(handles.FetsH);
        clf('reset'); % FIX
        set(gcf,'name','Features');
        waitbar(6/numPlots,waitH,'Plotting features');
        PlotFeatures(handles.SelClusters,handles.rez, handles.FetsH, ...
            handles.mmfFil, handles.mmfDat)
    end
    figure(handles.figure1);
    close(waitH);
    handles.CurrPlotClusters = handles.SelClusters;
    guidata(hObject, handles);
end

function PlotTraces(selClusters, rez, figH, mmfFil, mmfDat)
    colorset = colormap('lines');
    chanMap = rez.ops.chanMap;
    peakInds = rez.SortaSort.DistData.PeakInds;
    peakChans = chanMap(peakInds(selClusters));
    numTPts = size(mmfFil.Data.file,2);
    uVScale = rez.SortaSort.Settings.uVScale;
    uVSep = rez.SortaSort.Settings.uVSep;
    tWindow = rez.SortaSort.Settings.TraceWindow;
    sampRate = rez.ops.fs;
    tWindowInds = round(tWindow * sampRate);
    for j = 1:length(selClusters)
        currCluster = selClusters(j);
        allSpikes = rez.SortaSort.GroupedSpikes{currCluster};
        numSpks = length(allSpikes);
        subInds = randperm(numSpks, min([10 numSpks]));
        subSpks = allSpikes(subInds);
        subSpks((subSpks<tWindowInds)|(subSpks>(numTPts-tWindowInds))) = [];
        if isempty(subSpks)
            figure(figH);
            subplot(1,length(selClusters),j);
            plot(0);
            ylim([-uVSep uVSep*10]);
            title(num2str(rez.SortaSort.ClusterID{currCluster}));
            continue;
        end
        
        filInds = (subSpks + [-tWindowInds:tWindowInds])';
        datInds = (subSpks + [-tWindowInds:tWindowInds])';
        
        filData = double(mmfFil.data.file(peakChans(j),filInds(:)))*uVScale;
        filData = reshape(filData,[(2*tWindowInds)+1 length(subSpks)]);
        datData = double(mmfDat.data.file(peakChans(j),datInds(:)))*uVScale;
        datData = reshape(datData,[(2*tWindowInds)+1 length(subSpks)]);
        filTData = (-tWindowInds:tWindowInds)/rez.ops.fs;
        datTData = (-tWindowInds:tWindowInds)/rez.ops.fs;
        
        figure(figH);
        subplot(1,length(selClusters),j);
        PlotStackedLines(datData', datTData, 0.5*[1 1 1], 'ABSSEP', uVSep);
        hold on;
        PlotStackedLines(filData',filTData,colorset(mod(j,64)+1,:), 'ABSSEP', uVSep);
        set(gca,'ytick',[]);
        title(num2str(rez.SortaSort.ClusterID{currCluster}));
        
    end
    h = zoom;
    h.Motion = 'horizontal';
    h.Enable = 'on';
end

function PlotAmplitudes(selClusters, rez, figH, mmfFil, mmfDat)
    colorset = colormap('lines');
    % plot firing rates across time
    for j = 1:length(selClusters)
        currClust = selClusters(j);
        if isempty(rez.SortaSort.GroupedSpikes{currClust})
            continue;
        end
        spkTimes = duration(0,0,rez.SortaSort.GroupedSpikes{currClust}/rez.ops.fs);
        
        [fRate, binTimes] = histcounts(spkTimes,(min(spkTimes)-duration(0,0,1)):duration(0,0,1):max(spkTimes));
        
        figure(figH);
        set(figH,'defaultAxesColorOrder',[0 0 0; 0 0 0]);
        yyaxis left;
        areaH(j)= area(binTimes(2:end)-duration(0,0,0.5),fRate, ...
            'FaceColor',colorset(mod(j,64)+1,:), 'linestyle', 'none', ...
            'FaceAlpha', 0.5);
        ylabel('Firing rate (Hz)');
        xtickformat('hh:mm:ss')
        hold on;
    end
    % plot scatter of time points and amplitudes
    for j = 1:length(selClusters)
        currClust = selClusters(j);
        if isempty(rez.SortaSort.GroupedSpikes{currClust})
            continue;
        end
        spkTimes = duration(0,0,rez.SortaSort.GroupedSpikes{currClust}/rez.ops.fs);
        spkAmps = rez.SortaSort.GroupedAmplitudes{currClust};
                
        figure(figH);
        set(figH,'defaultAxesColorOrder',[0 0 0; 0 0 0]);
        yyaxis right;
        userData = {rez selClusters mmfFil mmfDat};
        scatter(spkTimes,spkAmps,'.','MarkerEdgeColor',colorset(mod(j,64)+1,:), ...
            'userdata',userData,'ButtonDownFcn',@Amplitudes_ButtonDownFunc);
        ylabel('Spike Amplitudes');
        
        hold on;
    end
    legend(rez.SortaSort.ClusterID(selClusters));
    set(gca,'tickdir','out');
    hold off;
    
    function Amplitudes_ButtonDownFunc(source,eventdata)
        colorset = colormap('lines');
        rez = source.UserData{1};
        selClusters = source.UserData{2};
        mmfFil = source.UserData{3};
        mmfDat = source.UserData{4};
        timePt = datevec(eventdata.IntersectionPoint(1));
        timePt = (timePt(4)*3600)+(timePt(5)*60)+timePt(6);
        timePt = ceil(timePt*30000);
        chanMap = rez.ops.chanMap;
    	peakInds = rez.SortaSort.DistData.PeakInds;
        peakChans = chanMap(peakInds(selClusters));
        
        
        numTPts = size(mmfFil.Data.file,2);
        
        filInds = timePt + [-3000:3000]';
        filInds(filInds<1) = 1;
        filInds(filInds>numTPts) = numTPts;
        
        datInds = timePt + [-15000:15000]';
        datInds(datInds<1) = 1;
        datInds(datInds>numTPts) = numTPts;
        filData = double(mmfFil.data.file(peakChans,filInds));
        datData = double(mmfDat.data.file(peakChans,datInds));
        
        filColors = colorset(mod(find(1:length(selClusters)),64)+1,:);
        filTData = (-3000:3000)/rez.ops.fs;    
        datTData = linspace(-0.1,0.1,30001);
        figure;
        axes;
        ylim([-1500 ((length(selClusters)-1)*1200)+1500]);
        hold on;
        line([0 0],ylim,'color',[0.5 0.5 0.5])
        PlotStackedLines(filData, filTData, filColors, 'ABSSEP', 1200);
        hold on;
        PlotStackedLines(datData/4, datTData, [0.5 0.5 0.5], 'ABSSEP', 1200);
        
        title('Traces around time point');
        h = zoom;
        h.Motion = 'horizontal';
        h.Enable = 'on';
    end
end
% 
% function PlotFeatures(selClusters, rez, figH, mmfFil, mmfDat)
%     colorset = colormap('lines');
%     
%     % plot scatter of time points and amplitudes
%     for j = 1:length(selClusters)
%         currClust = selClusters(j);
%         if isempty(rez.SortaSort.GroupedPCs{currClust})
%             continue;
%         end
%         spkFets = rez.SortaSort.GroupedPCs{currClust};
%                 
%         figure(figH);
%         set(figH,'defaultAxesColorOrder',[0 0 0; 0 0 0]);
%         
%         userData = {rez selClusters mmfFil mmfDat};
%         scatter3(spkFets(:,1),spkFets(:,2),spkFets(:,3),'.',...
%             'MarkerEdgeColor',colorset(mod(j,64)+1,:), ...
%             'userdata',userData,'ButtonDownFcn',@Features_ButtonDownFunc);
%         hold on;
%     end
%     xlabel('PC 1');
%     ylabel('PC 2');
%     zlabel('PC 3');
%     legend(rez.SortaSort.ClusterID(selClusters));
%     set(gca,'tickdir','out');
%     axis tight;
%     hold off;
%     % ADD GROUPED PCS TO MERGE
%     % HANDLE SPLITTING
%     function Features_ButtonDownFunc(source,eventdata)
%         colorset = colormap('lines');
%         rez = source.UserData{1};
%         selClusters = source.UserData{2};
%         mmfFil = source.UserData{3};
%         mmfDat = source.UserData{4};
%         pcPt = eventdata.IntersectionPoint;
%         spkPCs = vertcat(rez.SortaSort.GroupedPCs{selClusters});
%         spkInd = find((spkPCs(:,1)==pcPt(1))&(spkPCs(:,2)==pcPt(2))&...
%             (spkPCs(:,3)==pcPt(3)),1);
%         spkTimes = vertcat(rez.SortaSort.GroupedSpikes{selClusters});
%         timePt = spkTimes(spkInd);
%         chanMap = rez.ops.chanMap;
%     	peakInds = rez.SortaSort.DistData.PeakInds;
%         peakChans = chanMap(peakInds(selClusters));
%         
%         
%         numTPts = size(mmfFil.Data.file,2);
%         
%         filInds = timePt + [-3000:3000]';
%         filInds(filInds<1) = 1;
%         filInds(filInds>numTPts) = numTPts;
%         
%         datInds = timePt + [-15000:15000]';
%         datInds(datInds<1) = 1;
%         datInds(datInds>numTPts) = numTPts;
%         filData = double(mmfFil.data.file(peakChans,filInds));
%         datData = double(mmfDat.data.file(peakChans,datInds));
%         
%         filColors = colorset(mod(find(1:length(selClusters)),64)+1,:);
%         filTData = (-3000:3000)/rez.ops.fs;    
%         datTData = linspace(-0.1,0.1,30001);
%         figure;
%         axes;
%         ylim([-1500 ((length(selClusters)-1)*1200)+1500]);
%         hold on;
%         line([0 0],ylim,'color',[0.5 0.5 0.5])
%         PlotStackedLines(filData, filTData, filColors, 'ABSSEP', 1200);
%         hold on;
%         PlotStackedLines(datData/4, datTData, [0.5 0.5 0.5], 'ABSSEP', 1200);
%         
%         title('Traces around time point');
%         h = zoom;
%         h.Motion = 'horizontal';
%         h.Enable = 'on';
%     end
% end


function PlotWaveforms(selClusters, rez, figH)
    numSamps = size(rez.SortaSort.DistData.MeanSpkWaves,1);
    uVScale = rez.SortaSort.Settings.uVScale;
    sampRate = rez.ops.fs;
    uVScaleBar = rez.SortaSort.Settings.uVScaleBar;
    msScaleBar = rez.SortaSort.Settings.msScaleBar;
    colorset = colormap('lines');
    xCoords = rez.xc;
    yCoords = rez.yc;
    diffX = abs(bsxfun(@minus, xCoords, xCoords'));
    diffY = abs(bsxfun(@minus, yCoords, yCoords'));
    medDiffX = min(diffX(diffX~=0));
    medDiffY = min(diffY(diffY~=0));
    xScale = medDiffX/(numSamps/sampRate);
    
    traces = rez.SortaSort.DistData.MeanSpkWaves(:,:,selClusters)*uVScale;
    traceRanges = max(traces,2)-min(traces,2);
    traceRange = max(traceRanges(:));
    yScale = medDiffY/traceRange;
    
    xVals = cumsum(ones(size(traces)),1)/sampRate;
    % NOTE, HAD xCoords and yCoords transposed originally, should repair this
    xVals = xVals*xScale+xCoords(:)';
    yVals = traces*yScale+yCoords(:)';
    
    for j = 1:length(selClusters)
        currColor = colorset(mod(j,64)+1,:);
        figure(figH);
        plot(xVals(:,:,j),yVals(:,:,j),'color',currColor);
        hold on;
    end
    
    % plot scale bar for time
    line((max(xVals(:))+(medDiffX)/2)+[0 msScaleBar*xScale*0.001], ...
         (min(yVals(:))-(medDiffY))*[1 1],'color','black','LineWidth',2);
    % plot scale bar for voltage
    line((max(xVals(:))+(medDiffX)/2)*[1 1], ...
         (min(yVals(:))+(medDiffY))+[0 uVScaleBar*yScale],'color','black','LineWidth',2);
    
    hold off;
end
    


% --- Executes when entered data in editable cell(s) in tblSort.
function tblSort_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to tblSort (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
    clusterTypes = hObject.Data(:,1);
    clusterIDs = hObject.RowName;
    clusterInds = zeros(size(clusterTypes));
    for j = 1:length(clusterIDs)
        currInd = find(strcmp(clusterIDs{j},handles.rez.SortaSort.ClusterID));
        handles.rez.SortaSort.ClusterType(currInd) = categorical(clusterTypes(j),...
            {'Unknown' 'Noise' 'MU' 'SU'});
    end
    
    % recenter the scroll bar on selected element. Adapted from undocumented
    % matlab
    [jScrollpane, levels, parentIdx, listing] = findjobj(hObject);
    scrollValue = jScrollpane.getVerticalScrollBar.getValue;    
    UpdateTable(hObject,handles);
    drawnow;
    jScrollpane.getVerticalScrollBar.setValue(scrollValue);
    guidata(hObject,handles);
end   

function PlotXCorr(selClusters, rez, figH, mmfFil)
    colorset = colormap('lines');
    numClusters = length(selClusters);
    chanMap = rez.ops.chanMap;
    sampRate = rez.ops.fs;
    peakInds = rez.SortaSort.DistData.PeakInds;
    xCorrWindow = rez.SortaSort.Settings.XCorrWindow;
    dt = round(sampRate*0.001);
    xCorrWindow = round((xCorrWindow*sampRate)/dt);
    
%     bins = -3000:dt:3000;
    binCenters = -xCorrWindow:xCorrWindow;
    for j = 1:numClusters
        trigCluster = selClusters(j);
        trigSpikes = rez.SortaSort.GroupedSpikes{trigCluster};
        
        if length(trigSpikes)>10000
            missRatio = length(trigSpikes)/10000;
            trigSpikes = sort(trigSpikes(randperm(length(trigSpikes),10000)));
        else
            missRatio = 1;
        end
        for k = j:numClusters
            
            otherCluster = selClusters(k);
            if k == j
                otherSpikes = trigSpikes;
                plotColor = colorset(mod(j,64)+1,:);
            else
                otherSpikes = rez.SortaSort.GroupedSpikes{otherCluster};
                plotColor = [0 0 0];
            end
            
            subplot(numClusters,numClusters,j+((k-1)*numClusters));
            
            spkHist = zeros(1,length(binCenters));
            trigSpkBins = cell(1,length(binCenters));
            for p = 1:length(trigSpikes)
                currBins = round((otherSpikes-trigSpikes(p))/dt);
                if isempty(currBins)
                    continue;
                end
                ltEdge = find(currBins >= -xCorrWindow,1,'first');
                if ltEdge > 1
                    otherSpikes(1:(ltEdge-1)) = [];
                end
                gtEdge = find(currBins > xCorrWindow,1,'first')-1;
                
                currBins = currBins(ltEdge:gtEdge)+xCorrWindow+1;
                if isempty(currBins)
                    continue;    
                end
                spkHist = spkHist + accumarray(currBins,1,[(2*xCorrWindow)+1 1])';
%                 [currSpkHist, ~, currBins] = histcounts(otherSpikes-trigSpikes(p),bins);
%                 currBins(currBins==0)=[];
%                 spkHist = spkHist + accumarray(ceil(currSpkHist/dt)+;
                for m = 1:length(currBins)
                    trigSpkBins{currBins(m)} = [trigSpkBins{currBins(m)} trigSpikes(p)];
                end
            end
            spkHist = missRatio*(spkHist/length(trigSpikes))/0.001;
            if k == j
                spkHist(xCorrWindow+1) = 0;
                fRate = sampRate/mean(diff(rez.SortaSort.GroupedSpikes{trigCluster}));
%                     length(rez.SortaSort.GroupedSpikes{trigCluster})/...
%                     ((max(rez.SortaSort.GroupedSpikes{trigCluster})-...
%                     min(rez.SortaSort.GroupedSpikes{trigCluster}))/sampRate);
                if isempty(fRate)
                    fRate = 0;
                end
                arrowStr = '<->';
                frString = [' ' num2str(fRate,3) 'Hz'];
            else
                arrowStr = '->';
                frString = '';
            end
            
            % create data to save for generting trace plots based on spike
            % ISIs
            isiData = {rez trigCluster otherCluster trigSpkBins mmfFil selClusters};
            
            figure(figH);
            
            barH = bar(binCenters,spkHist,1, ...
                'EdgeColor','none','FaceColor',plotColor, ...
                'UserData', isiData, 'ButtonDownFcn',@XCorr_ButtonDownFunc);
            hold on;
            if j == k
                line(xlim, fRate*[1 1], 'color', 'black', 'LineStyle', '--');
            end
            hold off;
            title([num2str(rez.SortaSort.ClusterID{trigCluster}) arrowStr ...
                num2str(rez.SortaSort.ClusterID{otherCluster}) frString]);
            set(gca,'xtick',[0]);
            set(gca,'ytick',[]);
            set(gca,'tickdir','out');
            xlim([min(binCenters) max(binCenters)]);
            drawnow;
        end
    end
    
    function XCorr_ButtonDownFunc(source,eventdata)
        colorset = colormap('lines');
        numBins = length(source.UserData{4});
        rez = source.UserData{1};
        trigCluster = source.UserData{2};
        otherCluster = source.UserData{3};
        mmfFil = source.UserData{5};
        selClusters = source.UserData{6};
        uVScale = rez.SortaSort.Settings.uVScale;
        uVSep = rez.SortaSort.Settings.uVSep;
        xCorrWindow = rez.SortaSort.Settings.XCorrWindow;
        xCorrWindow = round(xCorrWindow*rez.ops.fs);
        chanMap = rez.ops.chanMap;
    	peakInds = rez.SortaSort.DistData.PeakInds;
        trigChan = chanMap(peakInds(trigCluster));
        otherChan = chanMap(peakInds(otherCluster));
        
        if trigCluster == otherCluster
            arrowStr = '<->';
        else
            arrowStr = '->';
        end
        
        centerBin = ceil(numBins/2);
        selBin = round(eventdata.IntersectionPoint(1))+centerBin;
        
        if selBin > numBins
            selBin = numBins;
        end
        subSpks = source.UserData{4}{selBin}';
        if length(subSpks)>10
            subSpks = subSpks(randperm(length(subSpks),10));
        end
        numTPts = size(mmfFil.Data.file,2);
        if isempty(subSpks)
            return;
        end
        filInds = (subSpks + [-xCorrWindow:xCorrWindow])';
        filInds(filInds<1) = 1;
        filInds(filInds>numTPts) = numTPts;
        trigData = double(mmfFil.data.file(trigChan,filInds(:)))*uVScale;
        otherData = double(mmfFil.data.file(otherChan,filInds(:)))*uVScale;
        trigData = reshape(trigData,[(2*xCorrWindow)+1 length(subSpks)]);
        otherData = reshape(otherData,[(2*xCorrWindow)+1 length(subSpks)]);
        trigColor = colorset(mod(find(trigCluster==selClusters),64)+1,:);
        otherColor = colorset(mod(find(otherCluster==selClusters),64)+1,:);
        filTData = (-xCorrWindow:xCorrWindow)/rez.ops.fs;        
        figure;
        axes;
        ylim([-uVSep ((length(subSpks)-1)*uVSep)+uVSep]);
        hold on;
        line([0 0],ylim,'color',[0.5 0.5 0.5])
        line([1 1]*((selBin-centerBin)/1000),ylim,'color',[0.5 0.5 0.5])
        PlotStackedLines(trigData', filTData, trigColor, 'ABSSEP', uVSep);
        
        if trigChan ~= otherChan
            hold on;
            PlotStackedLines(otherData', filTData, otherColor, 'ABSSEP', uVSep);
        end
        
        title([num2str(rez.SortaSort.ClusterID{trigCluster}) arrowStr ...
                num2str(rez.SortaSort.ClusterID{otherCluster}) ...
                ' lagged by ' num2str(selBin-centerBin) ' ms']);
        h = zoom;
        h.Motion = 'horizontal';
        h.Enable = 'on';
    end
end



        


% --- Executes on button press in pbReplotGraphs.
function pbReplotGraphs_Callback(hObject, eventdata, handles)
% hObject    handle to pbReplotGraphs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    UpdatePlots(hObject, handles);
end


% --- Executes on button press in pbSplit.
function pbSplit_Callback(hObject, eventdata, handles)
% hObject    handle to pbSplit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    selClusters = handles.SelClusters;
    rez = handles.rez;
    if length(selClusters)~=1
        errordlg('Only one cluster can be split at a time');
        return;
    end

    if ~isfield(handles,'AmpsH')
        handles.AmpsH = figure();
    elseif ~ishghandle(handles.AmpsH)
        handles.AmpsH = figure();
    end
    figure(handles.AmpsH);
    clf('reset'); % FIX
    set(gcf,'name','Amplitudes');
    PlotAmplitudes(selClusters,rez, handles.AmpsH, ...
            handles.mmfFil, handles.mmfDat);
    yyaxis right;
    polyH = impoly();
    if isempty(polyH)
        return;
    end
    polyPos = polyH.getPosition;
    polyH.delete;
    origSpks = rez.SortaSort.GroupedSpikes{selClusters};
%     origPCs = rez.SortaSort.GroupedPCs{selClusters};
    spkTimes = datenum(duration(0,0,origSpks/rez.ops.fs));
    spkAmps = rez.SortaSort.GroupedAmplitudes{selClusters};

    inSpks = inpolygon(spkTimes,spkAmps,polyPos(:,1),polyPos(:,2));

    
    rez.SortaSort.ClusterID{end+1} = num2str(max(cellfun(@str2num,rez.SortaSort.ClusterID))+1);
    rez.SortaSort.ClusterType([selClusters end+1]) = categorical({'Unknown' 'Unknown'},...
    {'Unknown' 'Noise' 'MU' 'SU'});

    rez.SortaSort.GroupedSpikes{selClusters} = origSpks(inSpks);
    rez.SortaSort.GroupedSpikes{end+1} = origSpks(~inSpks);

    rez.SortaSort.GroupedAmplitudes{selClusters} = spkAmps(inSpks);
    rez.SortaSort.GroupedAmplitudes{end+1} = spkAmps(~inSpks);

%     rez.SortaSort.GroupedPCs{selClusters} = origPCs(inSpks,:);
%     rez.SortaSort.GroupedPCs{end+1} = origPCs(~inSpks,:);
    % retrieve a new set of spike waveforms using the new clusters
    handles.rez = rez;
    for j = [selClusters length(rez.SortaSort.GroupedSpikes)]
        rez.SortaSort.SpkWaves{j} = ExtractSpkWaveforms(j, 100, handles);
    end
    rez.SortaSort.DistData = CalculateDistances(rez);
    handles.SelClusters = [selClusters numel(rez.SortaSort.ClusterID)];

    handles.rez = rez;

    guidata(hObject, handles);
    UpdateTable(hObject, handles);
    UpdatePlots(hObject, handles);
end

% 
% % --- Executes on key press with focus on any of the data plots
% function Plots_WindowKeyPressFcn(source, eventdata)
% % source    handle to current figure (see GCBO)
% % eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
% %	Key: name of the key that was pressed, in lower case
% %	Character: character interpretation of the key(s) that was pressed
% %	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% 	if isempty(handles.TblSortInds)
%         return;
%     end
%     currRow = handles.TblSortInds(:,1);
%     currCol = handles.TblSortInds(:,2);
%    
%     % taken from: https://www.mathworks.com/matlabcentral/answers/348905-programatically-selecting-cells-in-a-uitable
%     jUIScrollPane = findjobj(handles.tblSort);
%     jUITable = jUIScrollPane.getViewport.getView;
%     jUITable.changeSelection(row-1,col-1, false, false);
% %     switch eventdata.Key
%         
% end


% --- Executes on button press in cbFeatures.
function cbFeatures_Callback(hObject, ~, handles)
% hObject    handle to cbFeatures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbFeatures
end
% 
% % --- Executes on button press in pbSplitFet.
% function pbSplitFet_Callback(hObject, eventdata, handles)
% % hObject    handle to pbSplitFet (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% colorset = colormap('lines');
%  selClusters = handles.SelClusters;
%     rez = handles.rez;
%     if length(selClusters)~=1
%         errordlg('Only one cluster can be split at a time');
%         return;
%     end
% 
%     if ~isfield(handles,'FetsH') || ~ishghandle(handles.FetsH)
%         handles.FetsH = figure();
%         set(gcf,'name','Features');
%         PlotFeatures(handles.SelClusters, handles.rez, handles.FetsH, ...
%             handles.mmfFil, handles.mmfDat)
%     end
%     figure(handles.FetsH);
%     set(gcf,'defaultAxesColorOrder',[0 0 0; 0 0 0]);
% 
%     viewAng = get(gca,'View');
%     projTrans = viewmtx(viewAng(1),viewAng(2));
%     origPts = rez.SortaSort.GroupedPCs{handles.SelClusters};
%     projPts = projTrans*[origPts'; ones(1,length(origPts))];
%     projH = figure;
%     set(projH,'Name','Projected Features');
%     scatter(projPts(1,:),projPts(2,:),'.','MarkerEdgeColor',colorset(mod(1,64)+1,:));
%     
%     polyH = impoly();
%     if isempty(polyH)
%         return;
%     end
%     polyPos = polyH.getPosition;
%     polyH.delete;
%     close(projH);
%     inSpks = inpolygon(projPts(1,:),projPts(2,:),polyPos(:,1),polyPos(:,2));
%     origSpks = rez.SortaSort.GroupedSpikes{selClusters};
%     origPCs = rez.SortaSort.GroupedPCs{selClusters};
%     spkAmps = rez.SortaSort.GroupedAmplitudes{selClusters};
% 
%     rez.SortaSort.dWU = cat(3,rez.SortaSort.dWU, rez.SortaSort.dWU(:,:,selClusters));
%     rez.SortaSort.ClusterID{end+1} = num2str(max(cellfun(@str2num,rez.SortaSort.ClusterID))+1);
%     rez.SortaSort.ClusterType([selClusters end+1]) = categorical({'Unknown' 'Unknown'},...
%     {'Unknown' 'Noise' 'MU' 'SU'});
% 
%     rez.SortaSort.GroupedSpikes{selClusters} = origSpks(inSpks);
%     rez.SortaSort.GroupedSpikes{end+1} = origSpks(~inSpks);
% 
%     rez.SortaSort.GroupedAmplitudes{selClusters} = spkAmps(inSpks);
%     rez.SortaSort.GroupedAmplitudes{end+1} = spkAmps(~inSpks);
% 
%     rez.SortaSort.GroupedPCs{selClusters} = origPCs(inSpks,:);
%     rez.SortaSort.GroupedPCs{end+1} = origPCs(~inSpks,:);
%     
%     rez.SortaSort.DistData = CalculateDistances(rez);
%     handles.SelClusters = [selClusters numel(rez.SortaSort.ClusterID)];
% 
%     handles.rez = rez;
% 
%     guidata(hObject, handles);
%     UpdateTable(hObject, handles);
%     UpdatePlots(hObject, handles);
% end


% --------------------------------------------------------------------
function menuSettings_Callback(hObject, eventdata, handles)
% hObject    handle to menuSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    settingsStruct = handles.rez.SortaSort.Settings;
    fNames = fieldnames(settingsStruct);
    oldSetVals = struct2cell(settingsStruct);
    oldSetVals = cellfun(@(x)num2str(x),oldSetVals,'UniformOutput',false);
    
    setOut = inputdlg(fNames,'Settings',[1 40],oldSetVals);
    
    for j = 1:length(setOut)
        settingsStruct.(fNames{j}) = str2num(setOut{j});
    end
    handles.rez.SortaSort.Settings = settingsStruct;
    guidata(hObject, handles);

end

function tblSort_KeyPressFcn(hObject, eventdata, handles)
    keyChar = eventdata.Character;
    switch keyChar
        case 'n'
            handles.rez.SortaSort.ClusterType(handles.SelClusters) = {'Noise'};
            guidata(hObject, handles);
            UpdateTable(hObject, handles);
        case 'm'
            handles.rez.SortaSort.ClusterType(handles.SelClusters) = {'MU'};
            guidata(hObject, handles);
            UpdateTable(hObject, handles);
        case 's'
            handles.rez.SortaSort.ClusterType(handles.SelClusters) = {'SU'};
            guidata(hObject, handles);
            UpdateTable(hObject, handles);
         case 'u'
            handles.rez.SortaSort.ClusterType(handles.SelClusters) = {'Unknown'};
            guidata(hObject, handles);
            UpdateTable(hObject, handles);
    end
end

function spkWaves = ExtractSpkWaveforms(clusterInd, numSpks, handles)
    [numChan, recLength] = size(handles.mmfFil.Data.file);
    spkTimes = handles.rez.SortaSort.GroupedSpikes{clusterInd};
    spkTimes(spkTimes<31) = [];
    spkTimes(spkTimes>(recLength-31)) = [];
    if isempty(spkTimes)
        spkWaves = nan(61,numChan);
    else
        subSpkTimes = spkTimes(randperm(length(spkTimes),min([length(spkTimes) numSpks])));
        spkInds = (subSpkTimes+(-30:30))';
        spkWaves = handles.mmfFil.Data.file(:,spkInds(:));
        spkWaves = reshape(spkWaves,numChan,61,length(subSpkTimes));
        spkWaves = permute(spkWaves, [2 1 3]);
        spkWaves = double(spkWaves)-mean(spkWaves,1);
    end
end
    