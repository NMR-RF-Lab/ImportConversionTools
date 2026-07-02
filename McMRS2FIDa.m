% Import data from McMRSGUI and convert to FID-A format

function [newStruct,oldStruct] = McMRS2FIDa(matnm,fname_ID)
% DESCRIPTION:
% Reads in Matlab structure output from McMRSGUI as well as original
% Siemens twix raw data (.dat file) using the io_loadspec_twix.m function.
% It then converts .mat data to a FID-A structure.
% 
% INPUTS:
% matnm     = filename of McMRSGUI .mat data in structure format, 
%             i.e. {'nameoffile.mat'}. If directory is not included in name, folder
%             where data is saved must be open.
% fname_ID  = twix ID (a number) or file name of data (.dat, .fid, or .rda). If
%             inputting a number, folder where data is saved must be the current directory.
%
% OUTPUTS:
% newStruct = Input dataset in FID-A structure format.
% oldStruct = Original dataset in FID-A structure format.

    % Load McMRSGUI mat 
    data = load(matnm, 'Processed');

    % Determine data type
    if isfield(data, 'McMRSData')
        data = data.McMRSData;
        dtype = 'McMRSData';
    elseif isfield(data, 'Starting')
        dtype = 'Starting';
        data = data.Starting;
    elseif isfield(data, 'Processed')
        dtype = 'Processed';
        data = data.Processed;
    end
    
    %  Check if fname_ID is a number of a string
    if ischar(fname_ID)
        % Load FID-A structure and create the new struct based on the old one
        EXT={'.fid','.dat','.rda'}; % the list of extensions wanted
        
        if fname_ID(end-4:end) ~= ('.fid' || '.dat' || '.rda')
            for e = 1:length(EXT)
                if EXT{e} == '.fid'
                    fname_ID = ls([fname_ID '*.fid']);
                else
                    D=dir([fname_ID EXT{e}]); % directory of the path with all extensions
                    fname_ID = D.name;
                end
            end
        end
    
        switch fname_ID(end-4:end)
            case '.dat'
                oldStruct = io_loadspec_twix(fname_ID, 'ZIP', true);
            case '.rda'
                oldStruct = io_loadspec_rda(fname_ID);
            case '.fid'
                oldStruct = io_loadspec_varian(fname_ID);
        end
    else % if not string, then number
        oldStruct = io_loadspec_twix(fname_ID, 'ZIP', true);
    end

    newStruct = oldStruct;

    % Adjust dimensions to match FID-A structures (Npts x Nchann x Navgs)
    % Note: All McMRSGUI matlab structures should only have 1 channel
    switch dtype
        case 'McMRSData'
            newStruct.fids = permute(data.TimeDomain,[3 1 2]);
            newStruct.specs = permute(data.Spectrum,[3 1 2]);
        case 'Starting'
            newStruct.fids = permute(data.Data,[3 1 2]);
            newStruct.specs = fftshift(ifft(newStruct.fids,[],1),1);
        case 'Processed'
            if isfield(data, 'Time_domain')
                newStruct.fids = permute(data.Time_domain,[2 1]);
            else
                newStruct.fids = permute(data.TimeDomain,[2 1]);
            end
            newStruct.specs = permute(data.Spectrum,[2 1]);
    end
    newStruct.sz = size(newStruct.fids);

    % Spectra must be flipped from left to right to match FID-A structure
    data.Spectrum = fliplr(data.Spectrum);

    % Determine dimensionality
    switch dtype
        case 'McMRSData'
            if data.Size.Channels == 1 % 1 coil or coils combined
                newStruct.flags.addedrcvrs = 1;
                newStruct.dims.coils = 0;
            end
            if data.Size.Averages == 1 % 1 average or averages combined
                newStruct.flags.averaged = 1;
                newStruct.dims.averages = 0;
                newStruct.averages = 1;
            end
        case 'Starting'
            if size(data,1) == 1 % 1 coil or coils combined
                newStruct.flags.addedrcvrs = 1;
                newStruct.dims.coils = 0;
            end
            if size(data,2) == 1 % 1 average or averages combined
                newStruct.flags.averaged = 1;
                newStruct.dims.averages = 0;
                newStruct.averages = 1;
            end
        case 'Processed' % data is always only 1 channel 
            newStruct.flags.addedrcvrs = 1;
            newStruct.dims.coils = 0;
            if size(data,1) == 1 % 1 average or averages combined
                newStruct.flags.averaged = 1;
                newStruct.dims.averages = 0;
                newStruct.averages = 1;
            else
                newStruct.dims.averages = 2;
            end
    end

    % Add in phase info
    switch dtype
        case 'McMRSData'
            if ~isempty(data.Operations.ZeroOrderPhase)
                newStruct.flags.phasecorrected = 1;
                newStruct.results.phasecorrected.phShift1 = data.Operations.FirstOrderPhase;
                if isscalar(data.Operations.ZeroOrderPhase)
                    newStruct.results.phasecorrected.phShift = data.Operations.ZeroOrderPhase;
                else
                    newStruct.flags.freqcorrected = 1;
                    newStruct.results.freqcorrected.phsCum = ones(newStruct.rawaverages,data.Size.Channels).*-data.Operations.ZeroOrderPhase;
                    if size(data.Operations.ZeroOrderPhase,2) > 1 || (data.Size.Coils == 1 && data.Channels > 1)
                        newStruct.results.freqcorrected.SeparateRCVRs = true;
                    else
                        newStruct.results.freqcorrected.SeparateRCVRs = false;
                    end
                    newStruct = AddFreqCorrFields(newStruct);
                end
            end
        case 'Processed'
            if isfield(data, 'a0') && ~isempty(data.a0) && data.a0~=0
                newStruct.flags.phasecorrected = 1;
                newStruct.results.phasecorrected.phShift1 = data.a1(1);
                if isscalar(data.Operations.ZeroOrderPhase)
                    newStruct.results.phasecorrected.phShift = data.a0;
                else
                    newStruct.flags.freqcorrected = 1;
                    newStruct.results.freqcorrected.phsCum = ones(newStruct.rawaverages,data.Size.Channels).*-data.a0;
                    newStruct = AddFreqCorrFields(newStruct);
                end
                newStruct.results.freqcorrected.SeparateRCVRs = true;
            elseif isfield(data, 'A0') && ~isempty(data.A0) && data.A0~=0
                newStruct.flags.phasecorrected = 1;
                newStruct.results.phasecorrected.phShift1 = data.A1(1);
                if isscalar(data.Operations.ZeroOrderPhase)
                    newStruct.results.phasecorrected.phShift = data.A0;
                else
                    newStruct.flags.freqcorrected = 1;
                    newStruct.results.freqcorrected.phsCum = ones(newStruct.rawaverages,data.Size.Channels).*-data.A0;
                    newStruct = AddFreqCorrFields(newStruct);
                end
                newStruct.results.freqcorrected.SeparateRCVRs = true;
            end
    end

    % Add in remaining flags
    if isfield(data, 'Operations') && isfield(data.Operations, 'ZeroPadding') && ~isempty(data.Operations.ZeroPadding)
        newStruct.flags.zeropadded = 1;
        newStruct.results.filtered.zpFactor = data.Operations.ZeroPadding + 1;
    end
    if isfield(data, 'Operations') && isfield(data.Operations, 'Apodization') && ~isempty(data.Operations.Apodization)
        newStruct.flags.filtered = 1;
        newStruct.results.filtered.lb = data.Operations.Apodization;
    end
    if isfield(data, 'Operations') && isfield(data.Operations, 'FrequencyShift') && ~isempty(data.Operations.FrequencyShift)
        newStruct.flags.freqcorrected = 1;
        switch dtype
            case 'McMRSData'
                newStruct.results.freqcorrected.fsCum = ones(newStruct.rawaverages,data.Size.Channels).*-data.Operations.FrequencyShift;
                if size(data.Operations.FrequencyShift,2) > 1 || (data.Size.Coils == 1 && data.Channels > 1)
                    newStruct.results.freqcorrected.SeparateRCVRs = true;
                else
                    newStruct.results.freqcorrected.SeparateRCVRs = false;
                end
            case 'Processed'
                newStruct.results.freqcorrected.fsCum = ones(newStruct.rawaverages,1).*-data.Operations.FrequencyShift;
                newStruct.results.freqcorrected.SeparateRCVRs = true;
        end
        newStruct = AddFreqCorrFields(newStruct);
    end

    function struct = AddFreqCorrFields(struct)
        struct.results.freqcorrected.MaxIters = NaN;
        struct.results.freqcorrected.MaxSubIters = NaN;
        struct.results.freqcorrected.fsPolyThresh = NaN;
        struct.results.freqcorrected.phsPolyThresh = NaN;
        struct.results.freqcorrected.SubIter = NaN;
        struct.results.freqcorrected.ndcFLAG = false;
        struct.results.freqcorrected.tmax = NaN;
        struct.results.freqcorrected.Iter = NaN;
    end
end