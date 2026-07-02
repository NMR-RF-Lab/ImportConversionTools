%LOADDAT Import raw data *.dat files & header data. Requires FID-A software package.
%
%   dat = loaddat(fID)
%
%   Inputs
%           MID: number ID for measurement (must be in current folder)
%        KeepOS: *optional* phase oversampling; default = false
%           ZIP: *optional* 2x zero-interpolation filling multiplier; default = false
%       CombAvg: *optional* combine averages; default = false
%          Save: *optional* save spectroscopy data for McMRSGUI, default = false
%
%   Outputs
%           dat: structure containing image & header
% 
% Jacob Degitz, Texas A&M University
% Created 11/7/2023
% Last edited 04/24/2026

function dat = loaddat(MID, varargin)

% Parse inputs
[rmosFLAG, zipFLAG, avgFLAG, svFLAG] = parseInputs(varargin{:});

% Load data using FID-A function
if rmosFLAG
    if avgFLAG
        keys = {'removeos', 'doaverage', 'ignseg'};
    else
        keys = {'removeos', 'ignseg'};
    end
else
    if avgFLAG
        keys = {'doaverage', 'ignseg'};
    else
        keys = {'ignseg'};
    end
end
twix = cell(length(MID),1);
for m = 1:length(MID)
    twix{m} = mapVBVD(MID(m), keys{:});
end

multiMID = true;
if isscalar(MID)
    multiMID = false;
end

% Iterate through MIDs
twix_all = twix;
dat_all = cell(size(twix_all));
for m = 1:length(MID)
    twix = twix_all{m};

    % Check which data type
    switch twix.hdr.Config.IceProgramName
        case {'%SiemensIceProgs%\IceProgram2D', '%SiemensIceProgs%\IceProgram3D', '%SiemensIceProgs%\IceProgramOnline2D', '%SiemensIceProgs%\IceProgramOnline3D', '%SiemensIceProgs%\IceProgRFMap', '%SiemensIceProgs%\IceDixon'}
            % Save header info
            hdr = struct( ...
                   'SoftwareVersion', twix.hdr.Dicom.SoftwareVersions, ...
                           'MeasUID', twix.hdr.Config.MeasUID, ...
                           'Patient', twix.hdr.Meas.PatientName, ...
                   'PatientPosition', twix.hdr.Dicom.tPatientPosition, ...
                                'B0', twix.hdr.Dicom.flMagneticFieldStrength, ... in T
                                'F0', twix.hdr.Dicom.lFrequency, ... in Hz
                           'Nucleus', twix.hdr.Spice.ResonantNucleus, ...
                'ReferenceAmplitude', twix.hdr.Dicom.flTransRefAmpl, ... in V
                          'Sequence', twix.hdr.Config.SequenceString, ... % Dicom.tScanningSequence, Dicom.tSequenceVariant
                      'TransmitCoil', twix.hdr.Dicom.TransmittingCoil, ...
                         'FlipAngle', twix.hdr.Dicom.adFlipAngleDegree, ...
                        'DataMatrix', struct('Size', [size(twix.image(:,:,:,:,:,:,:,:,:,:,:),1) twix.image.NLin twix.image.NPar twix.image.NSli twix.image.NEco twix.image.NAve twix.image.NRep twix.image.NCha], ...
                                   'NumberOfEchoes', twix.image.NEco, ...
                                 'NumberOfAverages', twix.image.NAve, ...
                                    'NumberOfSlabs', twix.image.NSli, ...
                              'NumberOfSliceGroups', twix.image.NPar, ...
                                   'NumberOfSlices', twix.image.NSli, ...
                             'NumberOfMeasurements', twix.image.NRep, ...
                                 'NumberOfChannels', twix.image.NCha), ...
                          'EchoTime', cell2mat(twix.hdr.MeasYaps.alTE(1:twix.image.NEco))./1000, ...
                    'RepetitionTime', twix.hdr.Config.TR/1000, ... in msec
                     'InversionTime', [], ...
                         'DwellTime', twix.hdr.Phoenix.sRXSPEC.alDwellTime{1}, ... in usec???
                         'Geometry', struct('SlabThickness', [], ...
                                            'SliceThickness', [], ...
                                            'DistanceFactor', [], ...
                                            'NormalPlane', "", ...
                                            'Slice', struct('VoxelSize', [], ...
                                                                  'FOV', [], ...
                                                               'Normal', [], ...
                                                               'Offset', [])), ...
                       'ReadoutMode', "", ...
                       'TurboFactor', twix.hdr.Config.TurboFactor, ...
                          'ScanTime', twix.hdr.MeasYaps.lTotalScanTimeSec, ... in sec
                    'RFPulseVoltage', NaN);

            hdr.FFTScale = twix.hdr.Phoenix.asCoilSelectMeas{1,1}.aFFT_SCALE{1,1}.flFactor; % FFT factor?

            % Remove turbo factor if equal to 1
            if hdr.TurboFactor == 1
                hdr = rmfield(hdr, 'TurboFactor');
            end

            % Extract pulse voltage (not present if ref amp = 0)
            if isfield(twix.hdr.Phoenix.sTXSPEC.aRFPULSE{1,1}, "flAmplitude")
                hdr.RFPulseVoltage = twix.hdr.Phoenix.sTXSPEC.aRFPULSE{1,1}.flAmplitude; % in V
            end

            % Extract inversion time (if applicable)
            if twix.hdr.Config.TI ~= 0
                hdr.InversionTime = twix.hdr.Config.TI/1000; % in msec
            else
                hdr = rmfield(hdr, 'InversionTime');
            end


            % Extract normal plane ---------- does not work for TrueFISP or 3D ----------
            if twix.hdr.Spice.VoiNormalCor == 1
                hdr.Geometry.NormalPlane = "Coronal";
            elseif twix.hdr.Spice.VoiNormalSag == 1
                hdr.Geometry.NormalPlane = "Sagittal";
            elseif twix.hdr.Spice.VoiNormalTra == 1
                hdr.Geometry.NormalPlane = "Transverse";
            end
    
            % Determine readout mode (if multiecho)
            if length(hdr.EchoTime) > 1
                if twix.hdr.Meas.ucReadOutMode == 1
                    hdr.ReadoutMode = "Monopolar";
                else
                    hdr.ReadoutMode = "Bipolar";
                end
            else
                hdr = rmfield(hdr, 'ReadoutMode');
            end
    
            % Isolate kspace data
            ksp = twix.image(:, ...  1) Columns
                             :, ...  2) Channels
                             :, ...  3) Lines (rows)
                             :, ...  4) Partitions (slices for 3D data)
                             :, ...  5) Slices (slabs for 3D data)
                             :, ...  6) Averages
                             1, ...  7) Cardiac Phases
                             :, ...  8) Contrasts
                             :, ...  9) Measurements
                             1, ... 10) Sets
                             :);  % 11) Segments

            % Swap columns and lines
            % dims = [3 2 1];
            % if ndims(ksp) > 3
            %     dims = [dims 4:ndims(ksp)]; %#ok<AGROW>
            % end
            % ksp = permute(ksp,dims);
    
            % Check dimensionality
            if strcmp(twix.hdr.Config.IceProgramName, '%SiemensIceProgs%\IceProgram3D') || strcmp(twix.hdr.Config.IceProgramName, '%SiemensIceProgs%\IceProgramOnline3D') || strcmp(twix.hdr.Config.SequenceFileName, '%SiemensSeq%\fl3d_vibe') % 3D
                dims = 3;
                hdr.DataMatrix = rmfield(hdr.DataMatrix, 'NumberOfSliceGroups');
                hdr.DataMatrix.NumberOfSlices = twix.image.NPar;
                hdr.DataMatrix.Size(3) = hdr.DataMatrix.NumberOfSlices;
            else
                hdr.DataMatrix = rmfield(hdr.DataMatrix, 'NumberOfSlabs');

                % Check for no partitions
                if hdr.DataMatrix.Size(3) == 1
                    hdr.DataMatrix.Size = [hdr.DataMatrix.Size(1:2) hdr.DataMatrix.Size(4:end)];

                    hdr.DataMatrix.NumberOfSliceGroups = 1;
                else
                    % Re-order partitions to later
                    hdr.DataMatrix.Size = [hdr.DataMatrix.Size(1:2) hdr.DataMatrix.Size(4:end-2) hdr.DataMatrix.Size(3) hdr.DataMatrix.Size(end-2:end)];
                end
                if twix.image.dataSize(5) > 1
                    dims = 2.5; % 2D multislice
                else
                    dims = 2;
                end
            end

            if ~isempty(twix.hdr.Config.Is3D)
                error('Maybe this can be used!')
            end
    
            % Check for os
            if rmosFLAG
                NCol = twix.image.NCol/twix.hdr.Dicom.flReadoutOSFactor;
            else
                NCol = twix.image.NCol;
            end
    
            % Check for zip
            NCol_og = NCol;
            if zipFLAG
                NCol = NCol*2;
                NLin = twix.image.NLin*2;
                hdr.DataMatrix.Size(1:2) = hdr.DataMatrix.Size(1:2)*2;
                if dims == 3
                    NPar = twix.image.NPar*2;
                    hdr.DataMatrix.Size(3) = hdr.DataMatrix.Size(3)*2;
                else
                    NPar = twix.image.NPar;
                end
            else
                NLin = twix.image.NLin;
                NPar = twix.image.NPar;
            end
    
            % Check for averaging
            if avgFLAG
                % Prevent loop that saves individual averages
                AVG = 1;
            else
                AVG = hdr.DataMatrix.NumberOfAverages;
            end
            hdr.DataMatrix.Size(5) = AVG;
    
            % Save data size info in header & preallocate data
            if dims == 3
                hdr.Geometry = rmfield(hdr.Geometry, 'DistanceFactor');
                hdr.DataMatrix.NumberOfSlices = NPar; 
                hdr.Geometry.SlabThickness = twix.hdr.MeasYaps.sSliceArray.asSlice{1, 1}.dThickness;
                hdr.Geometry.SliceThickness = hdr.Geometry.SlabThickness/hdr.DataMatrix.NumberOfSlices;
            else
                hdr.Geometry = rmfield(hdr.Geometry, 'SlabThickness');
                hdr.Geometry.SliceThickness = twix.hdr.MeasYaps.sSliceArray.asSlice{1, 1}.dThickness;

                % Check for distance factor
                if isfield(twix.hdr.MeasYaps.sGroupArray.asGroup{1},'dDistFact') && hdr.DataMatrix.NumberOfSlices > 1
                    hdr.Geometry.DistanceFactor = twix.hdr.MeasYaps.sGroupArray.asGroup{1}.dDistFact*hdr.Geometry.SliceThickness; % mm
                else
                    hdr.Geometry.DistanceFactor = 0;
                end
            end

            img_all = zeros(hdr.DataMatrix.Size, 'like', ksp);
            ksp_all = zeros(hdr.DataMatrix.Size, 'like', ksp);
            if hdr.DataMatrix.NumberOfChannels > 1
                img = zeros(hdr.DataMatrix.Size(1:end-1), 'like', ksp);
            else
                img = zeros(hdr.DataMatrix.Size, 'like', ksp);
            end
    
            % Obtain order for multislice/multislab data
            orderChron = twix.hdr.Config.chronSliceIndices(1:twix.image.NSli) + 1;
            orderSlice = twix.hdr.Config.relSliceNumber(1:twix.image.NSli) + 1;
            
            % Re-order twix data
            if dims == 3
                sSliceIdx = 1:hdr.DataMatrix.NumberOfSlabs;
                sfield = 'Slab';
                hdr.Geometry.Slab = hdr.Geometry.Slice;
                hdr.Geometry = rmfield(hdr.Geometry,'Slice');
            else
                sSliceIdx = 1:hdr.DataMatrix.NumberOfSlices;
                sfield = 'Slice';
            end
            asSlice = twix.hdr.MeasYaps.sSliceArray.asSlice{1};
            for sSl = sSliceIdx
                try
                    asSlice(orderSlice(sSl)) = twix.hdr.MeasYaps.sSliceArray.asSlice{orderChron(sSl)};
                catch
                    if ~isempty(twix.hdr.MeasYaps.sSliceArray.asSlice{orderChron(sSl)})
                        subfields = fieldnames(twix.hdr.MeasYaps.sSliceArray.asSlice{orderChron(sSl)});
                        for ff = 1:numel(subfields)
                            asSlice(orderSlice(sSl)).(subfields{ff}) = twix.hdr.MeasYaps.sSliceArray.asSlice{orderChron(sSl)}.(subfields{ff});
                        end
                    else
                        asSlice(orderSlice(sSl)) = NaN;
                    end
                end
                hdr.Geometry.(sfield)(orderSlice(sSl)).Offset = twix.image.slicePos(1:3,twix.image.NLin*(sSl-1) + 1);
                hdr.Geometry.(sfield)(orderSlice(sSl)).Quaternion = twix.image.slicePos(4:end,twix.image.NLin*(sSl-1) + 1);
            end

            % Grab geometry
            pcsDict = dictionary(["dSag", "dCor", "dTra"], 1:3);
            for sSl = sSliceIdx
                hdr.Geometry.(sfield)(sSl).FOV = [asSlice(sSl).dReadoutFOV asSlice(sSl).dPhaseFOV asSlice(sSl).dThickness];
                if dims == 3
                    hdr.Geometry.(sfield)(sSl).VoxelSize = hdr.Geometry.(sfield)(sSl).FOV./hdr.DataMatrix.Size(1:3);
                else
                    hdr.Geometry.(sfield)(sSl).VoxelSize = hdr.Geometry.(sfield)(sSl).FOV./[hdr.DataMatrix.Size(1:2) 1];
                end
                hdr.Geometry.(sfield)(sSl).Normal = [false false false];

                field = [];
                if isfield(asSlice(sSl).sNormal, "dSag")
                    field = [field "dSag"]; %#ok<AGROW>
                end
                if isfield(asSlice(sSl).sNormal, "dCor")
                    field = [field "dCor"]; %#ok<AGROW>
                end
                if isfield(asSlice(sSl).sNormal, "dTra")
                    field = [field "dTra"]; %#ok<AGROW>
                end
                
                % Try another attempt to find normal plane
                if strcmp(hdr.Geometry.NormalPlane,"")
                    switch field
                        case "dCor"
                            hdr.Geometry.NormalPlane = "Coronal";
                        case "dSag"
                            hdr.Geometry.NormalPlane = "Sagittal";
                        case "dTra"
                            hdr.Geometry.NormalPlane = "Transverse";
                    end
                end
                for f = 1:numel(field)
                    hdr.Geometry.(sfield)(sSl).Normal(pcsDict(field(f))) = true;

                    % % If field is not present, then it is NaN?
                    % if isfield(asSlice(sSl), "sPosition")
                    %     if ~isempty(asSlice(sSl).sPosition)
                    %         % If field is not present, then it is 0
                    %         if ~isfield(asSlice(sSl).sPosition, field(f))
                    %             hdr.Geometry.(sfield)(sSl).Offset(pcsDict(field(f))) = 0;
                    %         else
                    %             hdr.Geometry.(sfield)(sSl).Offset(pcsDict(field(f))) = asSlice(sSl).sPosition.(field(f));
                    %         end
                    %     end
                    % else
                    %     hdr.Geometry.(sfield)(sSl).Offset(pcsDict(field(f))) = NaN;
                    % end
                end
            end

            fprintf('Formatting data for MID%i... ',MID(m))
            denom = hdr.DataMatrix.NumberOfMeasurements*AVG*hdr.DataMatrix.NumberOfEchoes;

            % Obtain images
            if dims == 3
                denom = denom*hdr.DataMatrix.NumberOfSlabs;
                reverseStr = UpdatePercent((1/denom)*100, '');
                for meas = 1:hdr.DataMatrix.NumberOfMeasurements
                    for avg = 1:AVG
                        for eco = 1:hdr.DataMatrix.NumberOfEchoes
                            for slb = 1:hdr.DataMatrix.NumberOfSlabs
                                % Determine slab plane
                                if hdr.DataMatrix.NumberOfSlabs == 1
                                    hdr.SlabPlane = string(fieldnames(asSlice.sNormal));
                                else
                                    hdr.SlabPlane(slb) = string(fieldnames(asSlice{slb}.sNormal));
                                end

                                % Read in data
                                if avgFLAG
                                    % Read in all averages
                                    temp = sum(ksp(:,:,:,:,slb,:,:,eco,meas), 6)./size(ksp, 6);
                                else
                                    % Read in single average
                                    temp = ksp(:,:,:,:,slb,avg,:,eco,meas);
                                end

                                % Permute to correct shape & preallocate data
                                temp_ksp = permute(temp, [1 3 4 2]);

                                % ZIP (if needed)
                                if zipFLAG
                                    temp_ksp = padarray(temp_ksp, [NCol/4, NLin/4, hdr.DataMatrix.NumberOfSlices/4, 0], 0, 'both');
                                end

                                % Take fft
                                temp_img = temp_ksp;
                                for f = 1:size(temp_img,4)
                                    temp_img(:,:,:,f) = fftshift(fftn(ifftshift(temp_ksp(:,:,:,f))));
                                end

                                % Save full data
                                img_all(:,:,:,orderSlice(slb),eco,avg,meas,:) = temp_img;
                                ksp_all(:,:,:,orderSlice(slb),eco,avg,meas,:) = temp_ksp;

                                % Sum-of-square coil combination
                                img(:,:,:,orderSlice(slb),eco,avg,meas) = SWC(temp_img, temp_ksp, 3, zipFLAG);

                                reverseStr = UpdatePercent(((meas*avg*eco*slb)/denom)*100, reverseStr);
                            end
                        end
                    end
                end

            else % 2D
                denom = denom*hdr.DataMatrix.NumberOfSlices;
                reverseStr = UpdatePercent((1/denom)*100, '');
                for meas = 1:hdr.DataMatrix.NumberOfMeasurements
                    for avg = 1:AVG
                        for eco = 1:hdr.DataMatrix.NumberOfEchoes
                            for sli = 1:hdr.DataMatrix.NumberOfSlices

                                % Read in data
                                if avgFLAG
                                    % Read in all averages
                                    temp = sum(ksp(:,:,:,:,sli,:,:,eco,meas), 6)./size(ksp, 6);
                                else
                                    % Read in single average
                                    temp = ksp(:,:,:,:,sli,avg,:,eco,meas);
                                end

                                % Preallocate image data
                                temp_ksp = zeros(NCol_og, twix.image.NLin, twix.image.NCha, 'like', temp);

                                % Save kspace
                                for f = 1:twix.image.NCha
                                    temp_ksp(:,:,f) = squeeze(temp(:,f,:));
                                end

                                % ZIP (if needed)
                                if zipFLAG
                                    temp_ksp = padarray(temp_ksp, [NCol/4, NLin/4, 0], 0, 'both');
                                end

                                % Take fft
                                temp_img = temp_ksp;
                                temp_img = fftshift(fft2(ifftshift(temp_img)));

                                % Sum-of-square coil combination
                                if hdr.DataMatrix.NumberOfSlices == 1 % only 1 slice
                                    % Save full data
                                    img_all(:,:,:,eco,avg,meas,:) = temp_img;
                                    ksp_all(:,:,:,eco,avg,meas,:) = temp_ksp;

                                    % Coil combination
                                    % img(:,:,:,eco,avg,meas) = SoS(temp_img, 2);
                                    img(:,:,:,eco,avg,meas) = SWC(temp_img, temp_ksp, 2, zipFLAG);
                                else
                                    % Save full data
                                    img_all(:,:,orderSlice(sli),eco,avg,meas,:) = temp_img;
                                    ksp_all(:,:,orderSlice(sli),eco,avg,meas,:) = temp_ksp;

                                    % Coil combination
                                    img(:,:,orderSlice(sli),eco,avg,meas) = SWC(temp_img, temp_ksp, 2, zipFLAG);
                                end

                                reverseStr = UpdatePercent(((meas*avg*eco*sli)/denom)*100, reverseStr);
                            end
                        end
                    end
                end
            end

            % Remove empty dimensions
            hdr.DataMatrix.Size = hdr.DataMatrix.Size(hdr.DataMatrix.Size > 1);

            % Combine all into structure
            dat.hdr = hdr;
            dat.img = squeeze(img);
            dat.ksp = squeeze(ksp_all);
            dat.twix = twix;
    
            % Only add img_all if Rx is array
            if hdr.DataMatrix.NumberOfChannels ~= 1
                dat.img_all = squeeze(img_all);
            end
            
            fprintf('\n')
    
        %% Spectroscopy
        case {'%SiemensIceProgs%\IcePrgSpectroscopy', '%CustomerIceProgs%\ejaIcePrgSpec'}
            % Load new data
            twix = io_loadspec_twix(MID, 'KeepOS', ~rmosFLAG, 'ZIP', zipFLAG, 'twix', twix); % OR using io_loadspec_rda function from FID-A
    
            % Convert data to double
            twix.fids = double(twix.fids);
            twix.specs = double(twix.specs);
            
            % Rename data to 'dat' for function
            dat = twix;
    
            % Determine RF amplitude & if a noise scan
            if ~isfield(twix.hdr.MeasYaps.sTXSPEC.aRFPULSE{1, 1}, 'flAmplitude')
                dat.amp = 0;
                dat.flags.isNoiseScan = true;
            else
                dat.amp = twix.hdr.MeasYaps.sTXSPEC.aRFPULSE{1, 1}.flAmplitude; % amplitude
                dat.flags.isNoiseScan = false;
            end
    
            % Check for water suppression
            if strcmp(dat.nucleus, '1H')
                if strcmp(twix.hdr.Dicom.tScanOptions, 'WS')
                    dat.flags.isWaterSuppressed = true;
                    twix.flags.isWaterSuppressed = true;
                else
                    dat.flags.isWaterSuppressed = false;
                    twix.flags.isWaterSuppressed = false;
                end
            end

            % Check if data needs to be circularly shifted (sometimes the
            % first index in the spec data is somehow at the last index)
            % spec_comb = sum(sum(dat.specs,3),2);
            % if abs(angle(spec_comb(end-1))-angle(spec_comb(end))) > abs(angle(spec_comb(1))-angle(spec_comb(end)))*10
                disp('Data will be shifted left by one point!')
                dat.specs = circshift(dat.specs,1,1);
                dat.fids = op_ifft(dat.specs,1);
            % end
    
            % Add additional fields
            dat.fa = twix.hdr.Meas.FlipAngle; % flip angle
            dat.refamp = twix.hdr.Meas.TransmitterReferenceAmplitude; % reference amplitude
            % dat.refamp = twix_old.hdr.MeasYaps.sTXSPEC.asNucleusInfo{1, 1}.flReferenceAmplitude; % reference amplitude
            % dat.refamp = twix_old.hdr.Meas.flTransRefAmpl; % default reference amplitude
            dat.scantime = twix.hdr.MeasYaps.lTotalScanTimeSec; % scan time, sec
            dat.pdur = twix.hdr.Meas.MixingTime(1); % pulse duration, usec
            % dat.pdur = cell2mat(twix_old.hdr.MeasYaps.alTD); % pulse duration, usec
            if isfield(twix.hdr.MeasYaps.sSpecPara, 'lPreparingScans')
                dat.prepscans = twix.hdr.MeasYaps.sSpecPara.lPreparingScans;
            else
                dat.prepscans = 0;
            end
    
            if length(dat.sz) < 4 % not CSI
                %%%%%%%%%%%%%%% OUTPUT for McMRSGUI %%%%%%%%%%%%%%%
                if svFLAG
                    % Convert to McMRS
                    McMRSData = FIDa2McMRS(dat);
                    McMRSData.Filename = ['MID', num2str(MID), '_raw'];

                    % Create path
                    if dat.flags.isNoiseScan
                        pth = [pwd, '\Noise\', McMRSData.Filename, '.mat'];
                    elseif strcmp(dat.nucleus, '1H')
                        if dat.flags.isWaterSuppressed
                            pth = [pwd, '\Suppressed\', McMRSData.Filename, '.mat'];
                        else
                            pth = [pwd, '\Unsuppressed\', McMRSData.Filename, '.mat'];
                        end
                    else
                        pth = [pwd, '\', McMRSData.Filename, '.mat'];
                    end
                    fprintf('\nSaving McMRSData...')
                    save(pth, "McMRSData")
                    fprintf('Done!')
                    twix = dat; 
                    fprintf('\nSaving twix data...')
                    save([pth(1:end-4), '_twix.mat'], "twix")
                    fprintf('Done!\n')
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
    end

    % Combine into cell array
    dat_all{m} = dat;
end

% Re-organize MIDs
dat = dat_all;
if multiMID
    % Ensure fieldnames are all correct
    fnames = string(fieldnames(dat{1}));
    fnames_hdr = string(fieldnames(dat{1}.hdr));
    for m = 2:length(MID)
        fnames2 = string(fieldnames(dat{m}));
        fnames_hdr2 = string(fieldnames(dat{m}.hdr));
        if length(fnames) ~= length(fnames2) || length(fnames_hdr) ~= length(fnames_hdr2)
            m = 1;
            break
        elseif any(fnames ~= fnames2) || any(fnames_hdr ~= fnames_hdr2)
            m = 1;
            break
        end
    end

    if m > 1
        sz = size(dat{1});
        
        % Re-organize based on data size
        dat = dat{1};
        for m = 2:length(MID)
            switch length(sz)
                case 2
                    if sz(1) == 1 && sz(2) == 1
                        dat(:,m) = dat_all{m};
                    else
                        dat(:,:,m) = dat_all{m};
                    end
                case 3
                    dat(:,:,:,m) = dat_all{m};
                case 4
                    dat(:,:,:,:,m) = dat_all{m};
            end
        end
    end
else
    dat = dat{1};
end
end

% Sum of squares combination & saving into arrays
function cc = SoS(ic, nd)
% Inputs
%   ic: individual channel data (image data)
%   nd: number of dimensions
% 
% Outputs
%   cc: combined channel data

    % Check if only one channel
    if nd == length(size(ic))
        cc = ic;
    else
        % Sum-of-square coil combination:
        temp_mag = squeeze(sqrt(sum(abs(ic).^2, nd+1)));
        temp_phs = squeeze(sqrt(sum(angle(ic).^2, nd+1)));

        % Change to complex
        cc = temp_mag.*exp(1i*temp_phs);
    end
end

% SNR weighted data combination
function cc = SWC(ici, ick, nd, zip)
% Inputs
%  ici: individual channel data (image)
%  ick: individual channel data (kspace)
%   nd: number of dimensions
%  zip: check if data was zipped
% 
% Outputs
%   cc: combined channel data

    % Check if only one channel
    if nd == length(size(ici))
        cc = ici;       
    else
        % Isolate size of each dimension
        dims = size(ici);
        nc = dims(end); % number of coils
        dims = dims(1:end-1);

        % Determine sampling info
        Sc = floor(dims./2); % center
        Sr = floor(dims./10); % range size

        % Check if data was zipped & update data accordingly
        if zip
            Sr = floor(dims./20); % range size
            Smod = floor(dims./4); % border range modifier
        else
            Smod = ones(length(dims), 1);
        end

        % Measure signal & noise regions
        if nd == 2 % 2D
            s = squeeze(mean(abs(ick(Sc(1)-Sr(1)+1:Sc(1)+Sr(1), Sc(2)-Sr(2)+1:Sc(2)+Sr(2), :)), [1 2]));
            Sr = Sr.*2;
            n(:,1) = squeeze(mean(abs(ick(        (1:Sr(1))+Smod(1),         (1:Sr(2))+Smod(2), :)), [1 2]));
            n(:,2) = squeeze(mean(abs(ick((end+1-Sr(1):end)-Smod(1),         (1:Sr(2))+Smod(2), :)), [1 2]));
            n(:,3) = squeeze(mean(abs(ick(        (1:Sr(1))+Smod(1), (end+1-Sr(2):end)-Smod(2), :)), [1 2]));
            n(:,4) = squeeze(mean(abs(ick((end+1-Sr(1):end)-Smod(1), (end+1-Sr(2):end)-Smod(2), :)), [1 2]));
        else % 2D multislice or 3D
            s = squeeze(mean(abs(ick(Sc(1)-Sr(1)+1:Sc(1)+Sr(1), Sc(2)-Sr(2)+1:Sc(2)+Sr(2), Sc(3)-Sr(3)+1:Sc(3)+Sr(3), :)), [1 2 3]));
            Sr = Sr.*2;
            n(:,1) = squeeze(mean(abs(ick(        (1:Sr(1))+Smod(1),         (1:Sr(2))+Smod(2),         (1:Sr(3))+Smod(3), :)), [1 2 3]));
            n(:,2) = squeeze(mean(abs(ick((end+1-Sr(1):end)-Smod(1),         (1:Sr(2))+Smod(2),         (1:Sr(3))+Smod(3), :)), [1 2 3]));
            n(:,3) = squeeze(mean(abs(ick(        (1:Sr(1))+Smod(1), (end+1-Sr(2):end)-Smod(2),         (1:Sr(3))+Smod(3), :)), [1 2 3]));
            n(:,4) = squeeze(mean(abs(ick((end+1-Sr(1):end)-Smod(1), (end+1-Sr(2):end)-Smod(2),         (1:Sr(3))+Smod(3), :)), [1 2 3]));
            n(:,5) = squeeze(mean(abs(ick(        (1:Sr(1))+Smod(1),         (1:Sr(2))+Smod(2), (end+1-Sr(3):end)-Smod(3), :)), [1 2 3]));
            n(:,6) = squeeze(mean(abs(ick((end+1-Sr(1):end)-Smod(1),         (1:Sr(2))+Smod(2), (end+1-Sr(3):end)-Smod(3), :)), [1 2 3]));
            n(:,7) = squeeze(mean(abs(ick(        (1:Sr(1))+Smod(1), (end+1-Sr(2):end)-Smod(2), (end+1-Sr(3):end)-Smod(3), :)), [1 2 3]));
            n(:,8) = squeeze(mean(abs(ick((end+1-Sr(1):end)-Smod(1), (end+1-Sr(2):end)-Smod(2), (end+1-Sr(3):end)-Smod(3), :)), [1 2 3]));
        end
        n = sum(n, 2)./nc;

        % Calculate SNR
        SNR = s./sqrt(n);

        % Combine data
        if nd == 2 % 2D
            num = ici(:,:,1).*SNR(1);
            for i = 2:nc; num = num + ici(:,:,i).*SNR(i); end
        else % 2D multislice or 3D
            num = ici(:,:,:,1).*SNR(1);
            for i = 2:nc; num = num + ici(:,:,:,i).*SNR(i); end
        end
        cc = num./sqrt(sum(SNR.^2));
    end
end

function [rmosFLAG, zipFLAG, avgFLAG, svFLAG] = parseInputs(varargin)

% Set defaults
rmosFLAG = true;
zipFLAG = false;
avgFLAG = false;
svFLAG = false;

% Check for each input
CHECK_os = true;
CHECK_zip = true;
CHECK_avg = true;
CHECK_sv = true;
for i = 1:numel(varargin)
    % Check for OS input
    if strcmpi(varargin{i}, 'KeepOS') && CHECK_os
        % Remove input from options
        CHECK_os = false;

        % Set flag
        rmosFLAG = false;

    % Check for ZIP input
    elseif strcmpi(varargin{i}, 'ZIP') && CHECK_zip
        % Remove input from options
        CHECK_zip = false;

        % Set flag
        zipFLAG = true;

    % Check for CombAvg input
    elseif strcmpi(varargin{i}, 'CombAvg') && CHECK_avg
        % Remove input from options
        CHECK_avg = false;

        % Set flag
        avgFLAG = true;

   % Check for save input
    elseif strcmpi(varargin{i}, 'Save') && CHECK_sv
        % Remove input from options
        CHECK_sv = false;

        % Set flag
        svFLAG = true;
    end
end
end