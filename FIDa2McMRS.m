function [McMRSData, McMRSDataw] = FIDa2McMRS(twix, varargin)
% FIDA2MCMRS converts twix structures from mapVBVD to McMRS structures
%
% Starting = FIDa2McMRS(twix) will convert the twix structure to a McMRS
% Starting structure and save in the current working directory. The filename
% will be automatically generated based on the measurement ID along with 
% flags present in the twix structure. 
%
% [McMRSData, McMRSDataw] = FIDa2McMRS(twix, 'Water', twixw) will convert the twix
% structure and corresponding water twix to McMRS data structures.
%
% McMRSData = FIDa2McMRS(twix, 'Save', saveFLAG) will convert the twix structure to 
% a McMRS data structure, but will only save the file if the boolean 
% saveFLAG is true.

% Generate input parser
p = inputParser;
addRequired(p, 'twix', @(x) isstruct(x))
addParameter(p, 'Water', [], @(x) isstruct(x))
addParameter(p, 'Save', true, @(x) (isscalar(x) && islogical(x)))
parse(p, twix, varargin{:})

% Parse inputs
twix = p.Results.twix;
saveFLAG = p.Results.Save;
if isstruct(p.Results.Water)
    waterFLAG = true;
else
    waterFLAG = false;
end

%{
Change Dimensions
fidA structure is generally [time, coils, averages, whatever else]
McMRS structure is [coils, averages, time]
%}
dimF = fields(twix.dims);
for d = 1:numel(dimF)
    if twix.dims.(dimF{d})~=0 && ~any(strcmp(dimF{d}, {'t', 'coils', 'averages'}))
        error('McMRS cannot handle data with extra or subSpec dimensions.')
    end
end

if twix.dims.coils==0 && twix.dims.averages==0 % single-channel, single average
    fids = reshape(permute(twix.fids, [2 twix.dims.t]), [1 1 twix.sz(twix.dims.t)]);
    specs = reshape(permute(twix.specs, [2 twix.dims.t]), [1 1 twix.sz(twix.dims.t)]);
    if waterFLAG
        fidsw = reshape(permute(twixw.fids, [2 twixw.dims.t]), [1 1 twixw.sz(twixw.dims.t)]);
        specsw = reshape(permute(twixw.specs, [2 twixw.dims.t]), [1 1 twixw.sz(twixw.dims.t)]);
    end
elseif twix.dims.coils~=0 && twix.dims.averages==0 % multi-channel, single average
    fids = reshape(permute(twix.fids, [twix.dims.coils twix.dims.t]), [twix.sz(twix.dims.coils) 1 twix.sz(twix.dims.t)]);
    specs = reshape(permute(twix.specs, [twix.dims.coils twix.dims.t]), [twix.sz(twix.dims.coils) 1 twix.sz(twix.dims.t)]);
    if waterFLAG
        fidsw = reshape(permute(twixw.fids, [twixw.dims.coils twixw.dims.t]), [twixw.sz(twixw.dims.coils) 1 twixw.sz(twixw.dims.t)]);
        specsw = reshape(permute(twixw.specs, [twixw.dims.coils twixw.dims.t]), [twixw.sz(twixw.dims.coils) 1 twixw.sz(twixw.dims.t)]);
    end
elseif twix.dims.coils==0 && twix.dims.averages~=0 % single-channel, multiple averages
    fids = reshape(permute(twix.fids, [twix.dims.averages twix.dims.t]), [1 twix.sz(twix.dims.averages) twix.sz(twix.dims.t)]);
    specs = reshape(permute(twix.specs, [twix.dims.averages twix.dims.t]), [1 twix.sz(twix.dims.averages) twix.sz(twix.dims.t)]);
    if waterFLAG
        fidsw = reshape(permute(twixw.fids, [twix.dims.averages twix.dims.t]), [1 twixw.sz(twixw.dims.averages) twixw.sz(twixw.dims.t)]);
        specsw = reshape(permute(twixw.specs, [twix.dims.averages twix.dims.t]), [1 twixw.sz(twixw.dims.averages) twixw.sz(twixw.dims.t)]);
    end
elseif twix.dims.coils~=0 && twix.dims.averages~=0 % multi-channel, multiple averages
    fids = permute(twix.fids, [twix.dims.coils twix.dims.averages twix.dims.t]);
    specs = permute(twix.specs, [twix.dims.coils twix.dims.averages twix.dims.t]);
    if waterFLAG
        fidsw = permute(twixw.fids, [twix.dims.coils twix.dims.averages twix.dims.t]);
        specsw = permute(twixw.specs, [twix.dims.coils twix.dims.averages twix.dims.t]);
    end
end
specs = flip(specs,3);
ppm = flip(twix.ppm);
if waterFLAG
    specsw = flip(specsw,3);
    ppmw = flip(twixw.ppm);
end

% Grab phase
phShift = [];
if waterFLAG
    phShiftw = [];
end
if isfield(twix,'results')
    if isfield(twix.results.phasecorrected,'phShift')
        phShift = twix.results.phasecorrected.phShift.';
        if waterFLAG
            phShiftw = twixw.results.phasecorrected.phShift.';
        end
    end
    if isfield(twix.results.freqcorrected,'phsCum') && ~isempty(twix.results.freqcorrected.phsCum) && any(twix.results.freqcorrected.phsCum ~= 0)
        phShift = phShift - twix.results.freqcorrected.phsCum.';
        if waterFLAG
            phShiftw = phShiftw - twixw.results.freqcorrected.phsCum.';
        end
    end
end

% Grab number of channels & averages
if twix.dims.coils ~= 0
    nc = twix.sz(twix.dims.coils);
    C = 1:nc;
else
    nc = 1;
    C = 1;
end
if twix.dims.averages ~= 0
    na = twix.sz(twix.dims.averages);
    A = 1:twix.rawAverages;
    if isfield(twix, 'results')
        Adx = true(size(A)); 
        Adx(twix.results.rmbadaverages.badaverages) = false;
        A = A(Adx);
    end
else
    na = 1;
    A = 1;
end
if isfield(twix, 'rawCoils')
    rC = twix.rawCoils;
else
    rC = twix.sz(twix.dims.coils);
end

% Grab apodization
if twix.flags.filtered
    lb = twix.results.filtered.lb;
else
    lb = [];
end

% Grab zero padding
if twix.flags.zeropadded
    zp = twix.results.zeropadded.zpFactor-1;
else
    zp = [];
end

% Grab frequency shift
if isfield(twix,'results') && isfield(twix.results.freqcorrected,'fsCum') && any(twix.results.freqcorrected.fsCum ~= 0,'all')
    fsh = twix.results.freqcorrected.fsCum.';
else
    fsh = [];
end

% Find peak location & noise limits
[~,ploc] = max(real(sum(specs,2)),[],3);
ploc = ppm(ploc);
switch twix.nucleus
    case '1H'
        nrloc = [-3 0];
    case '31P'
        nrloc = [-30 -20];
    otherwise
        nrloc = [];
end
if waterFLAG
    [~,plocw] = max(real(sum(specsw,2)),[],3);
    plocw = ppm(plocw);
end

% Set ref ppm point
if strcmp(twix.nucleus,'1H')
    rppm = 4.65;
else
    rppm = 0;
end

% Grab measurement ID
fname = strsplit(twix.fname, '\');
fname = fname{end};
fname = strsplit(fname,'_');
fname = fname{2};
MID = '';
MID_text = 'MID';
idx_fname = 0;
idx_MID = 0;
while idx_fname < numel(fname)
    idx_fname = idx_fname + 1;

    if idx_MID < 3
        % Move through 'MID'
        if strcmp(fname(idx_fname), MID_text(idx_MID+1))
            idx_MID = idx_MID + 1;
        else
            idx_MID = 0; % reset
        end
    else
        % Append numeric character
        if ~isnan(str2double(fname(idx_fname)))
            MID = [MID fname(idx_fname)]; %#ok<AGROW>
        else
            break
        end
    end
end

% Trim off excess zeros
MID = char(string(str2double(MID)));

% Generate name
ch = "";
nameend = "";
if twix.flags.leftshifted || twix.flags.freqcorrected || any(twix.flags.phasecorrected)
    nameend = append(nameend,"_preproc");
end
if twix.flags.filtered
    nameend = append(nameend,"_filt");
end
if twix.flags.zeropadded || (isfield(twix.flags, 'zeropaddedspec') && twix.flags.zeropaddedspec)
    nameend = append(nameend,"_zp");
end
if twix.flags.averaged
    nameend = append(nameend,"_avg");
end
if twix.flags.addedrcvrs
    nameend = append(nameend,"_cc");
end

% Save to 'Starting' structure
McMRSData = struct('Filename', append("MID",MID,ch,nameend), ...
                'SourceFiles', "", ... 
                   'Spectrum', specs, ...
                 'TimeDomain', fids, ... 
                    'PPMAxis', ppm, ... 
              'MagneticField', twix.Bo, ... 
                    'Nucleus', twix.nucleus, ...
                  'Frequency', twix.txfrq/1e6, ... 
              'SpectralWidth', twix.spectralwidth, ... 
               'ReferencePPM', rppm, ...
                   'Channels', C, ...
                   'Averages', A, ...
                'RawChannels', rC, ... 
                'RawAverages', twix.rawAverages, ... 
                       'Size', struct('Channels', nc, ... 
                                      'Averages', na, ... 
                                        'Points', twix.sz(twix.dims.t)), ...
               'PeakLocation', ploc, ... 
          'NoiseRegionLimits', nrloc, ... 
                      'Flags', struct('WaterSuppressed', twix.flags.isWaterSuppressed), ...
                 'Operations', struct('ZeroOrderPhase', phShift, ... 
                                     'FirstOrderPhase', [], ... 
                                         'Apodization', lb, ...
                                         'ZeroPadding', zp, ... 
                                      'FrequencyShift', fsh, ...
                                  'BaselineCorrection', [], ... 
                                  'ChannelCombination', struct('Method', ""), ...
                                          'SourceFile', struct('ZeroOrderPhase', [], ... 
                                                              'FirstOrderPhase', [], ... 
                                                                  'Apodization', [], ...
                                                                  'ZeroPadding', [], ... 
                                                               'FrequencyShift', [], ...
                                                           'BaselineCorrection', [])));

if waterFLAG
    McMRSData = struct('Filename', append("MID",MID,ch,nameend,"_w"), ...
                    'SourceFiles', "", ... 
                       'Spectrum', specsw, ...
                     'TimeDomain', fidsw, ... 
                        'PPMAxis', ppmw, ... 
                  'MagneticField', twixw.Bo, ... 
                        'Nucleus', twixw.nucleus, ...
                      'Frequency', twixw.txfrq/1e6, ... 
                  'SpectralWidth', twixw.spectralwidth, ... 
                       'Channels', C, ...
                       'Averages', A, ...
                    'RawChannels', twixw.rawCoils, ... 
                    'RawAverages', twixw.rawAverages, ... 
                           'Size', struct('Channels', nc, ... 
                                          'Averages', na, ... 
                                            'Points', twixw.sz(twixw.dims.t)), ...
                   'PeakLocation', plocw, ... 
              'NoiseRegionLimits', nrloc, ... 
                     'Operations', struct('ZeroOrderPhase', phShiftw, ... 
                                         'FirstOrderPhase', [], ... 
                                             'Apodization', lb, ...
                                             'ZeroPadding', zp, ... 
                                          'FrequencyShift', fsh, ...
                                      'BaselineCorrection', [], ... 
                                      'ChannelCombination', struct('Method', ""), ...
                                              'SourceFile', struct('ZeroOrderPhase', [], ... 
                                                                  'FirstOrderPhase', [], ... 
                                                                      'Apodization', [], ...
                                                                      'ZeroPadding', [], ... 
                                                                   'FrequencyShift', [], ...
                                                               'BaselineCorrection', [])));
else
    McMRSDataw = [];
end

% Split channel data
if saveFLAG
    save(McMRSData.Filename, "McMRSData")
    if waterFLAG
        Starting_tmp = McMRSData;
        McMRSData = McMRSDataw;
        save(McMRSData.Filename, "McMRSData")
        McMRSData = Starting_tmp;
    end
end