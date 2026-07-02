function [data, pth] = loadmrui(varargin)
%LOADMRUI Load .mrui data from jMRUI spectroscopy software package
% 
%   [data, pth] = loadmrui(file)
% 
%   Inputs
%          file: *optional* .mrui file, with or without full path
%                   - if full path is not included, file must be in current directory
% 
%   Outputs
%          data: structure containing spectroscopy data & header information
%           pth: source .mrui file path
%
% Jacob Degitz, Texas A&M University
% Created 12/17/2024
% Last edited 12/18/2024

% Parse inputs
pth = parseInputs(varargin{:});

% Open file
fid = fopen(pth, 'r', 'ieee-be.l64');

% Read header
hdr = fread(fid,13,'float64',0,'ieee-be.l64');

% Isolate number of points
N = hdr(2);
fseek(fid,2,'cof'); % skip two bytes from current position
n_frames = -1*round(floor((ftell(fid) - 512)/(N*16))*N*2); % determine number of bytes corresponding to data
fseek(fid,512,'bof'); % move to start of data

% Loop through datasets
multFLAG = false;
idx = 0;
while true
    % Save location
    loc = ftell(fid);

    % Read data
    tmp = fread(fid,n_frames,'float64',0,'ieee-be.l64');

    % Check if next section is binary or string
    % - this is the case if the length of tmp is less than the length of n_frames
    if length(tmp) < n_frames
        % Reset position
        fseek(fid,loc,'bof');

        % Extract string data
        if ~multFLAG
            name = fgetl(fid);
            date = fgetl(fid);
        else
            name = NaN;
            date = NaN;
            fgetl(fid);
            fgetl(fid);
            fgetl(fid);
            fgetl(fid);
            mets = strsplit(string(fgetl(fid)),';');
            if mets(end) == "", mets = mets(1:end-1); end
        end

        break
    else
        % Append
        idx = idx+1;
        data(:,idx) = tmp;

        if idx > 1 && ~multFLAG, multFLAG = true; end
    end
end

% Extract data from header & string array
sz = [hdr(2), 1];
dwelltime = hdr(3)/1000; % seconds
spectralwidth = 1/dwelltime; % Hz
txfrq = hdr(6); % Hz
dims = struct('t', 1, 'coils', 0, 'averages', 0, 'subSpecs', 0, 'extras', 0);
Bo = hdr(7); % T
[pth, fname, ext] = fileparts(pth); fname = char(append(fname, ext));
flags = struct('apodized', false, 'zeropadded', false);
if hdr(12) ~= 0, flags.apodized = true; end
if hdr(13) ~= 0, flags.zeropadded = true; end
a1 = hdr(5);
a0 = hdr(4);
reffrq = hdr(9);
refppm = hdr(10);

% Instead, find nucleus via B0 & F0
gamma = floor((txfrq/Bo)*1e-6);
switch gamma
    case 43
        nucleus = '1H';
        gamma = 42.5774780505984; % MHz/T
        ppm_off = 4.65;
    case 17
        nucleus = '31P';
        gamma = 17.2514528352478;
        ppm_off = 0;
    case 10
        nucleus = '13C';
        gamma = 10.7083987615955;
        ppm_off = 0;
    case 40
        nucleus = '19F';
        gamma = 40.0775824603147;
        ppm_off = 0;
    case 11
        nucleus = '23Na';
        gamma = 11.2688453499836;
        ppm_off = 0;
    otherwise
        nucleus = 'unknown';
        gamma = NaN;
        ppm_off = 0;
end

% Determine ID, if applicable
if ~isnan(double(string(fname(4:6))))
    ID = fname(4:6);
else
    if ~isnan(double(string(fname(4:6))))
        ID = fname(4:5);
    else
        ID = NaN;
    end
end

% Convert data to complex
if hdr(11) == 0 % fid
    fids = data(1:2:end,:) + 1i.*data(2:2:end,:);
    specs = zeros(size(fids));
    for i = 1:size(specs,2)
        specs(:,i) = fftshift(fft(fids(:,i)));
    end
else % echo
    specs = data(1:2:end,:) + 1i.*data(2:2:end,:);
    fids = zeros(size(specs));
    for i = 1:size(fids,2)
        fids(:,i) = ifft(ifftshift(specs(:,i)));
    end
end

% For mult data, normalize so they all start at 1
if multFLAG
    for i = 1:size(data,2)
        % Isolate real & imaginary
        specs_re = real(specs(:,i));
        specs_im = imag(specs(:,i));
        
        % Normalize
        specs_re = specs_re - specs_re(1);
        specs_im = specs_im - specs_im(1);

        % Combine & grab fid
        specs(:,i) = specs_re + 1i.*specs_im;
        fids(:,i) = ifft(ifftshift(specs(:,i)));
        clear specs_re specs_im
    end
end

% Create x-axis arrays
t = 0:dwelltime:(sz(1)-1)*dwelltime;
f = (-spectralwidth/2):spectralwidth/(sz(1)):(spectralwidth/2)-1*spectralwidth/sz(1); 
ppm = f/(Bo*gamma) - ppm_off;

% Combine into structure
data = struct('fids', fids(:,1), ...
              'specs', specs(:,1), ...
              'sz', sz, ...
              'ppm', ppm, ...
              't', t, ...
              'spectralwidth', spectralwidth, ...
              'dwelltime', dwelltime, ...
              'txfrq', txfrq, ...
              'date', date, ...
              'dims', dims, ...
              'Bo', Bo, ...
              'nucleus', nucleus, ...
              'gamma', gamma, ...
              'fname', fname, ...
              'flags', flags, ...
              'ID', ID, ...
              'name', name, ...
              'a1', a1, ...
              'a0', a0, ...
              'reffrq', reffrq, ...
              'refppm', refppm, ...
              'metabolite', char(mets(1)));
dtmp = data;
for i = 2:length(mets)
    data(i) = dtmp;
    data(i).fids = fids(:,i);
    data(i).specs = specs(:,i);
    data(i).metabolite = char(mets(i));
end

fclose(fid);

    function pth = parseInputs(varargin)
        % Check if no inputs
        if nargin == 0   
            % Extract path via UI input
            [file, pth] = uigetfile('*.mrui');
            pth = [pth file];

            % Check if input is a .mrui file
            if ~strcmp(pth(end-4:end),'.mrui'), error('Input must be an .mrui file'); end
        else
            % Check if input is string or character
            if isstring(varargin{1}) || ischar(varargin{1})
                in = char(varargin{1});
            else
                error('Input must be a string or character vector')
            end

            % Check if full path or just file
            if strcmp(in(1:2),'C:')
                pth = in;
            else
                % Check if input has file type
                if any(in == '.')
                    % Check if file is correct type
                    if strcmp(in(end-4:end),'.mrui')
                        pth = [pwd '\' in];
                    else
                        error('Input must be an .mrui file')
                    end
                else 
                    pth = [pwd '\' in '.mrui'];
                end
            end

            % Check if file exists
            if ~exist(pth, 'file'), error('File could not be identified'); end
        end
    end
end