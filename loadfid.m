function s = loadfid(name,blocks,traces)

%% Description
%Usage [RE,IM,np,nb] = loadfid(name);

%Input:

%name   = name of FID directory without the .fid extension

%Output:

%RE  = real part of data, blocks are ordered columnwise
%IM  = imaginary part
%np  = number of points (rows in RE/IM); optional
%nb  = number of blocks (columns); optional
%nt  = number of traces; optional

%% Examples:

%[RE IM] = loadfid('filename');
%[RE IM np nb] = loadfid('filename');

%----------------------------------------
%       Xiaojun Chen Sept. 2007
%----------------------------------------

%% Format input name

fullname = sprintf('%s.fid%cfid',name,'/');
% fullname=sprintf(name);
% Compatibility with UNIX file system
% You don't need this step if you are given just one fid file without any
% extension
% Put fullname to that file location with .fid extension
fid = fopen(fullname,'r','ieee-be');
if fid == -1
    str = sprintf('Can not open file %s',fullname);
    error(str);
end

%% Read data file header
nblocks   = fread(fid,1,'int32');
ntraces   = fread(fid,1,'int32');
np        = fread(fid,1,'int32');
ebytes    = fread(fid,1,'int32');
tbytes    = fread(fid,1,'int32');
bbytes    = fread(fid,1,'int32');
vers_id   = fread(fid,1,'int16');
status    = fread(fid,1,'int16');
nbheaders = fread(fid,1,'int32');
s_data    = bitget(status,1);
s_spec    = bitget(status,2);
s_32      = bitget(status,3);
s_float   = bitget(status,4);
s_complex = bitget(status,5);
s_hyper   = bitget(status,6);

%% Reset output structures
% RE = [];
% IM = [];
RE=zeros(np/2,nblocks);
IM=zeros(np/2,nblocks);

if exist('blocks','var') == 0
    blocks = 1:nblocks;
end
outblocks   = max(size(blocks));

if exist('traces','var') == 0
    traces = 1:ntraces;
end
outtraces   = max(size(traces));

inx = 1;
B   = 1;

for b = 1:nblocks
    sprintf('read block %d\n',b);
    % Read a block header
    scale     = fread(fid,1,'int16');
    bstatus   = fread(fid,1,'int16');
    index     = fread(fid,1,'int16');
    mode      = fread(fid,1,'int16');
    ctcount   = fread(fid,1,'int32');
    lpval     = fread(fid,1,'float32');
    rpval     = fread(fid,1,'float32');
    lvl       = fread(fid,1,'float32');
    tlt       = fread(fid,1,'float32');
    T = 1;
    update_B = 0;
    for t = 1:ntraces
        if s_float == 1
            data = fread(fid,np,'float32');
            str = 'reading floats';
        elseif s_32 == 1
            data = fread(fid,np,'int32');
            str = 'reading 32bit';
        else
            data = fread(fid,np,'int16');
            str = 'reading 16bit';
        end
        if (blocks(B) == b)
            if T <= outtraces
                if (traces(T) == t)
                    sprintf('Reading block %d, trace %d, inx = %d\n',b,T,inx);
                    % Read k-space data in real and imaginary part
                    RE(:,inx) = data(1:2:np);
                    IM(:,inx) = data(2:2:np);
                    inx = inx + 1;
                    T = T + 1;
                    update_B = 1;
                end
            end 
        end 
    end 
    if update_B, B = B + 1; end
    if B > outblocks
        break
    end
end  
if nargout > 2
    NP = np/2;
end
if nargout > 3
    NB = nblocks;
end
if nargout > 4
    NT = ntraces;
end

fclose(fid);

%% Image processing
% K-space data 
s = complex(RE,IM);
% DC correction
% 10/20/2007
corr = mean(mean(s(end-round(0.2*length(s)):end)));
s = s - corr;