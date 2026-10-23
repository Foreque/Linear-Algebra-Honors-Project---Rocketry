function [tc, Fc, m_prop] = readEng(filename)
% READENG  Read a RASP (.eng) motor thrust-curve file, e.g. from ThrustCurve.org
%
% Inputs:
%   filename - path to the .eng file (text), e.g. 'data/test_flat.eng'
%
% Outputs:
%   tc     - curve times, column vector [s], starts at 0, strictly increasing
%   Fc     - thrust at each time in tc, column vector [N]
%   m_prop - propellant mass from the header (field 5) [kg]
%
% File format (RASP):
%   ; lines starting with a semicolon are comments
%   NAME dia_mm len_mm delays prop_kg total_kg MFR   <- one header line
%   time thrust                                      <- data lines
%   ...
% The curve implicitly starts at (0,0); if the file's first time is > 0,
% that point is added here.

%% Open the file
fid = fopen(filename, 'r');          % fid = file ID (an integer handle)
if fid == -1
    error('readEng: could not open %s', filename)
end

%% Read line by line
gotHeader = false;                   % has the header line been read yet?
data      = zeros(0, 2);             % empty array, 2 columns: [time thrust]

while true
    ln = fgetl(fid);                 % next line as text, or -1 at end of file
    if ~ischar(ln), break, end       % -1 is a number, not text -> end of file
    ln = strtrim(ln);                % remove leading/trailing spaces

    % skip blank lines and comment lines
    if isempty(ln) || startsWith(ln, ';'), continue, end

    if ~gotHeader
        % ---- first real line: the header ----
        parts  = strsplit(ln);           % split on spaces -> cell array of text
        m_prop = str2double(parts{5});   % field 5 = propellant mass, kg
        gotHeader = true;
    else
        % ---- every line after: one data point ----
        row = sscanf(ln, '%f')';         % numbers in the line, as a row [t F]
        data(end+1, :) = row(1:2);       %#ok<AGROW> append as a new row
    end
end
fclose(fid);                         % release the file

%% Checks and cleanup
if ~gotHeader || isempty(data)
    error('readEng: no header or no data found in %s', filename)
end
if isnan(m_prop)
    error('readEng: header field 5 (propellant mass) is not a number')
end

tc = data(:, 1);                     % times, column
Fc = data(:, 2);                     % thrusts, column

% add the implied (0,0) starting point
if tc(1) > 0
    tc = [0; tc];                    % semicolon = stack vertically (columns)
    Fc = [0; Fc];
end

% interp1 needs strictly increasing times; catch problems here
if any(diff(tc) <= 0)
    error('readEng: times in %s are duplicated or out of order', filename)
end
