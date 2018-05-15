function [ F_behr, behr_dir ] = list_behr_files( start_date, end_date, prof_mode, region, version_str )
%LIST_BEHR_FILES List all BEHR .mat files between two dates
%   [ F, BEHR_DIR ] = LIST_BEHR_FILES( START_DATE, END_DATE ) will list all
%   monthly profile, US BEHR files between START_DATE and END_DATE, which
%   must be date specifications understood by Matlab (strings or numbers).
%   Returns F, a structure from DIR(), and BEHR_DIR, the directory it
%   searched.
%
%   [ ___ ] = LIST_BEHR_FILES( ___, PROF_MODE )
%   [ ___ ] = LIST_BEHR_FILES( ___, PROF_MODE, REGION ) allows you to use
%   different profile mode BEHR files ('daily' or 'monthly') and different
%   regions ('us', or 'hk').

E = JLLErrors;

start_date = validate_date(start_date);
end_date = validate_date(end_date);

if numel(start_date) ~= numel(end_date)
    E.badinput('Numel of elements in START_DATE and END_DATE must be equal');
end

if ~exist('prof_mode', 'var')
    prof_mode = 'monthly';
elseif ~ischar(prof_mode)
    E.badinput('PROF_MODE must be a string');
end
if ~exist('region', 'var')
    region = 'us';
elseif ~ischar(region)
    E.badinput('REGION must be a string');
end
if ~exist('version_str', 'var')
    version_str = BEHR_version();
end

F_behr = [];
file_pattern = behr_filename('*', prof_mode, region, '.mat', version_str);

for a=1:numel(start_date)
    if a == 1
        behr_dir = behr_paths.BEHRMatSubdir(region, prof_mode);
    else
        if ~strcmp(behr_dir, behr_paths.BEHRMatSubdir(region, prof_mode))
            E.notimplemented('LIST_BEHR_FILES cannot handle BEHR files across multiple directories');
        end
    end
    F = dir(fullfile(behr_dir, file_pattern));
    
    if isempty(F)
        continue
    end
    
    dvec = datenum(regexp({F.name}, '\d\d\d\d\d\d\d\d', 'match', 'once'), 'yyyymmdd');
    dd = dvec >= start_date(a) & dvec <= end_date(a);
    
    F(~dd) = [];
    
    F_behr = veccat(F_behr, F);

end

end

