function [] = check_v3B_wrf_profs(unit_test_data_dir)
%CHECK_V3B_WRF_PROFS Check that all swaths with changes in VCDs also use different WRF files

F = dir(fullfile(unit_test_data_dir, '*.mat'));

small_changes = [];
big_changes_diff_file = [];
big_changes_same_file = [];

for i_file = 1:numel(F)
    DataNew = load(fullfile(unit_test_data_dir, F(i_file).name),'Data');
    DataNew = DataNew.Data;
    behr_date = datenum(regexp(F(i_file).name, '\d\d\d\d\d\d\d\d', 'match', 'once'),'yyyymmdd');
    DataOld = load_behr_file(behr_date, 'daily', 'us', 'v3-0A');
    
    for i_swath = 1:numel(DataNew)
        max_del_vcd = max(abs(DataNew(i_swath).BEHRColumnAmountNO2Trop(:) - DataOld(i_swath).BEHRColumnAmountNO2Trop(:)));
        [~,new_wrf_file] = fileparts(DataNew(i_swath).BEHRWRFFile);
        [~,old_wrf_file] = fileparts(DataOld(i_swath).BEHRWRFFile);
        
        new_wrf_date = date_from_wrf_filenames(new_wrf_file);
        old_wrf_date = date_from_wrf_filenames(old_wrf_file);
        
        if max_del_vcd < 1e14
            small_changes = add_substruct(small_changes, behr_date, i_swath, max_del_vcd, old_wrf_file, new_wrf_file);
        elseif new_wrf_date == old_wrf_date
            big_changes_same_file = add_substruct(big_changes_same_file, behr_date, i_swath, max_del_vcd, old_wrf_file, new_wrf_file);
        else
            big_changes_diff_file = add_substruct(big_changes_diff_file, behr_date, i_swath, max_del_vcd, old_wrf_file, new_wrf_file);
        end
    end
end

fprintf('The following date/swaths have only small changes in VCDs:\n');
print_substructs(small_changes);
fprintf('The following date/swaths have big changes in VCDs but use different WRF files:\n');
print_substructs(big_changes_diff_file);
fprintf('Yikes! The following date/swaths have big changes in VCDs but use the same WRF files:\n');
print_substructs(big_changes_same_file)
end

function S = add_substruct(S, behr_date, behr_swath, max_diff, old_wrf_file, new_wrf_file)
Stemp = struct('date', behr_date, 'swath', behr_swath, 'del', max_diff, 'oldwrf', old_wrf_file, 'newwrf', new_wrf_file);
if isempty(S)
    S = Stemp;
else
    S(end+1) = Stemp;
end
end

function print_substructs(S)
for i = 1:numel(S)
    fprintf('    %s swath %d: max change = %g, old wrf = %s, new wrf = %s\n', datestr(S(i).date, 'yyyy-mm-dd'), S(i).swath, S(i).del, S(i).oldwrf, S(i).newwrf);
end
end