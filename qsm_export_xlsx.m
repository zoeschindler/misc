function qsm_export( ...
    input_path, ... % path to the QSM *.mat-files, necessary
    output_path, ... % path to the result folder, optional
    savename, ... % name of the file, optional
    txt_files, ... % boolean, should default text files be saved?
    xlsx_file, ... % boolean, should data be saved in xlsx file?
    xlsx_cylinder, ... % boolean, should binned cylinder data be saved?
    xlsx_branch) % boolean, should binned branch data be saved?

    %% CHECK INPUT & READ QSM

    % check for output path
    if ~exist('output_path','var') || output_path == ""
      output_path = fileparts(input_path);
    end
    
    % load QSM matlab file
    load(input_path, "QSM");  
    
    % check for savename
    if ~exist('savename','var') || output_path == ""
      savename = QSM.rundata.inputs.name;
    end

    % check for txt_files
    if ~exist('txt_files','var'); txt_files = false; end

    % check for xlsx_file
    if ~exist('xlsx_file','var'); xlsx_file = true; end

    % check for xlsx_cylinder
    if ~exist('xlsx_cylinder','var'); xlsx_cylinder = true; end

    % check for xlsx_branch
    if ~exist('xlsx_branch','var'); xlsx_branch = true; end

    
    %% PREPARE DATA
    % as in TreeQSM/src/tools/save_model_text.m

    % extract data from QSM object
    cylinder = QSM.cylinder;
    branch = QSM.branch;
    treedata = QSM.treedata;
    
    % form & round cylinder data
    Rad = round(10000*cylinder.radius)/10000; % radius (m)
    Len = round(10000*cylinder.length)/10000; % length (m)
    Sta = round(10000*cylinder.start)/10000; % starting point (m)
    Axe = round(10000*cylinder.axis)/10000; % axis (m)
    CPar = single(cylinder.parent); % parent cylinder
    CExt = single(cylinder.extension); % extension cylinder
    Added = single(cylinder.added); % is cylinder added to fil a gap
    Rad0 = round(10000*cylinder.UnmodRadius)/10000; % unmodified radius (m)
    B = single(cylinder.branch); % branch index of the cylinder
    BO = single(cylinder.BranchOrder); % branch order of the branch
    PIB = single(cylinder.PositionInBranch); % position of the cyl. in the branch
    Mad = single(round(10000*cylinder.mad)/10000); % mean abso. distance (m)
    SC = single(round(10000*cylinder.SurfCov)/10000); % surface coverage
    CylData = [Rad Len Sta Axe CPar CExt B BO PIB Mad SC Added Rad0];
    NamesC = ["radius_m","length_m","start_X", "start_Y","start_Z", ...
      "axis_X","axis_Y","axis_Z","parent","extension","branch", ...
      "branch_order","position_in_branch","mad","SurfCov","added", ...
      "UnmodRadius_m"];
    
    % form & round branch data
    BOrd = single(branch.order); % branch order
    BPar = single(branch.parent); % parent branch
    BDia = round(10000*branch.diameter)/10000; % diameter (m)
    BVol = round(10000*branch.volume)/10000; % volume (L)
    BAre = round(10000*branch.area)/10000; % area (m^2)
    BLen = round(1000*branch.length)/1000; % length (m)
    BAng = round(10*branch.angle)/10; % angle (deg)
    BHei = round(1000*branch.height)/1000; % height (m)
    BAzi = round(10*branch.azimuth)/10; % azimuth (deg)
    BZen = round(10*branch.zenith)/10; % zenith (deg)
    BranchData = [BOrd BPar BDia BVol BAre BLen BHei BAng BAzi BZen];
    NamesB = ["order","parent","diameter_m","volume_L","area_m^2", ...
      "length_m","height_m","angle_deg","azimuth_deg","zenith_deg"];
    
    % extract the field names of treedata
    Names = fieldnames(treedata);
    
    % get indices of relevant data
    n = 1;
    % increases index n until "location" is reached
    while ~strcmp(Names{n},'location')
        n = n + 1;
    end
    n = n - 1; % exclude "location", as it has three numbers, not one
    Names = Names(1:n);
    
    % create empty variable for storage
    TreeData = zeros(n,1); 
    
    % store data in new variable
    for i = 1:n
        TreeData(i) = treedata.(Names{i,:});
    end
    TreeData = change_precision(TreeData); % use less decimals
    NamesD = string(Names);
    
    %% SAVE DATA AS TXT-FILES
    % as in TreeQSM/src/tools/save_model_text.m
    
    % only if txt-files should be saved
    if txt_files
    
        % save cylinder data
        str = output_path + '/cylinder_' + savename +'.txt';
        fid = fopen(str, 'wt');
        fprintf(fid, [repmat('%s\t', 1, size(NamesC,2)-1) '%s\n'], NamesC.');
        fprintf(fid, [repmat('%g\t', 1, size(CylData,2)-1) '%g\n'], CylData.');
        fclose(fid);
        
        % save branch data
        str = output_path + '/branch_' + savename +'.txt';
        fid = fopen(str, 'wt');
        fprintf(fid, [repmat('%s\t', 1, size(NamesB,2)-1) '%s\n'], NamesB.');
        fprintf(fid, [repmat('%g\t', 1, size(BranchData,2)-1) '%g\n'], BranchData.');
        fclose(fid);
        
        % save tree data
        str = output_path + '/treedata_' + savename + '.txt';
        fid = fopen(str, 'wt');
        NamesD(:,2) = TreeData;
        fprintf(fid,'%s\t %g\n',NamesD.');
        fclose(fid);
    
    end
    
    %% SAVE DATA AS XLSX-FILE
    
    % only if xlsx-file should be saved
    if xlsx_file
    
        % set path to xlsx file
        str = output_path + '/QSM_data_' + savename + '.xlsx';
        
        % export input variables
        input_values = struct2cell(QSM.rundata.inputs);
        input_names = fieldnames(QSM.rundata.inputs);
        writetable(table(input_names,input_values, ...
            'VariableNames',{'Parameter','Value'}),str, ...
            'WriteRowNames',true, ...
            'Sheet','input_parameters');
        
        % export cylinder data
        cylinder_table = struct2table(QSM.cylinder);
        cylinder_table = splitvars(cylinder_table,{'start','axis'}, ...
            'NewVariableNames',{{'start_X','start_Y','start_Z'}, ...
            {'axis_X','axis_Y','axis_Z'}}); % split 'start' & 'axis' multicolumns
        writetable(cylinder_table,str, ...
            'Sheet','cylinder');
        
        % export branch data
        branch_table = struct2table(QSM.branch);
        writetable(branch_table,str, ...
            'Sheet','branch');
        
        % extract treedata
        treedata_values = struct2cell(QSM.treedata);
        treedata_names = fieldnames(QSM.treedata);
        
        % export treedata - overview
        treedata_values{n+1} = strjoin(replace(string(QSM.treedata.location), ...
            '.',','),"; "); % location -> string
        treedata_table = cell2table(treedata_values, ...
            'RowNames', treedata_names, ...
            'VariableNames',{'Value'});
        writetable(treedata_table(1:n+1,:), str, ...
            'WriteRowNames',true, ...
            'Sheet','treedata_overview');
        
        % export treedata - BranchOrder
        treedata_branchorder_table = table((1:numel(QSM.treedata.NumBranchOrd))', ...
            QSM.treedata.NumBranchOrd',QSM.treedata.LenBranchOrd', ...
            QSM.treedata.AreBranchOrd',QSM.treedata.VolBranchOrd', ...
            'VariableNames',{'BranchOrder','NumBranchOrd','LenBranchOrd', ...
            'AreBranchOrd','VolBranchOrd'});
        writetable(treedata_branchorder_table, str, ...
            'Sheet','treedata_BranchOrder');
        
        % export treedata - StemTaper
        % treedata.StemTaper + Cylinder.start + Cylinder.axis
        stem_taper_indices = 1:(width(QSM.treedata.StemTaper)-1);
        stem_taper_table = table( ...
            QSM.treedata.StemTaper(1,:)',QSM.treedata.StemTaper(2,:)', ...
            vertcat(QSM.cylinder.start(stem_taper_indices,1), 0), ...
            vertcat(QSM.cylinder.start(stem_taper_indices,2), 0), ...
            vertcat(QSM.cylinder.start(stem_taper_indices,3), 0), ...
            vertcat(QSM.cylinder.axis(stem_taper_indices,1), 0), ...
            vertcat(QSM.cylinder.axis(stem_taper_indices,2), 0), ...
            vertcat(QSM.cylinder.axis(stem_taper_indices,3), 0), ...
            'VariableNames', {'distance_m', 'diameter_m', ...
            'start_X','start_Y','start_Z','axis_X','axis_Y','axis_Z'});
        writetable(stem_taper_table, str, ...
            'Sheet','treedata_StemTaper');
        
        % only if cylinder classes should be saved
        if xlsx_cylinder
    
            % export treedata - Diameter classes, all cylinders
            treedata_cyl_dia_table = table(QSM.treedata.LenCylDia', ...
                QSM.treedata.AreCylDia', QSM.treedata.VolCylDia', ...
                'VariableNames',{'LenCylDia', 'AreCylDia', 'VolCylDia'});
            writetable(treedata_cyl_dia_table, str, ...
                'Sheet','treedata_CylDia');
            
            % export treedata - Height classes, all cylinders
            treedata_cyl_hei_table = table(QSM.treedata.LenCylHei', ...
                QSM.treedata.AreCylHei', QSM.treedata.VolCylHei', ...
                'VariableNames',{'LenCylHei', 'AreCylHei', 'VolCylHei'});
            writetable(treedata_cyl_hei_table, str, ...
                'Sheet','treedata_CylHei');
            
            % export treedata - Zenith classes, all cylinders
            treedata_cyl_zen_table = table(QSM.treedata.LenCylZen', ...
                QSM.treedata.AreCylZen', QSM.treedata.VolCylZen', ...
                'VariableNames',{'LenCylZen', 'AreCylZen', 'VolCylZen'});
            writetable(treedata_cyl_zen_table, str, ...
                'Sheet','treedata_CylZen');
            
            % export treedata - Azimuth classes, all cylinders
            treedata_cyl_azi_table = table(QSM.treedata.LenCylAzi', ...
                QSM.treedata.AreCylAzi', QSM.treedata.VolCylAzi', ...
                'VariableNames',{'LenCylAzi', 'AreCylAzi', 'VolCylAzi'});
            writetable(treedata_cyl_azi_table, str, ...
                'Sheet','treedata_CylAzi');
            
        end
    
        % only if branch classes should be saved
        if xlsx_branch
    
            % export treedata - Diameter classes, branches
            treedata_cyl_dia_table = table( ...
                QSM.treedata.LenBranchDia',QSM.treedata.LenBranch1Dia', ...
                QSM.treedata.AreBranchDia',QSM.treedata.AreBranch1Dia', ...
                QSM.treedata.VolBranchDia',QSM.treedata.VolBranch1Dia', ...
                QSM.treedata.NumBranchDia',QSM.treedata.NumBranch1Dia', ...
                'VariableNames',{ ...
                'LenBranchDia','LenBranch1Dia', 'AreBranchDia','AreBranch1Dia', ...
                'VolBranchDia','VolBranch1Dia', 'NumBranchDia','NumBranch1Dia'});
            writetable(treedata_cyl_dia_table, str, ...
                'Sheet','treedata_BranchDia');
            
            % export treedata - Height classes, branches
            treedata_cyl_hei_table = table( ...
                QSM.treedata.LenBranchHei',QSM.treedata.LenBranch1Hei', ...
                QSM.treedata.AreBranchHei',QSM.treedata.AreBranch1Hei', ...
                QSM.treedata.VolBranchHei',QSM.treedata.VolBranch1Hei', ...
                QSM.treedata.NumBranchHei',QSM.treedata.NumBranch1Hei', ...
                'VariableNames',{ ...
                'LenBranchHei','LenBranch1Hei', 'AreBranchHei','AreBranch1Hei', ...
                'VolBranchHei','VolBranch1Hei', 'NumBranchHei','NumBranch1Hei'});
            writetable(treedata_cyl_hei_table, str, ...
                'Sheet','treedata_BranchHei');
            
            % export treedata - Zenith classes, branches
            treedata_cyl_zen_table = table( ...
                QSM.treedata.LenBranchZen',QSM.treedata.LenBranch1Zen', ...
                QSM.treedata.AreBranchZen',QSM.treedata.AreBranch1Zen', ...
                QSM.treedata.VolBranchZen',QSM.treedata.VolBranch1Zen', ...
                QSM.treedata.NumBranchZen',QSM.treedata.NumBranch1Zen', ...
                'VariableNames',{ ...
                'LenBranchZen','LenBranch1Zen', 'AreBranchZen','AreBranch1Zen', ...
                'VolBranchZen','VolBranch1Zen', 'NumBranchZen','NumBranch1Zen'});
            writetable(treedata_cyl_zen_table, str, ...
                'Sheet','treedata_BranchZen');
            
            % export treedata - Azimuth classes, branches
            treedata_cyl_azi_table = table( ...
                QSM.treedata.LenBranchAzi',QSM.treedata.LenBranch1Azi', ...
                QSM.treedata.AreBranchAzi',QSM.treedata.AreBranch1Azi', ...
                QSM.treedata.VolBranchAzi',QSM.treedata.VolBranch1Azi', ...
                QSM.treedata.NumBranchAzi',QSM.treedata.NumBranch1Azi', ...
                'VariableNames',{ ...
                'LenBranchAzi','LenBranch1Azi', 'AreBranchAzi','AreBranch1Azi', ...
                'VolBranchAzi','VolBranch1Azi', 'NumBranchAzi','NumBranch1Azi'});
            writetable(treedata_cyl_azi_table, str, ...
                'Sheet','treedata_BranchAzi');
        
        end
    
        % export pmdistance - overview
        pmdist_values = struct2cell(QSM.pmdistance);
        pmdist_names = fieldnames(QSM.pmdistance);
        pmdist_table = cell2table(pmdist_values, ...
            'RowNames', pmdist_names, ...
            'VariableNames',{'Value'});
        writetable(pmdist_table(2:end,:), str, ...
            'WriteRowNames',true, ...
            'Sheet','pmdist_overview');
        
        % export pmdistance data - cylinder distance
        pmdist_cyldist_table = array2table(QSM.pmdistance.CylDist, ...
            'VariableNames',{'CylDist'});
        writetable(pmdist_cyldist_table, str, ...
            'Sheet','pmdist_CylDist');
    
    end

end

%% EXAMPLES

% % using treeqsm output
% % load point cloud
% myTree = 'C:\Daten\Arbeit\Test_TreeQSM\tree.txt';
% P = readtable(myTree); P = table2array(P(:,1:3));
% % compute QSM
% treeqsm(P, inputs)
% % export QSM
% qsm_export("results\QSM_dummy_t1_m1.mat")

% % using select_optimum output
% % load point cloud
% myTree = 'C:\Daten\Arbeit\Test_TreeQSM\tree';
% % compute QSMs
% TreeQSMs = make_models_parallel(myTree, 'multiple_trees', 5, inputs);
% % select best QSM
% select_optimum(TreeQSMs, 'all_mean_dis', 'dummy');
% % save QSM in seperate file
% load("results\OptimalQSMs_dummy.mat",'OptQSM')
% QSM = OptQSM;
% save("results\OptimalQSMs_dummy_best.mat",'QSM')
% % export QSM
% qsm_export("results\OptimalQSMs_dummy_best.mat")
