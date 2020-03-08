%%%% function to correct naming of multi-simulation on real data output
%%%% files. Can be run after each PACE run to rename results so that
%%%% remaining simulations can be run again.

function correct_multi_simulation_ouput_names(inp_dir)
%     inp_dir = ['C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\PACE\RealData\20190414_AML70_unc47_gcy32_CyOFP\Results_multi_data',num2str(i),'_numLand3'];
%     inp_dir = ['C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\PACE\AtlasTest\Results\Results_130rand_206label_randLandmarks\Results_numLabelRemove_50_LB'];
    fileList = dir(inp_dir);
    rng_fix_parameter = [];
    sim_number = [];
    for i = 1:size(fileList,1)
        name = fileList(i).name;
        if ~strcmp(name,'.') && ~strcmp(name,'..')
            underscore_pos = strfind(name,'_');
            rng_fix_parameter = [rng_fix_parameter;str2num(name(1,underscore_pos(1,end-1)+1:underscore_pos(1,end)-1)),i];
            sim_number = [sim_number;str2num(name(1,underscore_pos(1,end)+1:size(name,2) - 4))];
        end
    end
    if ~isempty(rng_fix_parameter)
        uniq_rng_fix_parameter = unique(rng_fix_parameter(:,1));

        for i = 1:size(uniq_rng_fix_parameter,1)
            cnt = 0;
            curr_rng_fix_parameter_index = find(uniq_rng_fix_parameter(i,1) == rng_fix_parameter(:,1));
            curr_rng_fix_parameter_sim = sim_number(curr_rng_fix_parameter_index,:);
            [sort_sim,sort_index] = sort(curr_rng_fix_parameter_sim);
            for j = 1:size(sort_index,1)
                name = fileList(rng_fix_parameter(curr_rng_fix_parameter_index(sort_index(j,1),1),2)).name;
                underscore_pos = strfind(name,'_');
                old_name = [inp_dir,'\',name];
                new_name = [inp_dir,'\',name(1,1:underscore_pos(1,end)),num2str(cnt+1),'.mat'];
                if ~strcmp(old_name,new_name)
                    movefile(old_name,new_name)
                end
                cnt = cnt + 1;
            end
        end
    end
end