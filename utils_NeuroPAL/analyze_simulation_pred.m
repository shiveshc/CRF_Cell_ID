%%%% function to analyze neuroPal run results

function pred_struct = analyze_simulation_pred()
%     addpath(genpath('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\'))
    in_direc = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos\Results_numLand_NotMultiLabel\Results_3D_data_20190706_5_numLand15\',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos\Results_numLand_NotMultiLabel\Results_3D_data_20190706_21_numLand15\',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos\Results_numLand_NotMultiLabel\Results_3D_data_20190710_2_numLand15\',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos\Results_numLand_NotMultiLabel\Results_3D_data_20190710_5_numLand15\',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos\Results_numLand_NotMultiLabel\Results_3D_data_20190710_8_numLand15\',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos\Results_numLand_NotMultiLabel\Results_3D_data_20190710_10_numLand15\',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos\Results_numLand_NotMultiLabel\Results_3D_data_20190710_22_numLand15\',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos\Results_numLand_NotMultiLabel\Results_3D_data_20190710_24_numLand15\',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos\Results_numLand_NotMultiLabel\Results_3D_data_20190710_27_numLand15\'};
    
    pred_struct = struct();
    
    
    for f = 1:size(in_direc,2)
        correct_fraction = [];
        fileList = dir(in_direc{1,f});
        for i = 1:size(fileList,1)
            name = fileList(i).name;
            if ~strcmp(name,'.') && ~strcmp(name,'..')
                load([in_direc{1,f},'\',name])
                correct_fraction = [correct_fraction;experiments(1).landmark_match_score(:,end)];
            end
        end
        pred_struct(f).data_name = data{1,f};
        pred_struct(f).correct_frac = correct_fraction;
    end

    copy_data = [];
    for i = 1:size(pred_struct,2)
        copy_data = [copy_data;cat(2,repmat(i,size(pred_struct(i).correct_frac,1),1),pred_struct(i).correct_frac)];
    end
end