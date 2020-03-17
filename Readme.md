# CRF Cell ID

<img src = "extra/readme_img_v2.jpg" width=45% align="right">

This software/tools were developed for automatically detecting cells in dense brain images, annotating biological identities, visualizing results, and building new atlases using ground-truth annotations. In contrast to popular registration based methods, `CRF_Cell_ID` formulates a structured prediction problem that is solved with second-order Conditional Random Fields model. Thus `CRF_Cell_ID` finds optimal labels for all cells that are most consistent with prior knowledge in terms of maintaining second-order features such as positional relationships. `CRF_Cell_ID` provides a computationally scalable method for building atlases from annotated data.

To learn about details of the method, please read the <a href="https://www.biorxiv.org/content/10.1101/2020.03.10.986356v1">paper</a>. A detailed guide on using the framework is __coming soon__

This repository contains - 
1. Code for predicting identities in new datasets
2. Code for building new atlas based on ground-truth annotated datasets
3. Raw datasets used in the paper, ground-truth annotations, and codes for reproducing results

<img src = "extra/video1_gif.gif" width=100% align="center">

# Installation
Please make sure to keep the name of the directory as `CRF_Cell_ID` so that walkthrough instructions can be run without changing any paths.
### Installing with Git Shell
```
git clone https://github.com/shiveshc/CRF_Cell_ID.git
```
### Installing online
Click on `clone or download` button on top right corner and `Download ZIP`

### Dependencies
1. MATLAB 2018b
2. <a href = "https://www.cs.ubc.ca/~schmidtm/Software/UGM.html">UGM Toolbox</a>
3. <a href = "https://sites.google.com/site/myronenko/research/cpd">CPD</a> registration (optional, required only for comparisons)
4. <a href = "https://github.com/Vaa3D">Vaa3D</a> (optional, required only for visualizing ground-truth annotations)

These dependencies except MATLAB and Vaa3D are provided in the repository so no need to download them individually. However please consider citing them.

### Description of respository contents
1. __Datasets__ : contains all raw datasets used in the <a href="https://www.biorxiv.org/content/10.1101/2020.03.10.986356v1">paper</a> and their ground-truth annotations.
   - __ComplexOrientations__ : contains 7 datasets of animals imaged with varying rotations along AP axis. `marker_names.xlsx` files in  folders record the names and IDs of manually annotated cells. `markers` file is a Vaa3D marker annotation file that records IDs and positions of all cells. Folder with suffix `_B`, `_C` or `_mN` denote BFP, CyOFP and mNeptune channels
   - __GeneExpressionAnalysis__ : contains 21 datasets of animals used in gene-expression analysis. `.tif` folders contain raw z-planes of 3D stacks acquired in two channels: pan-neuronal RFP and GFP
   - __MultiCellCalciumImaging__ : contains 31 datasets of animals used in multi-cell calcium imaging analysis. `.tif` folders contain raw z-planes of 3D stacks acquired in one channel: GFP
   - __NeuroPAL__ : contains 9 datasets of NeuroPAL animals.`marker_names.xlsx` files in these folders record the names and IDs of manually annotated cells. `markers` file is a Vaa3D marker annotation file that records IDs and positions of all cells in image stack.
2. __UGM__ : contains the undirected graphical models library used to formulate and optimize CRF model
3. __utils__ : contains auxillary functions for -
   - reading images and marker files : `read_marker_files`, `readCyOFP_SD`, `readtif_Brucker_for_preprocess_data`
   - detecting cells : `em_segmentation_3d_gpu`, `segment_red_first_half`, `segment_red_second_half`, `remove_close_neurons`
   - preprocessing data : `preprocess_alternate_strain_data`, `preprocess_alternate_strain_data_brucker`, `define_axes_specifying_neurons`
   - registration based prediction : `registration_based_matching_multi`, `registration_based_matching_marker_only`, `registration_based_matching_marker_only_missing_neurons`
   - visualizing intermediate steps and results : `visualize_label_results_marker_only`, `visualize_multi_label_results`, `plot_all_labels_on_img`, `visualize_duplicates_with_scores`, `plot_node_pot`, `plot_edge_pot`
4. __utils_NeuroPAL__ : contains auxillary functions for -
   - reading marker files : `read_marker_files`, `read_marker_files_wo_marker_name`
   - preprocessing data : `preprocess_NeuroPAL_data`, `defining_axes_specifying_neurons`
   - registration based prediction : `registration_based_matching_multi`, `registration_based_matching_2D`
   - comparing absolute position vs relative position variability : `absolute_positions_of_annotated_neurons`, `relative_position_consistency_of_annotated_neurons`
   - building data-driven atlas : `compile_annotated_data`, `color_intensity_distributions`, `update_atlas_based_on_annotation`, `update_atlas_based_on_annotation_with_rotated`
   - visualizing data and prediction : `rgb_img_all_files_in_folder`, `rgb_img_from_indev_channels`, `rgb_img_annotated_data`, `plot_labels_on_rgb_img`



	
	
	
  

  
