# CRF Cell ID

<img src = "extra/readme_img_v2.jpg" width=45% align="right">

This software/tools were developed for automatically detecting cells in dense brain images, annotating biological identities, visualizing results, and building new atlases using ground-truth annotations. In contrast to popular registration based methods, `CRF_Cell_ID` formulates a structured prediction problem that is solved with second-order Conditional Random Fields model. Thus `CRF_Cell_ID` finds optimal labels for all cells that are most consistent with prior knowledge in terms of maintaining second-order features such as positional relationships. `CRF_Cell_ID` also provides a computationally scalable method for building atlases from thousands of annotated images.

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
Click on `Clone or download` button on top right corner and `Download ZIP`

### Dependencies
1. MATLAB 2018b
2. <a href = "https://www.cs.ubc.ca/~schmidtm/Software/UGM.html">UGM Toolbox</a>
3. <a href = "https://sites.google.com/site/myronenko/research/cpd">CPD</a> registration (optional, required only for comparisons)
4. <a href = "https://github.com/Vaa3D">Vaa3D</a> (optional, required only for visualizing ground-truth annotations)

UGM and CPD are provided in the repository so there is no need to download them separately. However please consider citing them.

# Description of respository contents
1. __Main__ : This folder contains the main code for predicting cell identities in new datasets. Includes codes for - 
   - reading image data : `read_tif`
   - preparing data for annotation : `preprocess_data`
   - predicting identities : `annotation_CRF_landmark`
2. __Datasets__ : contains all raw datasets used in the <a href="https://www.biorxiv.org/content/10.1101/2020.03.10.986356v1">paper</a> and their ground-truth annotations.
   - __ComplexOrientations__ : contains 7 datasets of animals imaged with varying rotations along AP axis. `marker_names.xlsx` files in  folders record the names and IDs of manually annotated cells. `markers` file is a Vaa3D marker annotation file that records IDs and positions of all cells. Folder with suffix `_B`, `_C` or `_mN` denote BFP, CyOFP and mNeptune channels
   - __GeneExpressionAnalysis__ : contains 21 datasets of animals used in gene-expression analysis. `.tif` folders contain raw z-planes of 3D stacks acquired in two channels: pan-neuronal RFP and GFP
   - __MultiCellCalciumImaging__ : contains 31 datasets of animals used in multi-cell calcium imaging analysis. `.tif` folders contain raw z-planes of 3D stacks acquired in one channel: GFP
   - __NeuroPAL__ : contains 9 datasets of NeuroPAL animals.`marker_names.xlsx` files in these folders record the names and IDs of manually annotated cells. `markers` file is a Vaa3D marker annotation file that records IDs and positions of all cells in image stack.
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
5. __Run_GeneExpression__ : contains codes for predicting cell identities in `GeneExpressionAnalysis` dataset, reproducing prediction comparison results for different methods.
6. __Run_MultiCellCalciumImaging__ : contains codes for predicting cell identities in `MultiCellCalciumImaging` dataset, reproducing prediction comparison results for different methods.
7. __Run_NeuroPAL__ : ccontains codes for predicting cell identities in `NeuroPAL` dataset, reproducing prediction comparison results for different methods. These methods include `Color_Registration`, `CRF`, `CRF_Color`, `CRF_Registration`, and `CRF_Color_Registration`. Additionally registration based matching can be performed using `registration_based_matching_multi` in `utils_NeuroPAL` folder.
8. __Run_SyntheticData__ : contains code for predicting identities in synthetic datasets, analyze the effects of position noise and count noise on accuracy, tuning model, simulating landmark choice and comparison with other methods on synthetic datasets
9. __UGM__ : contains the undirected graphical models library used to formulate and optimize CRF model

# Step-by-Step Walkthrough
Here we provide a detailed step-by-step walkthrough of the framework.

First change current working directory in MATLAB to the cloned `CRF_Cell_ID` folder and add `Main` to path
```
cd '\CRF_Cell_ID' % provide full path of CRF_Cell_ID
addpath('Main')
```

### 1. Read image data
The main function for reading image data is `read_tif` in `Main` folder.
```
Inputs:
 style		'one' if the 3D image stack to be read is saved as one .tif file or 
 		'separate' if individual z-planes are saved as separate images.
 input		image path if style is 'one' e.g. 'CRF_Cell_ID\test_data\img1.tif' or
 		image directory path if style is 'separate' e.g. 'CRF_Cell_ID\test_data\img2'
		
Outputs:
 img		3d image stack variable
```
e.g. we'll read `sample_data1` in `sample_run` folder. Two different channels are present in `sample_data1` as `img1` and `img2`. Each folder contains z-planes saved separately.

```
image1 = read_tif('separate','sample_run\sample_data1\img1');
image2 = read_tif('separate','sample_run\sample_data1\img2');
```
visualizing the max-projection of the two channels, output should look like below
```
figure,imshow(max(mat2gray(image1),[],3),'border','tight')
caxis([0,0.3])
figure,imshow(max(mat2gray(image2),[],3),'border','tight')
caxis([0,0.3])
```
<img src = "extra/sample_data1_img1.jpg" width=50% ><img src = "extra/sample_data1_img2.jpg" width=50% >

We'll predict identities of cells in the channel on the right and will use the channel on the left as landmark channel i.e. channel in which identities of cells are known. In practice `CRF_Cell_ID` supports any number of landmarks channels as will be shown below. Further, it is not necessary to know the identities of all cells in a landmark channel. Thus `CRF_Cell_ID` easily integrates available information from multiple channels.

### 2. Preprocess data for annotation
In this step we'll 1) detect all cells in each channel read in previous step. 2) specify landmark information 3) generate a coordinate axes in head, and 4) generate a data file that will be used as input in the next step. 

The main script to do this is `preprocess_data` in `Main`
```
Inputs :
 out_folder		directory path where the output data file be saved
 data_name		name of the output data file
 img1			image channel in which identities of cells are to be predicted
 img1_marker		Optional Vaa3d marker file if manually detected cells in img1. If not available set as []
 img1_marker_name	Optional .xlsx file storing the names of manually detected cells in img1. If not available set as []
 varargin		variable number of image channels that can be used as landmark channels. Each varargin is triplet combination of 
 			image, marker file and marker name file e.g. {img,\img_marker,\img_marker_name.xlsx} or {img,[],[]} if manual
			detection is not available.

Outputs :
 data file		.mat data file that will be used as input in next step
```

For us, `image2` is the image channel in which cell identities are to be predicted thus our first argument. Further we are not going to use any manual detection, there we run
```
preprocess_data('sample_run\sample_data1',...
	'data_annotation_sample_data1',...
	image2,...
	[],...
	[],...
	{image1,[],[]})
```
#### 1. Detecting cells in image channels
The output on terminal should look like
```
Segmenting main image channel ... 
Enter thresh_parameter value between 0-100. Higher values for detecting few but brightest cells - 
```

Here the `thresh_parameter` is an input parameter for segmenting cells which is used to initialize the means of gaussian mixture components. We'll start with trying 99.95 as the value. This will generate an output image showing detected cells

<img src = "extra/thresh_param_1.jpg" width=70% >

Clearly the parameeter value is too high as many cells are not detected thus we'd like to specify the parameter value again. On the terminal you shoud see the message asking to accept/reject the current parameter value

```
Enter thresh_parameter value between 0-100. Higher values for detecting few but brightest cells - 99.95
Enter 'y' if accept thresh_parameter else enter 'n' -
```

We'll enter 'n' becuase we are not satisfied with the chosen `thresh_parameter` and lower the value to detect more cells. After a couple of trial and errors we select 99.83 as the parameter value. In this case, the image should look like

<img src = "extra/thresh_param_2.jpg" width=70% >

Since we're satisfied with this parameter, we'll enter 'y' on terminal this time. Now the GM model will be fitted to image data to segment cells. The terminal should look like
```
Enter 'y' if accept thresh_parameter else enter 'n' -'y'

ans =

      228643           3


ans =

    'iter:1, likelihood:-12.5171, error:0.34564'


ans =

    'iter:2, likelihood:-12.5379, error:0.32442'


ans =

    'iter:3, likelihood:-12.5495, error:0.31976'


ans =

    'iter:4, likelihood:-12.5559, error:0.3181'


ans =

    'iter:5, likelihood:-12.5598, error:0.31733'
```
Next the same process will be repeated for all landmark channels. In this case, the output should look like
```
Segmenting image channel 2
Enter thresh_parameter value between 0-100. Higher values for detecting few but brightest cells - 99.98
Enter 'y' if accept thresh_parameter else enter 'n' -'n'
Enter thresh_parameter value between 0-100. Higher values for detecting few but brightest cells - 99.97
```
<img src = "extra/thresh_param_3.jpg" width=70% >

```
Enter 'y' if accept thresh_parameter else enter 'n' -'y'

ans =

       36865           3


ans =

    'iter:1, likelihood:-10.4939, error:0.37325'


ans =

    'iter:2, likelihood:-10.4993, error:0.36659'


ans =

    'iter:3, likelihood:-10.5018, error:0.36522'


ans =

    'iter:4, likelihood:-10.5028, error:0.36471'
```

_Skipping cell detection in specific channels (optional)_

Note that automatic cell segmentation step can be skipped for any image channel if automatic segmentation does not generate good results and manual detections for the channel are available such as marker annotation files generated with Vaa3D. e.g. instead of automatic detetction for `image1` we can use `img1_markers` and `img1_marker_names.xlsx` in `sample_data_1`. In this case, we run `precprocess_data` as follows

```
preprocess_data('sample_run\sample_data1',...
	'data_annotation_sample_data1',...
	image2,...
	[],...
	[],...
	{image1,'sample_run\sample_data1\img1_markers','sample_run\sample_data1\img1_marker_names.xlsx'})
```
In the above, automatic segmentation will be performed only for `image2` and not for `image1`.

#### 2. Specifying landmarks in image channels
Next we'll specify landmark information in channels by going through each landmark channel one-by-one, specifying landmark cell in each channel (via user input) and entering the selected cell's name. The prompt at terminal will ask for which image channels should be used for specifying landmarks

```
Enter which channels to use for specifying landmarks e.g [2,4] else enter blank (single quotes) -
```

We'll use `image1` channel to specify identities of easily identified cells in whole-brain stack as well as the landmark channel i.e `image2`, thus, enter [1,2]. Note, here 1 always denotes the `img1` argument and channel 2 onwards denote images provided in `varargin`. If you do not want to specify landmarks in any channel, enter `[]`

```
Enter which channels to use for specifying landmarks e.g [2,4] else enter blank (single quotes) -[1,2]
```

In this case we'll first see `image2` for specifying landmarks. With the cursors on the image, users can click on any cell whose identity they want to specify. We also get the following prompt on terminal

<img src = "extra/landmark_1.jpg" width=70% > 

```
Enter name of the selected landmark e.g. 'RMEL' -
```

After clicking on the cell on the image e.g. the one highlighted in image above we enter its name on the terminal 'RMEV'.

```
Enter name of the selected landmark e.g. 'RMEL' -'RMEV'
```

Next we get the follwoing prompt
```
If done with this channel, enter 'y' -
```
We'll enter 'n' since we want to specify more landmark names. When done, we'll specify 'y'. Thus, following the promts on terminal, users can easily specify identities of as many landmarks as they want e.g.

<img src = "extra/landmark_2.jpg" width=70% >

```
Enter name of the selected landmark e.g. 'RMEL' -'RMED'
If done with this channel, enter 'y' -'y'
```

Next the same process will be performed one-by-one for all channels specified above, `[1,2]` in our case. Thus we specify landmarks in `image1` next.

<img src = "extra/landmark_3.jpg" width=70% >

```
Enter name of the selected landmark e.g. 'RMEL' -'ASGL'
If done with this channel, enter 'y' -'n'

```

<img src = "extra/landmark_4.jpg" width=70% >

```
Enter name of the selected landmark e.g. 'RMEL' -'AIBL'
If done with this channel, enter 'y' -'y'
```

#### 3. Define axes specifying cells
Next we'll specify axes specifying cells via user input. These cells enable defining a consistent coordinate system in head which is crucical for accurately extracting features from image. To do so, the terminal prompt asks to sequentially click on cells in the image.

First we will click on two cells in the image, one in the anterior and one in the posterior region of head. These cells will be used to define anterior-posterior (AP) axis. Here, users can click on any cell; the specific identity of selected cell doesn't matter as these cells are used only to check the consistency of direction of axis. First, we click on a cell in anterior of head ganglion.
```
ans =

    'select (A)-P neuron. If not present than right click.'
```
<img src = "extra/A_neuron.jpg" width=70% >

Next, we click on any cell in the posterior of head. Again, any cell can be clicked on as long as it is posterior to the cell selected in anterior region in previous step.
```
ans =

    'select A-(P) neuron. If not present than right click.'
```
<img src = "extra/P_neuron.jpg" width=70% >

Next, we will click on two cells that will help in defining left-right (LR) axis. Since distinguishing between cells in the left and right region of head in current view of the image is difficult, we'll skip selecting these by right-clicking anywhere in the image
```
ans =

    'select (L)-R neuron. If not present than right click.'


ans =

    'select L-(R) neuron. If not present than right click.'
```
Lastly, we will click on two cells that will help in defining dorsal-ventral (DV) axis.
```
ans =

    'select (D)-V neuron. If not present than right click.'
```
<img src = "extra/D_neuron.jpg" width=70% >
```
ans =

    'select D-(V) neuron. If not present than right click.'
```
<img src = "extra/V_neuron.jpg" width=70% >





	
	
	
  

  
