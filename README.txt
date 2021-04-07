
* NOTE: These steps up to the last have already been done (right up to and including running 
	VPAdataSelection.m or VSdataSelection.m for all EEG data files). 
	I can't include all of the EEG data (.edf files) as they exceed the 100MB submission limit.
	Instead I have included a sample file in the "Sample EDF File" folder.
	However all of the exported data from this stage has been stored in "ERP Files" in their
	respective task folders, under their respective groupings and electrode type.

* NOTE: I have included all steps of the pipeline for completeness, however, only the
	Python files need to be run to classify the ERP data. The entire pipeline can
	be run with the sample file at the discretion of the user.

* IMPORTANT: The only files which need to be run to observe the results discussed in
	     this project are VPAClassifier.py and VSClassifier.py as they are (without changes).

-------------------------------------------------------------------------------------------------------------------------------------------

To run EEG data (.edf) through the entire pipeline:

- EEGLab must be downloaded from https://sccn.ucsd.edu/eeglab/download.php to complete the signal processing
  and metric extraction stage of the pipeline.
- Open MATLAB and make sure the EEGLab folder is in MATLAB's path (you must add this to the MATLAB path variables if not).
- Run EEGLab in MATLAB by typing "eeglab" into the command window.
- Import your wanted files:
	+ Go to File -> Import data -> Using the BIOSIG Interface, 
	  then navigate to the .edf file you want to load (sample .edf file included in "Sample EDF File").
	  do not fill out any of the fields and click OK.
	+ Go to File -> Import event info -> from E-Prime ASCII (text) file.
	  Then select the corresponding event info text file in the same folder as the 
	  data you've just imported.
	  Fill in the fields as follows:
	  	- Input field (column) names: type "latency type" without the quotes.
          	- Number of file header lines: type the value 1.
	  	- Time unit (sec): type the value 1.
	  	- Align event latencies to data events: type "NaN" without the quotes.
	  	- Click OK to import the event info.
- Then run VPAdataSelection.m or VSdataSelection.m depending on the .edf file loaded in.
  (VPAdataSelection.m is to be run on the sample file included.)
  Export to whatever folder you want by changing destination as indicated in the comment
  at the top of the file.
- Then run the relevant Python script, VPAClassifier.py or VSClassifier.py on the output
  from the previous step.

-------------------------------------------------------------------------------------------------------------------------------------------