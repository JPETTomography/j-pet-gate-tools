README for root2lst_multi
(root2lst_multi_old assumes incorrectly there is a block at the top and rotates all LORs by 11.25 degrees).
Written by Geron Bindseil July 2011

This program processes Gate root files containing coincidences and outputs a Siemens Inveon list mode data format ".lst" file that can be used by Inveon Acquisition Workplace for making sinograms and reconstruction.

The directory where root files are stored should contain files named with the following convention: filename_0.root, filename_1.root, filename_2.root, etc.

You must specify the number of files to process in the directory.

From the command line, the usage for root2lst_multi is as follows:


Usage: ./root2lst_multi <trues> <num_root_files> <rootfile_path> <output_file_name>

<trues> set to 'trues' if you want only trues and 'prompts' if you want all prompts and 'prompts_plus_delays' if you want to include random coincidences as delayed event packets

<num_root_files> set the number of root files you have in the directory to process

<first_root_file_name> set the path to the first root file. For example: /Data/inveonPET/inveonPET_0.root

<output_file_name> set the output file name. For example: '/Data/inveonPET/filename.lst'
