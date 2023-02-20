#!/bin/bash


#
# bash script_run.sh category pole
#

# vertex : dipole, point
# category : noRegge, Regge
# pole : quark, diquark, both
# kzrange : kzall, kzpos


#^#^#

vertex=$1

category=$2

pole=$3

kzrange=$4


# * * * * * * * * * * * * *


s=${vertex}Vertex_${category}_${pole}Pole_${kzrange}



# # # # # # # # # # # # # #

# # #   FILE PATHS

# # # # # # # # # # # # # # 


#   HEADER PATH   .................................................................


#home_dir=$HOME/Documents/G2code2/pseudo_PDF_2019/model_xzT/pseudo_pdf_final

#home_dir=$HOME/Desktop/Abha/pseudo_PDF_2019/model_xzT/pseudo_pdf_final

home_dir=$HOME/Documents/Research/pseudo_pdf_final

header_path=$home_dir/program_files

export CPLUS_INCLUDE_PATH=$header_path:$CPLUS_INCLUDE_PATH



#   DIRECTORIES   .................................................................


program_dir=$header_path


data_dir=$home_dir/data


plots_dir=$home_dir/plots

kz_plt_dir=$plots_dir/kz_dependence_${vertex}_${kzrange}/${category}_${pole}

pseudo_plt_dir=$plots_dir/pseudo_pdfs_${vertex}_${kzrange}/${category}_${pole}


lc_dir=$program_dir



#   PROGRAM FILES  ................................................................


main_lc=$lc_dir/ioffe_lc.cpp


main_kz=$program_dir/kz_dependence.cpp


main_pseudo=$program_dir/ioffe_dependence.cpp



#   DATA FILES   ..................................................................  


lc_data=$data_dir/ioffe_lc.dat   ##### model on the light cone


kz_data=${data_dir}/${s}_kz.dat


pseudo_data=${data_dir}/${s}_pseudo.dat



#   PLOT SCRIPT FILES   ...........................................................


kz_plt=$plots_dir/kz_dependence_plots.sh


pseudo_plt=$plots_dir/pseudo_plots.sh


#..................................................................................




# # # # # # # # # # # # # # 

# # #   PROGRAM RUNS      

# # # # # # # # # # # # # # 



#   MODEL ON LIGHT CONE   .........................................................


#g++ -Wall $main_lc $header_path/pseudo_pdf_lattice.cpp -lfftw3 -lm -o go

#./go $lc_data



#   KZ DEPENDENCE    ..............................................................


#g++ -Wall $main_kz $header_path/pseudo_pdf_lattice.cpp -lfftw3 -lm -o go_kz

#./go_kz $vertex $category $pole $kzrange $kz_data



#   PSEUDO PDF   ..................................................................


g++ -Wall $main_pseudo $header_path/pseudo_pdf_lattice.cpp -lfftw3 -lm -o go_pseudoPz

#./go_pseudoPz $vertex $category $pole $kzrange $pseudo_data


#..................................................................................




# # # # # # # # # # # # # # 

# # #   PLOTS WITH GNUPLOT     

# # # # # # # # # # # # # # 



#---**********************--------------------------#
#       KZ DEPENDENCE                               #
#---**********************--------------------------#


# columns in dataFile:

# 1 Pz GeV  #2 kz GeV  #3 kz/Pz  #4 fu(kz, Pz) norm 1  #5 fu norm Pz 


#*#-------------------------------------------------------------------------#*#
#*#--- Plot 1a vs kz, fu normalization 1 -----------------------------------#*#
#*#-------------------------------------------------------------------------#*#


filnm="1a_kz_fu_Pz_Norm1_${s}.eps"  

#bash  $kz_plt $kz_data "$kz_plt_dir/$filnm"  


#*#-------------------------------------------------------------------------#*#
#*#--- Plot 1b vs kz, fu normalization Pz ----------------------------------#*#
#*#-------------------------------------------------------------------------#*#


filnm="1b_kz_fu_Pz_NormPz_${s}.eps" 

#bash  $kz_plt $kz_data "$kz_plt_dir/$filnm" 


#*#-------------------------------------------------------------------------#*#
#*#--- Plot 2a vs x, fu normalization 1 ------------------------------------#*#
#*#-------------------------------------------------------------------------#*#


filnm="2a_x_fu_Pz_Norm1_${s}.eps"  

#bash  $kz_plt $kz_data "$kz_plt_dir/$filnm" 


#*#-------------------------------------------------------------------------#*#
#*#--- Plot 2b vs x, fu normalization Pz -----------------------------------#*#
#*#-------------------------------------------------------------------------#*#


filnm="2b_x_fu_Pz_NormPz_${s}.eps"  

#bash  $kz_plt $kz_data "$kz_plt_dir/$filnm" 


#*#-------------------------------------------------------------------------#*#
#*#--- Plot 3a vs x, x times fu normalization 1 ----------------------------#*#
#*#-------------------------------------------------------------------------#*#


filnm="3a_x_xfu_Pz_Norm1_${s}.eps"  

#bash  $kz_plt $kz_data "$kz_plt_dir/$filnm" 


#*#-------------------------------------------------------------------------#*#
#*#--- Plot 3b vs x, x times fu normalization Pz ---------------------------#*#
#*#-------------------------------------------------------------------------#*#


filnm="3b_x_xfu_Pz_NormPz_${s}.eps"  

#bash  $kz_plt $kz_data "$kz_plt_dir/$filnm" 



#*#-------------------------------------------------------------------------#*#


#echo 'kz dependence plots done'


#---**********************--------------------------#
#       PSEUDO PDF                                  #
#---**********************--------------------------#


# columns in dataFile:

#1 nu = P.z  #2 Pz (GeV)  #3 z (1/GeV)  #4 cos (GeV)  #5 sin (GeV)  #6 cos/cos0  #7 sin/cos0  


#*#-----------------------------------------------------------------------------#*#
#*#--- Plot 1a cos vs z for different Pz, P.z (ioffe time) tagged --------------#*#
#*#-----------------------------------------------------------------------------#*#                           


filnm="1a_z_cos_Pz_ioffeTag_${s}.eps"

bash $pseudo_plt $pseudo_data "$pseudo_plt_dir/$filnm" 


#*#-----------------------------------------------------------------------------#*#
#*#--- Plot 1b cos vs z for different P.z (ioffe time), Pz tagged --------------#*#
#*#-----------------------------------------------------------------------------#*#                           


filnm="1b_z_cos_ioffe_PzTag_${s}.eps"

bash  $pseudo_plt $pseudo_data "$pseudo_plt_dir/$filnm"


#*#-----------------------------------------------------------------------------#*#
#*#--- Plot 1c cos ratio vs z for different Pz, P.z (ioffe time) tagged --------#*#
#*#-----------------------------------------------------------------------------#*#                           


filnm="1c_z_cosratio_Pz_ioffeTag_${s}.eps"

bash  $pseudo_plt $pseudo_data "$pseudo_plt_dir/$filnm" 


#*#-----------------------------------------------------------------------------#*#
#*#--- Plot 1d cos ratio vs z for different P.z (ioffe time), Pz tagged --------#*#
#*#-----------------------------------------------------------------------------#*#


filnm="1d_z_cosratio_ioffe_PzTag_${s}.eps"

bash  $pseudo_plt $pseudo_data "$pseudo_plt_dir/$filnm" 


#*#-----------------------------------------------------------------------------#*#
#*#--- Plot 2a cos vs P.z (ioffe time) for different Pz, z tagged --------------#*#
#*#-----------------------------------------------------------------------------#*#                           


filnm="2a_ioffe_cos_Pz_zTag_${s}.eps"

bash  $pseudo_plt $pseudo_data "$pseudo_plt_dir/$filnm" $lc_data


#*#----------------------------------------------------------------------------#*#
#*#--- Plot 2b cos vs P.z (ioffe time) for different z, Pz tagged -------------#*#
#*#----------------------------------------------------------------------------#*#


filnm="2b_ioffe_cos_z_PzTag_${s}.eps"

bash  $pseudo_plt $pseudo_data "$pseudo_plt_dir/$filnm" $lc_data



#*#-----------------------------------------------------------------------------#*#
#*#--- Plot 2c cos ratio vs P.z (ioffe time) for different Pz, z tagged --------#*#
#*#-----------------------------------------------------------------------------#*#                           


filnm="2c_ioffe_cosratio_Pz_zTag_${s}.eps"

bash  $pseudo_plt $pseudo_data "$pseudo_plt_dir/$filnm" $lc_data


#*#-----------------------------------------------------------------------------#*#
#*#--- Plot 2d cos ratio vs P.z (ioffe time) for different z, Pz tagged --------#*#
#*#-----------------------------------------------------------------------------#*#


filnm="2d_ioffe_cosratio_z_PzTag_${s}.eps"

bash  $pseudo_plt $pseudo_data "$pseudo_plt_dir/$filnm" $lc_data


#*#-----------------------------------------------------------------------------#*#
#*#--- Plot 3a cos vs Pz for different z, P.z (ioffe time) tagged --------------#*#
#*#-----------------------------------------------------------------------------#*#


filnm="3a_Pz_cos_z_ioffeTag_${s}.eps"

bash  $pseudo_plt $pseudo_data "$pseudo_plt_dir/$filnm" 


#*#-----------------------------------------------------------------------------#*#
#*#--- Plot 3b cos vs Pz for different P.z (ioffe time), z tagged --------------#*#
#*#-----------------------------------------------------------------------------#*#


filnm="3b_Pz_cos_ioffe_zTag_${s}.eps"

bash  $pseudo_plt $pseudo_data "$pseudo_plt_dir/$filnm" 


#*#-----------------------------------------------------------------------------#*#
#*#--- Plot 3c cos ratio vs Pz for different z, P.z (ioffe time) tagged --------#*#
#*#-----------------------------------------------------------------------------#*#


filnm="3c_Pz_cosratio_z_ioffeTag_${s}.eps"

bash  $pseudo_plt $pseudo_data "$pseudo_plt_dir/$filnm" 


#*#----------------------------------------------------------------------------#*#
#*#--- Plot 3d cos ratio vs Pz for different P.z (ioffe time), z tagged -------#*#
#*#----------------------------------------------------------------------------#*#


filnm="3d_Pz_cosratio_ioffe_zTag_${s}.eps"

bash  $pseudo_plt $pseudo_data "$pseudo_plt_dir/$filnm" 


#*#-----------------------------------------------------------------------------#*#
#*#--- Plot 4a sin vs z for different Pz, P.z (ioffe time) tagged --------------#*#
#*#-----------------------------------------------------------------------------#*#                           


filnm="4a_z_sin_Pz_ioffeTag_${s}.eps"

bash $pseudo_plt $pseudo_data "$pseudo_plt_dir/$filnm" 


#*#-----------------------------------------------------------------------------#*#
#*#--- Plot 4b sin vs z for different P.z (ioffe time), Pz tagged --------------#*#
#*#-----------------------------------------------------------------------------#*#                           


filnm="4b_z_sin_ioffe_PzTag_${s}.eps"

bash  $pseudo_plt $pseudo_data "$pseudo_plt_dir/$filnm"


#*#-----------------------------------------------------------------------------#*#
#*#--- Plot 4c sin ratio vs z for different Pz, P.z (ioffe time) tagged --------#*#
#*#-----------------------------------------------------------------------------#*#                           


filnm="4c_z_sinratio_Pz_ioffeTag_${s}.eps"

bash  $pseudo_plt $pseudo_data "$pseudo_plt_dir/$filnm" 


#*#-----------------------------------------------------------------------------#*#
#*#--- Plot 4d sin ratio vs z for different P.z (ioffe time), Pz tagged --------#*#
#*#-----------------------------------------------------------------------------#*#


filnm="4d_z_sinratio_ioffe_PzTag_${s}.eps"

bash  $pseudo_plt $pseudo_data "$pseudo_plt_dir/$filnm" 


#*#-----------------------------------------------------------------------------#*#
#*#--- Plot 5a sin vs P.z (ioffe time) for different Pz, z tagged --------------#*#
#*#-----------------------------------------------------------------------------#*#                           


filnm="2a_ioffe_sin_Pz_zTag_${s}.eps"

bash  $pseudo_plt $pseudo_data "$pseudo_plt_dir/$filnm" $lc_data


#*#----------------------------------------------------------------------------#*#
#*#--- Plot 5b sin vs P.z (ioffe time) for different z, Pz tagged -------------#*#
#*#----------------------------------------------------------------------------#*#


filnm="5b_ioffe_sin_z_PzTag_${s}.eps"

bash  $pseudo_plt $pseudo_data "$pseudo_plt_dir/$filnm" $lc_data


#*#-----------------------------------------------------------------------------#*#
#*#--- Plot 5c sin ratio vs P.z (ioffe time) for different Pz, z tagged --------#*#
#*#-----------------------------------------------------------------------------#*#                           


filnm="5c_ioffe_sinratio_Pz_zTag_${s}.eps"

bash  $pseudo_plt $pseudo_data "$pseudo_plt_dir/$filnm" $lc_data


#*#-----------------------------------------------------------------------------#*#
#*#--- Plot 5d sin ratio vs P.z (ioffe time) for different z, Pz tagged --------#*#
#*#-----------------------------------------------------------------------------#*#


filnm="5d_ioffe_sinratio_z_PzTag_${s}.eps"

bash  $pseudo_plt $pseudo_data "$pseudo_plt_dir/$filnm" $lc_data


#*#-----------------------------------------------------------------------------#*#
#*#--- Plot 6a sin vs Pz for different z, P.z (ioffe time) tagged --------------#*#
#*#-----------------------------------------------------------------------------#*#


filnm="6a_Pz_sin_z_ioffeTag_${s}.eps"

bash  $pseudo_plt $pseudo_data "$pseudo_plt_dir/$filnm" 


#*#-----------------------------------------------------------------------------#*#
#*#--- Plot 6b sin vs Pz for different P.z (ioffe time), z tagged --------------#*#
#*#-----------------------------------------------------------------------------#*#


filnm="6b_Pz_sin_ioffe_zTag_${s}.eps"

bash  $pseudo_plt $pseudo_data "$pseudo_plt_dir/$filnm" 


#*#-----------------------------------------------------------------------------#*#
#*#--- Plot 6c sin ratio vs Pz for different z, P.z (ioffe time) tagged --------#*#
#*#-----------------------------------------------------------------------------#*#


filnm="6c_Pz_sinratio_z_ioffeTag_${s}.eps"

bash  $pseudo_plt $pseudo_data "$pseudo_plt_dir/$filnm" 


#*#----------------------------------------------------------------------------#*#
#*#--- Plot 6d sin ratio vs Pz for different P.z (ioffe time), z tagged -------#*#
#*#----------------------------------------------------------------------------#*#


filnm="6d_Pz_sinratio_ioffe_zTag_${s}.eps"

bash  $pseudo_plt $pseudo_data "$pseudo_plt_dir/$filnm" 


#*#---------------------------------------------------------------------------#*#


#echo 'pseudo plots done'


#^#^#




