#!/bin/bash


#
# bash script_run_loop.sh


#^#^#


home_dir=$HOME/Documents/Research/pseudo_pdf_final

script=${home_dir}/script_run.sh


# * * * * * * * 

# set doLoop y to run loop

doLoop=n            


# * * * * * * * 

# single case variables :

category=Regge
pole=diquark


# * * * * * * * 


declare -a category_ar=(noRegge Regge)
declare -a pole_ar=(quark diquark both)


if [ $doLoop == y ]
then
    for(( i=0 ; i < 2 ; i++))
    do
	for(( j=0 ; j < 3 ; j++))
	do
	    category=${category_ar[i]}
             pole=${pole_ar[j]}

	     bash $script $category $pole
		
	done
    done
else

    bash $script $category $pole 	

fi



#^#^#




