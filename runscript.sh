#!/bin/bash


while test $# -gt 0; do
           case "$1" in
                -name)
                    shift
                    name=$1
                    shift
                    ;;
                -input)
                    shift
                    input=$1
                    shift
                    ;;
				-match)
                    shift
                    match=$1
                    shift
                    ;;
               
			     -gap)
                    shift
                    gap=$1
                    shift
                    ;;	
				-missmatch)
                    shift
                    missmatch=$1
                    shift
                    ;;
				-threads)
                    shift
                    threads=$1
                    shift
                    ;;
                *)
                   echo "$1 is not a recognized flag!"
                   return 1;
                   ;;
		  esac

		 
  done 

generatedFilesDir=${input}
srcDir=${input}

for program in "seq" "omp_fine" "omp_coarse" "pthreads_fine" "pthreads_coarse"
do	
	if [ $program = "seq" ]; then 
		gcc   "${srcDir}"${program}.c -o "${generatedFilesDir}"${program}
        elif [ $program = "omp_fine" -o $program = "omp_coarse" ]; then 
		gcc   "${srcDir}"${program}.c -o "${generatedFilesDir}"${program} -fopenmp
	else
		gcc   "${srcDir}"${program}.c -o "${generatedFilesDir}"${program} -pthread 
	fi 
done

for program in "seq" "omp_coarse" "omp_fine" "pthreads_coarse" "pthreads_fine"
do
    for doc in "D1.txt" "D2.txt" "D3.txt" "D4.txt" "D5.txt" "D9.txt"
    do
		if  [ $program = "seq" ] && [ $doc="D1.txt" -o $doc="D2.txt" -o $doc="D3.txt" -o $doc="D4.txt" -o $doc="D5.txt" -o $doc="D9.txt" ]; then #
                       report=${input}"Report_"${name}".txt"
		       path=${input}${doc}
		    echo   
                    echo --------- ${program} $report $path $match $missmatch $gap ---------------
		        ./${program} $report $path $match $missmatch $gap 
		elif [ -z "$threads" ]; then
			         echo "Not enough arguments"

                elif [ $program = "omp_fine" ] && [ $doc="D1.txt" -o $doc="D2.txt" -o $doc="D3.txt" -o $doc="D4.txt" -o $doc="D5.txt" -o $doc="D9.txt" ]; then #
                     report=${input}"Report_"${name}"_OMP_"$threads".txt"
			         path=${input}${doc}
			echo
                       echo --------- ${program} $report $path $match $missmatch $gap $threads --------- 
                      ./${program} $report $path $match $missmatch $gap $threads 
                elif [ $program = "omp_coarse" ] && [ $doc="D9.txt" ]; then 
                     report=${input}"Report_"${name}"_OMP_"$threads".txt"
		     path=${input}${doc}
		     echo
                       echo --------- ${program} $report $path $match $missmatch $gap $threads ---------  
          	        ./${program} $report $path $match $missmatch $gap $threads 
		elif [ $program = "pthreads_fine" ] && [ $doc="D1.txt" -o $doc="D2.txt" -o $doc="D3.txt" -o $doc="D4.txt" -o $doc="D5.txt" -o $doc="D9.txt" ]; then #
			 report=${input}"Report_"${name}"_PTH_"$threads".txt"
			 path=${input}${doc}
			echo 
                        echo --------- ${program} $report $path $match $missmatch $gap $threads ---------  
			./${program} $report $path $match $missmatch $gap $threads	
		elif [ $program = "pthreads_coarse" ] && [ $doc="D9.txt" ]; then 
		          report=${input}"Report_"${name}"_PTH_"$threads".txt"
			  path=${input}${doc}
                          echo
			  echo --------- ${program} $report $path $match $missmatch $gap $threads --------- 
			   ./${program} $report $path $match $missmatch $gap $threads           		
		fi
    done
done
