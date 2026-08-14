#!/bin/bash

BASE="/mnt/d/Test_MTLCIRCUIT"
SCR_QUANT_PY="/mnt/d/MTLCIRCUITS/quantify_pathways.py"
timepoints=( "bl" "fu" )
versions=( "1" "2" ) 
hemis=( "left" "right" )


for v in ${versions[@]}; do 
	for t in ${timepoints[@]}; do
		for h in ${hemis[@]}; do
			printf "\nProcessing ${t}_${v} ${h}\n"	
			
			if [ $v -eq 1 ]; then 
				P_WB_TCK="${BASE}/${t}_${v}/Tracts/rWB_PROB_SIFT_2M.tck"
			elif [ $v -eq 2 ]; then 
				P_WB_TCK="${BASE}/${t}_${v}/Tracts/WB_ACT_SIFT_2M.tck"
		     	fi
					
			P_NODE="${BASE}/${t}_1/subfields/NODES_${h}_${t}_1.nii.gz"
			P_NODE_MIF="${P_NODE%.nii.gz}.mif"								
			if [ ! -f ${P_NODE} ]; then 
				printf "\nMissing NODE Files. Starting...\n"
				python3 << PY 
						
import nibabel as nib
import numpy as np
import os 

print("\n-----# Node Creator #-----\n") 

lbl_ashs = {"ERC": 10,
	    "SUB": 8,
	    "CA1": 1,
            "CA3":  4, 
	    "DG":  3
	    }	 

lbl_node = {1 : 1,
	    4 : 2, 
	    3 : 3, 
            8 : 4,
	    10 : 5
	    }		

print(f"Loading ${h} ashs...") 

p_ashs = os.path.join("${BASE}", f"${t}_1", "coregistration", f"rASHS_${h}.nii.gz") 
i_ashs = nib.load(p_ashs)
aff_ashs = i_ashs.affine
arr_ashs = i_ashs.get_fdata()

print("Creating Nodes...") 

node_space = np.zeros_like(arr_ashs)
for r,l_orig in lbl_ashs.items(): 
	p_subf_dir = os.path.join("${BASE}", "${t}_${v}", "subfields")
	node = nib.load(os.path.join(p_subf_dir, f"{r}_${h}_${t}_1.nii.gz")).get_fdata()
	l_new = lbl_node[l_orig]	
	node_space[node == l_orig] = l_new 
	
	print(f"Added node {r}...")
	
print(f"\nSaving nodes ${t}_${v}...")

i_nodes = nib.Nifti1Image(node_space, aff_ashs) 
nib.save(i_nodes, "${P_NODE}")

print("Saved file ${P_NODE}\n")


print("\n-----# Created all Nodes. #-----\n")

PY

			else 
				echo Nodes already exist. Proceeding..
			fi
		
			if [ ! -f $P_NODE_MIF ]; then 
				echo Nodes only in nifti format. Converting to mif...
				mrconvert $P_NODE $P_NODE_MIF
			fi
			
			
			P_OUT_DIR="${BASE}/${t}_${v}/Tracts/MTL_tracts"
			CSV_CONN="${P_OUT_DIR}/connectome.csv"
			ASSIGNMENTS="${P_OUT_DIR}/assignments.txt"
			P_MSP="${P_OUT_DIR}/MSP_ERC_CA1.tck"
			P_TSP_1="${P_OUT_DIR}/TSP_DG_CA3.tck"
			P_TSP_2="${P_OUT_DIR}/TSP_CA3_CA1.tck"

			mkdir -p $P_OUT_DIR
			
			printf "\nCreating Connectomes..."
		      	
			if [ ! -f $CSV_CONN ];	then		
				tck2connectome \
					$P_WB_TCK \
					$P_NODE_MIF \
					$CSV_CONN \
					-out_assignments $ASSIGNMENTS \
					-assignment_radial_search 1 \
					-symmetric \
					-zero_diagonal \
					-force
						
			else
				printf "Connectome .csv already exists. Skipping..."
			fi

			if [ ! -f $P_MSP ]; then
				echo Extracting MSP connectome...
				
				connectome2tck \
					$P_WB_TCK \
					$ASSIGNMENTS \
					$P_MSP \
					-nodes 5,1 \
					-exclusive \
					-files single \
					-force
		
			else
				echo MPS connectome already exists. Skipping... 
			fi
					
			if [ ! -f $P_TSP_1 ]; then
				echo Extracting TSP_1 cnnectome...

				connectome2tck \
					$P_WB_TCK \
					$ASSIGNMENTS \
					$P_TSP_1 \
					-nodes 5,3,2 \
					-exclusive \
					-files single \
					-force

			else
				echo TSP_1 connectome already exists. Skipping...
			fi
			
			if [ ! -f $P_TSP_2 ]; then
				echo Extracting TSP_2 cnnectome...

				connectome2tck \
					$P_WB_TCK \
					$ASSIGNMENTS \
					$P_TSP_2 \
					-nodes 2,1 \
					-exclusive \
					-files single \
					-force

			else
				echo TSP_2 connectome already exists. Skipping...
			fi
			
			P_QUANT="${P_OUT_DIR}/${h}_${t}_${v}_pathway_stats.csv"
			python3	${SCR_QUANT_PY} ${CSV_CONN} ${t} ${v} ${h} ${P_QUANT}



		done
	done
done		



	
		
