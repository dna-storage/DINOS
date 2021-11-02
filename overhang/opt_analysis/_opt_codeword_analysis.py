'''
Filename: _opt_codeword_analysis.py

Description: Method defintions that perform analysis to illustrate an optimal codeword size - overhang count

Author: Kevin Volkel

'''

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import overhang.tree as tree
import overhang.reaction_node as node
import logging
from overhang.dnastorage_utils.system.dnafile import *
import os
import shutil
import math
import numpy as np
from scipy import stats

import overhang.plot_utils.plot_utils as plt_util
import time
import pickle as pi

tlogger=logging.getLogger('dna.overhang.tools.tree_analysis')
tlogger.addHandler(logging.NullHandler())

def _sweep_overhangs_codewordsize(self,data_buffer,workloadID,debug_fname=None): #workloadID[0] == category, workloadID[1] == filename
    output_filename=debug_fname
    overhang_list=[3,5,9,17,65]
    codeword_list=[1,2,3,4,5,8]
    overhang_index_mapping={} #map number of overhangs to an index (used to index into result array)
    codewordsize_index_mapping={}
    strand_length_bytes=509 #508 bytes to account for the extra 4 bytes of indexing
    index_bytes=3
    root_prefix=os.path.normpath(self._out_dir[workloadID[0]]["root"])+'/'
    
    #set up result dictionary for the file we are analyzing 
    if workloadID[0] not in self._opt_codeword_results:
        self._opt_codeword_results[workloadID[0]]={}
        if workloadID[1] not in self._opt_codeword_results[workloadID[0]]:
            self._opt_codeword_results[workloadID[0]][workloadID[1]]={}
    else:
        if workloadID[1] not in self._opt_codeword_results[workloadID[0]]:
            self._opt_codeword_results[workloadID[0]][workloadID[1]]={}

    
    self._opt_codeword_results[workloadID[0]][workloadID[1]]["opt_reaction_count"]=np.zeros((len(codeword_list),len(overhang_list))) #2d array for each workload file
    self._opt_codeword_results[workloadID[0]][workloadID[1]]["ideal_reaction_count"]=np.zeros((len(codeword_list),len(overhang_list))) 
    self._opt_codeword_results[workloadID[0]][workloadID[1]]["overhang_list"]=overhang_list
    self._opt_codeword_results[workloadID[0]][workloadID[1]]["codeword_list"]=codeword_list
    self._opt_codeword_results[workloadID[0]][workloadID[1]]["no_opt_reaction_count"]=np.zeros((len(codeword_list),len(overhang_list)))
    self._opt_codeword_results[workloadID[0]][workloadID[1]]["codewordsize_x_overhangcount"]=np.zeros((len(codeword_list),len(overhang_list)))
    self._opt_codeword_results[workloadID[0]][workloadID[1]]["ideal_height_map"]={} #going to have a height map for each codeword
    self._opt_codeword_results[workloadID[0]][workloadID[1]]["opt_height_map"]={}

    
    for codeword_index, codewordsize in enumerate(codeword_list):
        strand_length_codewords=int(math.ceil((strand_length_bytes*8)/codewordsize))
        index_length_codewords=int(math.ceil((index_bytes*8)/codewordsize))
        total_strand_length=strand_length_codewords+index_length_codewords
        
        print("---- starting codeword size {} ----".format(codewordsize))
        self._opt_codeword_results[workloadID[0]][workloadID[1]]["ideal_height_map"][codewordsize]=np.zeros((len(overhang_list),int(math.ceil(math.log(strand_length_codewords,2)))),dtype=np.uint) #build np array to be used as heat map
        self._opt_codeword_results[workloadID[0]][workloadID[1]]["opt_height_map"][codewordsize]=np.zeros((len(overhang_list),int(math.ceil(math.log(strand_length_codewords,2)))),dtype=np.uint)

        for overhang_index,overhang_count in enumerate(overhang_list):
            print("---- starting codeword size {} with overhang count {} ----".format(codewordsize,overhang_count))
            
            self._opt_codeword_results[workloadID[0]][workloadID[1]]["codewordsize_x_overhangcount"][codeword_index,overhang_index]=codewordsize*overhang_count
            dna_file=OverhangBitStringWriteDNAFile(primer5=self._primer5, formatid=self._format_ID, primer3=self._primer3, out_fd=output_filename,fsmd_abbrev='OH_BITSTRING_XXX',\
                                                   bits_per_block=codewordsize,strand_length=strand_length_codewords,num_overhangs=overhang_count)
            dna_file.write(data_buffer)
            dna_file.header_flush()
            strand_list=dna_file.get_strands() #get strands after encoding
            start_time=time.time()
            optimized_tree=tree.construct_tree_baseopt_lite(strand_list,overhang_count,int(math.ceil(math.log(overhang_count,4))),codewordsize,
                                                       total_strand_length, repair_strategy=0,
                                                       h_array=self._opt_codeword_results[workloadID[0]][workloadID[1]]["opt_height_map"][codewordsize][overhang_index][:],) #+2 for the number of bytes used for indexing
            self._opt_codeword_results[workloadID[0]][workloadID[1]]["opt_reaction_count"][codeword_index,overhang_index]=optimized_tree.order()
            del optimized_tree  #limit the amount of time trees are buffered to free up memory space for subsequent tree builds
            print("---- optimized tree build on {} took {} seconds ---".format(workloadID[1],time.time()-start_time))
            start_time=time.time()
            unoptimized_tree=tree.construct_tree_unoptimized_lite(strand_list,overhang_count,int(math.ceil(math.log(overhang_count,4))),codewordsize,
                                                             total_strand_length, repair_strategy=0)
            self._opt_codeword_results[workloadID[0]][workloadID[1]]["no_opt_reaction_count"][codeword_index,overhang_index]=unoptimized_tree.order()
            del unoptimized_tree
            print("---- no optimized tree build on {} took {} seconds ---".format(workloadID[1],time.time()-start_time))
            start_time=time.time()
            ideal_tree=tree.construct_tree_ideal_lite(strand_list,overhang_count,int(math.ceil(math.log(overhang_count,4))),codewordsize,total_strand_length, repair_strategy=0,
            h_array=self._opt_codeword_results[workloadID[0]][workloadID[1]]["ideal_height_map"][codewordsize][overhang_index][:])
            self._opt_codeword_results[workloadID[0]][workloadID[1]]["ideal_reaction_count"][codeword_index,overhang_index]=ideal_tree.order()
            del ideal_tree
            print("---- ideal tree build on {} took {} seconds ---".format(workloadID[1],time.time()-start_time))
            tlogger.debug('Finished building trees for '+output_filename)
            #collect the number of nodes in the constructed graphs, this equals the number of reactions the have to be performed

    #checkpoint the results here 
    picklefile=open(root_prefix+'1_results_codeword_opt','wb')
    pi.dump(self._opt_codeword_results,picklefile)
    picklefile.close()#store the ultimate results file
    
def analyze_opt_codewordsize(self):
    #analyze all workloads across different overhangs and data-in-block sizes
    for category in self._workloadDict:
        for work_file in self._workloadDict[category]:
            data_buffer=self._workloadDict[category][work_file]
            output_filename=category+'_codewordsize_'+work_file+'.output'
            self._sweep_overhangs_codewordsize(data_buffer,(category,work_file),output_filename)


def draw_opt_codewordsize(self):
    assert len(self._opt_codeword_results)>0
    #figure font settings
    font = {'family' : 'serif',
             'weight' : 'normal',
             'size'   : 6}
    matplotlib.rc('font',**font)
    for category in self._opt_codeword_results:
        #self._out_dir and self._1_bit_results should both have the same keys and key structure
        root_prefix=os.path.normpath(self._out_dir[category]["root"])+'/'

    

        #gather together data for each category and normalize results when necessary
        for _file in self._opt_codeword_results[category]: 
            file_prefix=os.path.normpath(self._out_dir[category][_file])+'/'
            resultsDict=self._opt_codeword_results[category][_file]

            overhang_list=resultsDict["overhang_list"]
            codeword_list=resultsDict["codeword_list"]
            
            #normalized data arrays
            opt_react_norm_array=np.zeros(resultsDict["opt_reaction_count"].shape)
            no_opt_react_norm_array=np.zeros(resultsDict["opt_reaction_count"].shape)
            ideal_react_norm_array=np.zeros(resultsDict["opt_reaction_count"].shape)
            
            #Normalized figures, data normalized internally to each type of optimization
            opt_react_norm_fig,opt_react_norm_axes=plt.subplots(nrows=1,ncols=1,figsize=(6,2.5),constrained_layout=True)
            no_opt_react_norm_fig,no_opt_react_norm_axes=plt.subplots(nrows=1,ncols=1,figsize=(6,2.5), constrained_layout=True)
            ideal_react_norm_fig,ideal_react_norm_axes=plt.subplots(nrows=1,ncols=1,figsize=(6,2.5), constrained_layout=True)

            #raw reaction counts for each optimization type 
            opt_react_raw_fig,opt_react_raw_axes=plt.subplots(nrows=1,ncols=1,figsize=(6,2.5), constrained_layout=True)
            no_opt_react_raw_fig,no_opt_react_raw_axes=plt.subplots(nrows=1,ncols=1,figsize=(6,2.5), constrained_layout=True)
            ideal_react_raw_fig,ideal_react_raw_axes=plt.subplots(nrows=1,ncols=1,figsize=(6,2.5), constrained_layout=True)

            #dump files for raw and normalized arrays
            opt_react_raw_dump=open(file_prefix+'codeword_over_opt_raw_'+category+'.csv','w+')
            no_opt_react_raw_dump=open(file_prefix+'codeword_over_no_opt_raw_'+category+'.csv','w+')
            opt_react_norm_dump=open(file_prefix+'codeword_over_opt_react_norm_'+category+'.csv','w+')
            ideal_react_norm_dump=open(file_prefix+'codeword_over_ideal_react_norm_'+category+'.csv','w+')
            no_opt_react_norm_dump=open(file_prefix+'codeword_over_no_opt_react_norm_'+category+'.csv','w+')
            ideal_react_raw_dump=open(file_prefix+'codeword_over_ideal_raw_'+category+'.csv','w+')

            #process the reaction count arrays
            for i in range(0,resultsDict["opt_reaction_count"].shape[0]):
                for j in range(0,resultsDict["opt_reaction_count"].shape[1]):
                    #going to normalize to [0,0] for
                    ideal_react_norm_array[i,j]=resultsDict["ideal_reaction_count"][i,j]/resultsDict["ideal_reaction_count"][0,0]
                    no_opt_react_norm_array[i,j]=resultsDict["no_opt_reaction_count"][i,j]/resultsDict["no_opt_reaction_count"][0,0]
                    opt_react_norm_array[i,j]=resultsDict["opt_reaction_count"][i,j]/resultsDict["opt_reaction_count"][0,0]

            #now we have the normalized arrays, and the non-normalized arrays, feed them to the heat map plotter
            title=category+" "+_file+" "

            #plot normalized arrays
            plt_util.plot_heatmap(ideal_react_norm_array,ideal_react_norm_axes,codeword_list,overhang_list,title,'Overhang Count','Codeword Size','Normalized Reactions',dumpFile=ideal_react_norm_dump)
            plt_util.plot_heatmap(opt_react_norm_array,opt_react_norm_axes,codeword_list,overhang_list,title,'Overhang Count','Codeword Size','Normalized Reactions',dumpFile=opt_react_norm_dump)
            plt_util.plot_heatmap(no_opt_react_norm_array,no_opt_react_norm_axes,codeword_list,overhang_list,title,'Overhang Count','Codeword Size','Normalized Reactions',dumpFile=no_opt_react_norm_dump)
            
            #plot the raw arrays
            plt_util.plot_heatmap(resultsDict["ideal_reaction_count"],ideal_react_raw_axes,codeword_list,overhang_list,'Codeword size vs Overhang Count','Overhang Count','Codeword Size','Reaction Count',dumpFile=ideal_react_raw_dump)
            plt_util.plot_heatmap(resultsDict["opt_reaction_count"],opt_react_raw_axes,codeword_list,overhang_list,'Codeword size vs Overhang Count','Overhang Count','Codeword Size','Reaction Count',dumpFile=opt_react_raw_dump)
            plt_util.plot_heatmap(resultsDict["no_opt_reaction_count"],no_opt_react_raw_axes,codeword_list,overhang_list,'Codeword size vs Overhang Count','Overhang Count','Codeword Size','Reaction Count',dumpFile=no_opt_react_raw_dump)
            
            
            #close dump files
            opt_react_norm_dump.close()
            no_opt_react_norm_dump.close()
            opt_react_raw_dump.close()
            no_opt_react_raw_dump.close()
            ideal_react_raw_dump.close()
            ideal_react_norm_dump.close()
            #save figures
            opt_react_norm_fig.savefig(file_prefix+'codeword_over_opt_react_norm_'+category+'.eps',format='eps')
            ideal_react_norm_fig.savefig(file_prefix+'codeword_over_ideal_react_norm_'+category+'.eps',format='eps')
            no_opt_react_norm_fig.savefig(file_prefix+'codeword_over_no_opt_react_norm_'+category+'.eps',format='eps')
            opt_react_raw_fig.savefig(file_prefix+'codeword_over_opt_react_raw_'+category+'.eps',format='eps')
            ideal_react_raw_fig.savefig(file_prefix+'codeword_over_opt_react_raw_'+category+'.eps',format='eps')
            no_opt_react_raw_fig.savefig(file_prefix+'codeword_over_no_opt_react_raw_'+category+'.eps',format='eps')
