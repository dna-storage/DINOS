import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import overhang.tree as tree
import overhang.reaction_node as node
import logging
from overhang.dnastorage_utils.system.dnafile import *
import os
import sys
import shutil
import math
import numpy as np
import overhang.plot_utils.plot_utils as plt_util
import time
import pickle as pi
import gc
import pandas as pd
import scipy


tlogger=logging.getLogger('dna.overhang.tools.tree_analysis')
tlogger.addHandler(logging.NullHandler())


def _sweep_overhangs_1_bit(self,data_buffer,workloadID,debug_fname=None): #workloadID is a tuple indicating the complete name of the workload
    output_filename=debug_fname
    overhang_list=[3,5,9,17,65]
    #overhang_list=[9,17,65]
    strand_length_bytes=509 #508 bytes to account for the extra 4 bytes of indexing
    index_bytes=3
    root_prefix=os.path.normpath(self._out_dir[workloadID[0]]["root"])+'/'
    #set up result dictionary for the file we are analyzing 
    if workloadID[0] not in self._1_bit_results:
        self._1_bit_results[workloadID[0]]={}
        if workloadID[1] not in self._1_bit_results[workloadID[0]]:
            self._1_bit_results[workloadID[0]][workloadID[1]]={}
    else:
        if workloadID[1] not in self._1_bit_results[workloadID[0]]:
            self._1_bit_results[workloadID[0]][workloadID[1]]={}

    self._1_bit_results[workloadID[0]][workloadID[1]]["opt_reaction_count"]=[] #array of integers
    self._1_bit_results[workloadID[0]][workloadID[1]]["transform_reaction_count"]=[]
    self._1_bit_results[workloadID[0]][workloadID[1]]["ideal_reaction_count"]=[] #array of integers
    self._1_bit_results[workloadID[0]][workloadID[1]]["overhang_array"]=overhang_list
    self._1_bit_results[workloadID[0]][workloadID[1]]["no_opt_reaction_count"]=[] #array of integers
    self._1_bit_results[workloadID[0]][workloadID[1]]["rotate_reaction_count"]=[]
    
    self._1_bit_results[workloadID[0]][workloadID[1]]["ideal_height_map"]=np.zeros((len(overhang_list),int(math.ceil(math.log((strand_length_bytes+index_bytes)*8,2)))),dtype=np.uint) #build np array to be used as heat map
    self._1_bit_results[workloadID[0]][workloadID[1]]["opt_height_map"]=np.zeros((len(overhang_list),int(math.ceil(math.log((strand_length_bytes+index_bytes)*8,2)))),dtype=np.uint)
    
    for overhang_index,overhang_count in enumerate(overhang_list):
        print("{} {}".format(workloadID,overhang_count))
        dna_file=OverhangBitStringWriteDNAFile(primer5=self._primer5, formatid=self._format_ID, primer3=self._primer3, out_fd=output_filename,fsmd_abbrev='OH_BITSTRING_XXX',\
                                                   bits_per_block=1,strand_length=strand_length_bytes*8,num_overhangs=overhang_count)
        dna_file.write(data_buffer)
        dna_file.header_flush()
        strand_list=dna_file.get_strands() #get strands after encoding
     
        start_time=time.time()
        transform_tree=tree.construct_tree_transform_lite(strand_list,overhang_count,int(math.ceil(math.log(overhang_count,4))),1,
                                                          (index_bytes+strand_length_bytes)*8,h_array=None) #+4 for the number of bytes used for indexing
        self._1_bit_results[workloadID[0]][workloadID[1]]["transform_reaction_count"].append(transform_tree.order())
        #print transform_tree.order()
        del transform_tree
        print("---- transform  tree build on {} took {} seconds ---".format(workloadID[1],time.time()-start_time))


        start_time=time.time()
        optimized_tree=tree.construct_tree_baseopt_lite_w(strand_list,overhang_count,int(math.ceil(math.log(overhang_count,4))),1,
                                                       (index_bytes+strand_length_bytes)*8,h_array=self._1_bit_results[workloadID[0]][workloadID[1]]["opt_height_map"][overhang_index][:]) #+4 for the number of bytes used for indexing
        self._1_bit_results[workloadID[0]][workloadID[1]]["opt_reaction_count"].append(optimized_tree.order())
        #print optimized_tree.order()
        opt_hash=optimized_tree.strand_hash_table #grab the hash table
        del optimized_tree
        
        print("---- optimized tree build on {} took {} seconds ---".format(workloadID[1],time.time()-start_time))

        
        start_time=time.time()
        rotate_tree=tree.construct_tree_rotate_lite(strand_list,overhang_count,int(math.ceil(math.log(overhang_count,4))),1,
                                                          (index_bytes+strand_length_bytes)*8,h_array=None,opt_dictionary=opt_hash) #+4 for the number of bytes used for indexing
        self._1_bit_results[workloadID[0]][workloadID[1]]["rotate_reaction_count"].append(rotate_tree.order())
        #print transform_tree.order()
        del rotate_tree
        print("---- rotate  tree build on {} took {} seconds ---".format(workloadID[1],time.time()-start_time))


        start_time=time.time()
        unoptimized_tree=tree.construct_tree_unoptimized_lite(strand_list,overhang_count,int(math.ceil(math.log(overhang_count,4))),1,
                                                         (index_bytes+strand_length_bytes)*8)
        self._1_bit_results[workloadID[0]][workloadID[1]]["no_opt_reaction_count"].append(unoptimized_tree.order())
        del unoptimized_tree
        gc.collect()
        print("---- no optimized tree build on {} took {} seconds ---".format(workloadID[1],time.time()-start_time))
        start_time=time.time()
        ideal_tree=tree.construct_tree_ideal_lite(strand_list,overhang_count,int(math.ceil(math.log(overhang_count,4))),1,
                                                 (index_bytes+strand_length_bytes)*8,h_array=self._1_bit_results[workloadID[0]][workloadID[1]]["ideal_height_map"][overhang_index][:])
        self._1_bit_results[workloadID[0]][workloadID[1]]["ideal_reaction_count"].append(ideal_tree.order())
        del ideal_tree
        gc.collect()
        print("---- ideal tree build on {} took {} seconds ---".format(workloadID[1],time.time()-start_time))
        tlogger.debug('Finished building trees for '+output_filename)
        #collect the number of nodes in the constructed graphs, this equals the number of reactions the have to be performed 
        sys.stdout.flush()
    #checkpoint results by pickling data for each file
    picklefile=open(root_prefix+'1_bit_results','wb')
    pi.dump(self._1_bit_results,picklefile)
    picklefile.close()#store the ultimate results file
    
def analyze_1_bit(self):
    #analyze all workloads across different overhangs and data-in-block sizes
    for category in self._workloadDict:
        for work_file in self._workloadDict[category]:
            data_buffer=self._workloadDict[category][work_file]
            output_filename=category+'_'+work_file+'.output'
            self._sweep_overhangs_1_bit(data_buffer,(category,work_file),output_filename)


#draw figures and dump results to csv 
def draw_1_bit(self):
    assert len(self._1_bit_results)>0
    
    #figure font settings
    font = {'family' : 'serif',
             'weight' : 'normal',
             'size'   : 6}
    matplotlib.rc('font',**font)
    #draw results from sweeping the number of overhangs for 1 bit building blocks  
    #create line graphs for optimized and unoptimized reaction counts 
    for category in self._1_bit_results:
        #self._out_dir and self._1_bit_results should both have the same keys and key structure
        root_prefix=os.path.normpath(self._out_dir[category]["root"])+'/'
        #create plots and axes for each category of data 
        opt_react_norm_fig,opt_react_norm_axes=plt.subplots(nrows=1,ncols=1,figsize=(5.2,2.5),constrained_layout=True)
        no_opt_react_norm_fig,no_opt_react_norm_axes=plt.subplots(nrows=1,ncols=1,figsize=(6,2.5), constrained_layout=True)
        opt_to_no_opt_normalized_fig, opt_to_no_opt_normalized_axes=plt.subplots(nrows=1,ncols=1,figsize=(5.2,2), constrained_layout=True)
        ideal_to_no_opt_fig, ideal_to_no_opt_axes=plt.subplots(nrows=1,ncols=1,figsize=(9,2.5),constrained_layout=True)
        opt_react_raw_fig,opt_react_raw_axes=plt.subplots(nrows=1,ncols=1,figsize=(3,2.5), constrained_layout=True)
        no_opt_react_raw_fig,no_opt_react_raw_axes=plt.subplots(nrows=1,ncols=1,figsize=(6,2.5), constrained_layout=True)
        ideal_react_norm_fig,ideal_react_norm_axes=plt.subplots(nrows=1,ncols=1,figsize=(5.2,2.5), constrained_layout=True)
        ideal_react_raw_fig,ideal_react_raw_axes=plt.subplots(nrows=1,ncols=1,figsize=(6,2.5), constrained_layout=True)

        #transform optimization graphs
        transform_react_norm_fig,transform_react_norm_axes=plt.subplots(nrows=1,ncols=1,figsize=(5.2,2.5), constrained_layout=True)
        transform_react_raw_fig,transform_react_raw_axes=plt.subplots(nrows=1,ncols=1,figsize=(6,2.5), constrained_layout=True)
        transform_to_no_opt_fig, transform_to_no_opt_axes=plt.subplots(nrows=1,ncols=1,figsize=(5.2,2.5),constrained_layout=True)

        #rotate optimization graphs
        rotate_react_norm_fig,rotate_react_norm_axes=plt.subplots(nrows=1,ncols=1,figsize=(5.2,2.5), constrained_layout=True)
        rotate_react_raw_fig,rotate_react_raw_axes=plt.subplots(nrows=1,ncols=1,figsize=(6,2.5), constrained_layout=True)
        rotate_to_no_opt_fig, rotate_to_no_opt_axes=plt.subplots(nrows=1,ncols=1,figsize=(5.2,2.5),constrained_layout=True)
                

        #dump files
        opt_react_raw_dump=open(root_prefix+'1_bit_opt_raw_'+category+'.csv','w+')
        no_opt_react_raw_dump=open(root_prefix+'1_bit_no_opt_raw_'+category+'.csv','w+')
        opt_react_norm_dump=open(root_prefix+'1_bit_opt_react_norm_'+category+'.csv','w+')
        ideal_react_norm_dump=open(root_prefix+'1_bit_ideal_react_norm_'+category+'.csv','w+')
        no_opt_react_norm_dump=open(root_prefix+'1_bit_no_opt_react_norm_'+category+'.csv','w+')
        opt_to_no_opt_dump=open(root_prefix+'1_bit_opt_to_no_opt_'+category+'.csv','w+')
        ideal_to_no_opt_dump=open(root_prefix+'1_bit_ideal_to_no_opt_'+category+'.csv','w+')
        ideal_react_raw_dump=open(root_prefix+'1_bit_ideal_raw_'+category+'.csv','w+')

        #transform dump files
        transform_to_no_opt_dump=open(root_prefix+'1_bit_transform_to_no_opt_'+category+'.csv','w+')
        transform_react_raw_dump=open(root_prefix+'1_bit_transform_raw_'+category+'.csv','w+')
        transform_react_norm_dump=open(root_prefix+'1_bit_transform_react_norm_'+category+'.csv','w+')

        #rotate dump files
        rotate_to_no_opt_dump=open(root_prefix+'1_bit_rotate_to_no_opt_'+category+'.csv','w+')
        rotate_react_raw_dump=open(root_prefix+'1_bit_rotate_raw_'+category+'.csv','w+')
        rotate_react_norm_dump=open(root_prefix+'1_bit_rotate_react_norm_'+category+'.csv','w+')

        
        #arrays to hold data to be plotted 
        file_name_array=[]
        opt_react_norm_data_array=[]
        no_opt_react_norm_data_array=[]
        opt_to_no_opt_data_array=[]
        ideal_react_norm_data_array=[]
        ideal_to_no_opt_data_array=[]
        opt_react_raw_data_array=[]
        ideal_react_raw_data_array=[]
        no_opt_react_raw_data_array=[]

        #transform arrays
        transform_react_norm_data_array=[]
        transform_to_no_opt_data_array=[]
        transform_react_raw_data_array=[]

        #transform arrays
        rotate_react_norm_data_array=[]
        rotate_to_no_opt_data_array=[]
        rotate_react_raw_data_array=[]

        

        
        #gather together data for each category and normalize results when necessary
        for _file in self._1_bit_results[category]:
            file_prefix=os.path.normpath(self._out_dir[category][_file])+'/'
            opt_heat_map_dump=open(file_prefix+'opt_heat_map'+category+"_"+_file+'.csv','w+')
            ideal_heat_map_dump=open(file_prefix+'ideal_heat_map'+category+'_'+_file+'.csv','w+')

            #heat maps
            ideal_heat_fig,ideal_heat_axes=plt.subplots(nrows=1,ncols=1,figsize=(6,2.5), constrained_layout=True)
            opt_heat_fig,opt_heat_axes=plt.subplots(nrows=1,ncols=1,figsize=(6,2.5),constrained_layout=True)
            resultsDict=self._1_bit_results[category][_file]
            #get data for heat maps
            ideal_heat_data=resultsDict["ideal_height_map"]
            opt_heat_data=resultsDict["opt_height_map"]
            #get data for line and group them together into arrays
            overhang_array=resultsDict["overhang_array"]
            opt_react_raw_data_array.append(resultsDict["opt_reaction_count"])
            no_opt_react_raw_data_array.append(resultsDict["no_opt_reaction_count"])
            ideal_react_raw_data_array.append(resultsDict["ideal_reaction_count"])
            
            opt_react_norm_data=[float(_)/float(resultsDict["opt_reaction_count"][0]) for _ in resultsDict["opt_reaction_count"]]
            ideal_react_norm_data=[float(_)/float(resultsDict["ideal_reaction_count"][0]) for _ in resultsDict["ideal_reaction_count"]]
            no_opt_react_norm_data=[float(_)/float(resultsDict["no_opt_reaction_count"][0]) for _ in resultsDict["no_opt_reaction_count"]]
            opt_to_no_opt_data=[float(opt)/float(no_opt) for opt,no_opt in zip(resultsDict["opt_reaction_count"],resultsDict["no_opt_reaction_count"])]
            ideal_to_no_opt_data=[float(ideal)/float(no_opt) for ideal,no_opt in zip(resultsDict["ideal_reaction_count"],resultsDict["no_opt_reaction_count"])]

            #transform calcs
            transform_react_raw_data_array.append(resultsDict["transform_reaction_count"])
            transform_react_norm_data=[float(_)/float(resultsDict["transform_reaction_count"][0]) for _ in resultsDict["transform_reaction_count"]]
            transform_to_no_opt_data=[float(ideal)/float(no_opt) for ideal,no_opt in zip(resultsDict["transform_reaction_count"],resultsDict["no_opt_reaction_count"])]

            #rotate calcs
            rotate_react_raw_data_array.append(resultsDict["rotate_reaction_count"])
            rotate_react_norm_data=[float(_)/float(resultsDict["rotate_reaction_count"][0]) for _ in resultsDict["rotate_reaction_count"]]
            rotate_to_no_opt_data=[float(ideal)/float(no_opt) for ideal,no_opt in zip(resultsDict["rotate_reaction_count"],resultsDict["no_opt_reaction_count"])]
            
            
            if "." in _file:
                file_name_array.append(_file.split('.')[0])
            else:
                file_name_array.append(_file)
            opt_react_norm_data_array.append(opt_react_norm_data)
            ideal_react_norm_data_array.append(ideal_react_norm_data)
            no_opt_react_norm_data_array.append(no_opt_react_norm_data)
            opt_to_no_opt_data_array.append(opt_to_no_opt_data)
            ideal_to_no_opt_data_array.append(ideal_to_no_opt_data)

            #append transform data
            transform_to_no_opt_data_array.append(transform_to_no_opt_data)
            transform_react_norm_data_array.append(transform_react_norm_data)

            #append rotate data
            rotate_to_no_opt_data_array.append(rotate_to_no_opt_data)
            rotate_react_norm_data_array.append(rotate_react_norm_data)

            
            #plot the heat maps, 1 for each file and 1 for ideal/opt (may want to compress opt and ideal into one: 1 array or 2 subplots)
            heat_fig_label=category+" "+_file+" "
            heat_xLabel="height in tree"
            heat_yLabel="number of overhangs"
            cbar_label="match count"
            plt_util.plot_heatmap(ideal_heat_data,ideal_heat_axes,overhang_array,range(1,len(ideal_heat_data[0][:])+1),heat_fig_label+" ideal",heat_xLabel,heat_yLabel,cbar_label,dumpFile=ideal_heat_map_dump, fontsize=4)
            plt_util.plot_heatmap(opt_heat_data,opt_heat_axes,overhang_array,range(1,len(ideal_heat_data[0][:])+1),heat_fig_label+" opt",heat_xLabel,heat_yLabel,cbar_label,dumpFile=opt_heat_map_dump,fontsize=4)
            #save a heat map for each file studied
            opt_heat_fig.savefig(file_prefix+'opt_heat_map_'+category+"_"+_file+'.eps',format='eps')
            ideal_heat_fig.savefig(file_prefix+'ideal_heat_map_'+category+"_"+_file+'.eps',format='eps')
            opt_heat_map_dump.close()
            ideal_heat_map_dump.close()
        markerSet=(None,None,'o','^','x','D',None,None,None,'H','+','X')
        linestyleSet=('-','-','-','-','--','-','-','-','-','--','--','--')
        markeverySet=[1]*12
        #draw line charts

        
        

        #optimized reactions raw graph
        plt_util.plot_components_wrapper(overhang_array,opt_react_raw_axes,opt_react_raw_data_array,category+" (opt raw reaction count)","Number of Overhangs","Number of Reactions",labelSet=file_name_array,linestyleSet=linestyleSet,markerSet=markerSet,dumpFile=opt_react_raw_dump,linewidth_=0.7,markeverySet=markeverySet)
        chartBox = opt_react_raw_axes.get_position()
        opt_react_raw_axes.set_position([chartBox.x0+0.1, chartBox.y0+0.1, chartBox.width*0.6, chartBox.height*0.9])
        opt_react_raw_axes.get_legend().remove()
        opt_react_raw_fig.legend(fontsize=4,loc='center',bbox_to_anchor=(0.835,0.7),ncol=1)



        #ideal reactions raw graph
        plt_util.plot_components_wrapper(overhang_array,ideal_react_raw_axes,ideal_react_raw_data_array,category+" (ideal raw reaction count)","Number of Overhangs","Number of Reactions",labelSet=file_name_array,linestyleSet=linestyleSet,markerSet=markerSet,dumpFile=ideal_react_raw_dump,linewidth_=0.7,markeverySet=markeverySet)
        chartBox = opt_react_raw_axes.get_position()
        ideal_react_raw_axes.set_position([chartBox.x0+0.1, chartBox.y0+0.1, chartBox.width*0.6, chartBox.height*0.9])
        ideal_react_raw_axes.get_legend().remove()
        ideal_react_raw_fig.legend(fontsize=4,loc='center',bbox_to_anchor=(0.835,0.7),ncol=1)

        #trasform raw graph
        plt_util.plot_components_wrapper(overhang_array,transform_react_raw_axes,transform_react_raw_data_array,category+" (transform raw reaction count)","Number of Overhangs","Number of Reactions",labelSet=file_name_array,linestyleSet=linestyleSet,markerSet=markerSet,dumpFile=transform_react_raw_dump,linewidth_=0.7,markeverySet=markeverySet)
        chartBox = transform_react_raw_axes.get_position()
        transform_react_raw_axes.set_position([chartBox.x0+0.1, chartBox.y0+0.1, chartBox.width*0.6, chartBox.height*0.9])
        transform_react_raw_axes.get_legend().remove()
        transform_react_raw_fig.legend(fontsize=4,loc='center',bbox_to_anchor=(0.835,0.7),ncol=1)


        
        #trasform raw graph
        plt_util.plot_components_wrapper(overhang_array,rotate_react_raw_axes,rotate_react_raw_data_array,category+" (rotate raw reaction count)","Number of Overhangs","Number of Reactions",labelSet=file_name_array,linestyleSet=linestyleSet,markerSet=markerSet,dumpFile=rotate_react_raw_dump,linewidth_=0.7,markeverySet=markeverySet)
        chartBox = rotate_react_raw_axes.get_position()
        rotate_react_raw_axes.set_position([chartBox.x0+0.1, chartBox.y0+0.1, chartBox.width*0.6, chartBox.height*0.9])
        rotate_react_raw_axes.get_legend().remove()
        rotate_react_raw_fig.legend(fontsize=4,loc='center',bbox_to_anchor=(0.835,0.7),ncol=1)
        


        

        #No optimized reactions raw graph
        plt_util.plot_components_wrapper(overhang_array,no_opt_react_raw_axes,no_opt_react_raw_data_array,category+" (no-opt raw reaction count)","Number of Overhangs","Number of Reactions",labelSet=file_name_array,linestyleSet=linestyleSet,markerSet=markerSet,dumpFile=no_opt_react_raw_dump,linewidth_=0.7,markeverySet=markeverySet)
        chartBox = no_opt_react_raw_axes.get_position()
        no_opt_react_raw_axes.set_position([chartBox.x0+0.1, chartBox.y0+0.1, chartBox.width*0.6, chartBox.height*0.9])
        no_opt_react_raw_axes.get_legend().remove()
        no_opt_react_raw_fig.legend(fontsize=4,loc='center',bbox_to_anchor=(0.4,0.7),ncol=1)


        #opt self normalized graph
        plt_util.plot_components_wrapper(overhang_array,opt_react_norm_axes,opt_react_norm_data_array,category+" (opt normalized)","Number of Overhangs","Normalized Reactions",labelSet=file_name_array,linestyleSet=linestyleSet,markerSet=markerSet,dumpFile=opt_react_norm_dump,linewidth_=0.7,markeverySet=markeverySet)
        chartBox = opt_react_norm_axes.get_position()
        opt_react_norm_axes.set_position([chartBox.x0+0.1, chartBox.y0+0.1, chartBox.width*0.6, chartBox.height*0.9])
        opt_react_norm_axes.get_legend().remove()
        opt_react_norm_fig.legend(fontsize=4,loc='center',bbox_to_anchor=(0.5,0.68),ncol=1)

        #no opt self normalized graph
        plt_util.plot_components_wrapper(overhang_array,no_opt_react_norm_axes,no_opt_react_norm_data_array,category+" (no opt normalized)","Number of Overhangs","Normalized Reactions",labelSet=file_name_array,linestyleSet=linestyleSet,markerSet=markerSet,dumpFile=no_opt_react_norm_dump,linewidth_=0.7,markeverySet=markeverySet)
        chartBox = no_opt_react_norm_axes.get_position()
        no_opt_react_norm_axes.set_position([chartBox.x0+0.1, chartBox.y0+0.1, chartBox.width*0.6, chartBox.height*0.9])
        no_opt_react_norm_axes.get_legend().remove()
        no_opt_react_norm_fig.legend(fontsize=4,loc='center',bbox_to_anchor=(0.61,0.43),ncol=1)



        #ideal self normalized graph
        plt_util.plot_components_wrapper(overhang_array,ideal_react_norm_axes,ideal_react_norm_data_array,category+" (ideal normalized)","Number of Overhangs","Normalized Reactions",labelSet=file_name_array,linestyleSet=linestyleSet,markerSet=markerSet,dumpFile=ideal_react_norm_dump,linewidth_=0.7,markeverySet=markeverySet)
        chartBox = ideal_react_norm_axes.get_position()
        ideal_react_norm_axes.set_position([chartBox.x0+0.1, chartBox.y0+0.1, chartBox.width*0.6, chartBox.height*0.9])
        ideal_react_norm_axes.get_legend().remove()
        #ideal_react_norm_fig.legend(fontsize=4,loc='center',bbox_to_anchor=(0.835,0.7),ncol=1)


        #transform self normalized
        plt_util.plot_components_wrapper(overhang_array,transform_react_norm_axes,transform_react_norm_data_array,category+" (transform normalized)","Number of Overhangs","Normalized Reactions",labelSet=file_name_array,linestyleSet=linestyleSet,markerSet=markerSet,dumpFile=transform_react_norm_dump,linewidth_=0.7,markeverySet=markeverySet)
        chartBox = transform_react_norm_axes.get_position()
        transform_react_norm_axes.set_position([chartBox.x0+0.1, chartBox.y0+0.1, chartBox.width*0.6, chartBox.height*0.9])
        transform_react_norm_axes.get_legend().remove()
        #ideal_react_norm_fig.legend(fontsize=4,loc='center',bb


        

        #rotate self normalized
        plt_util.plot_components_wrapper(overhang_array,rotate_react_norm_axes,rotate_react_norm_data_array,category+" (rotate normalized)","Number of Overhangs","Normalized Reactions",labelSet=file_name_array,linestyleSet=linestyleSet,markerSet=markerSet,dumpFile=rotate_react_norm_dump,linewidth_=0.7,markeverySet=markeverySet)
        chartBox = rotate_react_norm_axes.get_position()
        rotate_react_norm_axes.set_position([chartBox.x0+0.1, chartBox.y0+0.1, chartBox.width*0.6, chartBox.height*0.9])
        rotate_react_norm_axes.get_legend().remove()
        #ideal_react_norm_fig.legend(fontsize=4,loc='center',bb
               

        #opt reactions normalized to no-opt reactions
        plt_util.plot_components_wrapper(overhang_array,opt_to_no_opt_normalized_axes,opt_to_no_opt_data_array,category + " (opt/no-opt)","Number of Overhangs","Normalized Reactions",labelSet=file_name_array,linestyleSet=linestyleSet,markerSet=markerSet,dumpFile=opt_to_no_opt_dump,linewidth_=0.7,markeverySet=markeverySet)
        chartBox = opt_to_no_opt_normalized_axes.get_position()
        opt_to_no_opt_normalized_axes.set_position([chartBox.x0+0.1, chartBox.y0+0.1, chartBox.width*0.6, chartBox.height*0.9])
        opt_to_no_opt_normalized_axes.get_legend().remove()
        #opt_to_no_opt_normalized_fig.legend(fontsize=4.3,loc='center',bbox_to_anchor=(0.3,0.68),ncol=1)
        opt_to_no_opt_normalized_fig.legend(fontsize=4.3,loc='center',bbox_to_anchor=(0.61,0.43),ncol=1)
        
        #ideal reactions normalized to no-opt reactions
        plt_util.plot_components_wrapper(overhang_array,ideal_to_no_opt_axes,ideal_to_no_opt_data_array,category + " (ideal/no-opt)","Number of Overhangs","Normalized Reactions",labelSet=file_name_array,linestyleSet=linestyleSet,markerSet=markerSet,dumpFile=ideal_to_no_opt_dump,linewidth_=0.7,markeverySet=markeverySet)
        chartBox = ideal_to_no_opt_axes.get_position()
        ideal_to_no_opt_axes.set_position([chartBox.x0+0.1, chartBox.y0+0.1, chartBox.width*0.6, chartBox.height*0.9])
        ideal_to_no_opt_axes.get_legend().remove()
        #ideal_to_no_opt_fig.legend(fontsize=4,loc='center',bbox_to_anchor=(0.729,0.7),ncol=1)

        #transform reactions normalized to no-opt reactions
        plt_util.plot_components_wrapper(overhang_array,transform_to_no_opt_axes,transform_to_no_opt_data_array,category + " (transform/no-opt)","Number of Overhangs","Normalized Reactions",labelSet=file_name_array,linestyleSet=linestyleSet,markerSet=markerSet,dumpFile=transform_to_no_opt_dump,linewidth_=0.7,markeverySet=markeverySet)
        chartBox = transform_to_no_opt_axes.get_position()
        transform_to_no_opt_axes.set_position([chartBox.x0+0.1, chartBox.y0+0.1, chartBox.width*0.6, chartBox.height*0.9])
        transform_to_no_opt_axes.get_legend().remove()
        #ideal_to_no_opt_fig.legend(fontsize=4,loc='center',bbox_to_anchor=(0.729,0.7),ncol=1)

        
        #rotate reactions normalized to no-opt reactions
        plt_util.plot_components_wrapper(overhang_array,rotate_to_no_opt_axes,rotate_to_no_opt_data_array,category + " (rotate/no-opt)","Number of Overhangs","Normalized Reactions",labelSet=file_name_array,linestyleSet=linestyleSet,markerSet=markerSet,dumpFile=rotate_to_no_opt_dump,linewidth_=0.7,markeverySet=markeverySet)
        chartBox = rotate_to_no_opt_axes.get_position()
        rotate_to_no_opt_axes.set_position([chartBox.x0+0.1, chartBox.y0+0.1, chartBox.width*0.6, chartBox.height*0.9])
        rotate_to_no_opt_axes.get_legend().remove()
        #ideal_to_no_opt_fig.legend(fontsize=4,loc='center',bbox_to_anchor=(0.729,0.7),ncol=1)

        
        
        
        #close dump files
        opt_react_norm_dump.close()
        no_opt_react_norm_dump.close()
        opt_to_no_opt_dump.close()
        ideal_to_no_opt_dump.close()
        opt_react_raw_dump.close()
        no_opt_react_raw_dump.close()
        ideal_react_raw_dump.close()
        ideal_react_norm_dump.close()

        #close transform dumps
        transform_react_raw_dump.close()
        transform_react_norm_dump.close()
        transform_to_no_opt_dump.close()

        #close rotate dumps
        rotate_react_raw_dump.close()
        rotate_react_norm_dump.close()
        rotate_to_no_opt_dump.close()

        
        
        #ideal to no opt bar graph
        df=pd.read_csv(root_prefix+'1_bit_ideal_to_no_opt_'+category+'.csv',index_col=0)
        #calculate the geomean on the dataframe
        geomean=[_ for _ in scipy.stats.gmean(df.iloc[:,:],axis=1)]
        _dict={"geomean":geomean}
        _df=pd.DataFrame(_dict,index=df.index)
        df_ideal=pd.concat([df, _df],axis=1, sort=False)
        ideal_to_no_opt_bar_fig=plt_util.pandas_barchart(df_ideal.transpose(),[3,5,9,17,65],"(ideal/no-opt)","Normalized Reactions",True,True,True,(0.007,1.1))
        ideal_to_no_opt_bar_fig.axes[0].get_legend().remove()
        legend=ideal_to_no_opt_bar_fig.legend(fontsize=8,loc='center',bbox_to_anchor=(0.935,0.5),ncol=1,title=r'$|\mathcal{O}|$',markerscale=1.2)
        legend.get_title().set_fontsize('10') #legend 'Title' fontsize



        #transform to no opt bar graph
        df=pd.read_csv(root_prefix+'1_bit_transform_to_no_opt_'+category+'.csv',index_col=0)
        #calculate the geomean on the dataframe
        geomean=[_ for _ in scipy.stats.gmean(df.iloc[:,:],axis=1)]
        _dict={"geomean":geomean}
        _df=pd.DataFrame(_dict,index=df.index)
        df_transform=pd.concat([df, _df],axis=1, sort=False)
        if 'zpaq' in category:
            transform_to_no_opt_bar_fig=plt_util.pandas_barchart(df_transform.transpose(),[3,5,9,17,65],"(transform/no-opt) zpaq","Normalized Reactions",True,True,True,(0.007,1.1))
        else:
            transform_to_no_opt_bar_fig=plt_util.pandas_barchart(df_transform.transpose(),[3,5,9,17,65],"transform/no-opt","Normalized Reactions",True,True,True,(0.007,1.1))

        transform_to_no_opt_bar_fig.axes[0].get_legend().remove()
        transform_to_no_opt_bar_fig.legend(fontsize=6,loc='center',bbox_to_anchor=(0.28,0.92),ncol=len(df.index))

        
        #rotate to no opt bar graph
        df=pd.read_csv(root_prefix+'1_bit_rotate_to_no_opt_'+category+'.csv',index_col=0)
        #calculate the geomean on the dataframe
        geomean=[_ for _ in scipy.stats.gmean(df.iloc[:,:],axis=1)]
        _dict={"geomean":geomean}
        _df=pd.DataFrame(_dict,index=df.index)
        df_rotate=pd.concat([df, _df],axis=1, sort=False)
        if 'zpaq' in category:
            rotate_to_no_opt_bar_fig=plt_util.pandas_barchart(df_rotate.transpose(),[3,5,9,17,65],"(rotate/no-opt) zpaq","Normalized Reactions",True,True,True,(0.007,1.1))
        else:
            rotate_to_no_opt_bar_fig=plt_util.pandas_barchart(df_rotate.transpose(),[3,5,9,17,65],"rotate/no-opt","Normalized Reactions",True,True,True,(0.007,1.1))
        rotate_to_no_opt_bar_fig.axes[0].get_legend().remove()
        rotate_to_no_opt_bar_fig.legend(fontsize=6,loc='center',bbox_to_anchor=(0.28,0.92),ncol=len(df.index))
 
        #opt to no opt bar graph
        df=pd.read_csv(root_prefix+'1_bit_opt_to_no_opt_'+category+'.csv',index_col=0)
        #calculate the geomean on the dataframe
        geomean=[_ for _ in scipy.stats.gmean(df.iloc[:,:],axis=1)]
        _dict={"geomean":geomean}
        _df=pd.DataFrame(_dict,index=df.index)
        df_opt=pd.concat([df, _df],axis=1, sort=False)
        if 'zpaq' in category:
            opt_to_no_opt_bar_fig=plt_util.pandas_barchart(df_opt.transpose(),[3,5,9,17,65],'(base-opt/no-opt) zpaq',"Normalized Reactions",True,True,True,(0.007,1.1))
        else:
            opt_to_no_opt_bar_fig=plt_util.pandas_barchart(df_opt.transpose(),[3,5,9,17,65],'base-opt/no-opt',"Normalized Reactions",True,True,True,(0.007,1.1))

        opt_to_no_opt_bar_fig.axes[0].get_legend().remove()
        legend=opt_to_no_opt_bar_fig.legend(fontsize=8,loc='center',bbox_to_anchor=(0.935,0.5),ncol=1,title=r'$|\mathcal{O}|$',markerscale=1.2)
        legend.get_title().set_fontsize('10') #legend 'Title' fontsize

        #ideal to opt bar graph
        if 'zpaq' in category:
            ideal_to_opt_bar_fig=plt_util.pandas_barchart((df_ideal/df_opt).transpose(),[3,5,9,17,65],"(ideal/base-opt) zpaq","Normalized Reactions",True,False,True,(0,1.05))
        else:
            ideal_to_opt_bar_fig=plt_util.pandas_barchart((df_ideal/df_opt).transpose(),[3,5,9,17,65],"ideal/base-opt","Normalized Reactions",True,False,True,(0,1.05))
        ideal_to_opt_bar_fig.axes[0].get_legend().remove()
        legend=ideal_to_opt_bar_fig.legend(fontsize=8,loc='center',bbox_to_anchor=(0.935,0.5),ncol=1,title=r'$|\mathcal{O}|$',markerscale=1.2)
        legend.get_title().set_fontsize('10') #legend 'Title' fontsize

        #rotate to opt - 1 bar graph
        rotate_to_opt=df_rotate/df_opt
        percent_change=rotate_to_opt.subtract(1).multiply(100.0)
        print(percent_change)
        percent_change.to_csv(root_prefix+"rotate_percent_change.csv",index=False)
        if 'zpaq' in category:
            rotate_to_opt_bar_fig=plt_util.pandas_barchart(percent_change.transpose(),[9,17,65],"(Alignment Percent Change) zpaq","Percent Change \n Relative to base-opt",True,False,True,(-42,1))
        else:
            rotate_to_opt_bar_fig=plt_util.pandas_barchart(percent_change.transpose(),[9,17,65],"Alignment Percent Change","Percent Change \n Relative to base-opt",True,False,True,(-37,6))
            rotate_to_opt_bar_fig.axes[0].yaxis.set_ticks([-30,-20,-10,0,5])

        rotate_to_opt_bar_fig.axes[0].get_legend().remove()
        rotate_to_opt_bar_fig.axes[0].axhline(linewidth=0.5, color='k')
        legend=rotate_to_opt_bar_fig.legend(fontsize=8,loc='center',bbox_to_anchor=(0.935,0.5),ncol=1,title=r'$|\mathcal{O}|$',markerscale=1.2)
        legend.get_title().set_fontsize('10') #legend 'Title' fontsize
        rotate_to_opt_bar_fig.axes[0].tick_params(axis='y', labelsize=7.0)
        rotate_to_opt_bar_fig.axes[0].tick_params(axis='x', labelsize=8.7)

        
        #ideal to transform bar graph
        if 'zpaq' in category:
            ideal_to_transform_bar_fig=plt_util.pandas_barchart((df_ideal/df_transform).transpose(),[3,5,9,17,65],"(ideal/transform) zpaq","Normalized Reactions",True,False,True,(0,1.05))
        else:
            ideal_to_transform_bar_fig=plt_util.pandas_barchart((df_ideal/df_transform).transpose(),[3,5,9,17,65],"ideal/transform","Normalized Reactions",True,False,True,(0,1.05))
        ideal_to_transform_bar_fig.axes[0].get_legend().remove()
        ideal_to_transform_bar_fig.legend(fontsize=6,loc='center',bbox_to_anchor=(0.28,0.92),ncol=len(df.index))
        

        #ideal to rotate bar graph
        if 'zpaq' in category:
            ideal_to_rotate_bar_fig=plt_util.pandas_barchart((df_ideal/df_rotate).transpose(),[3,5,9,17,65],"(ideal/rotate) zpaq","Normalized Reactions",True,False,True,(0,1.05))
        else:
            ideal_to_rotate_bar_fig=plt_util.pandas_barchart((df_ideal/df_rotate).transpose(),[3,5,9,17,65],"ideal/rotate","Normalized Reactions",True,False,True,(0,1.05))
        ideal_to_rotate_bar_fig.axes[0].get_legend().remove()
        ideal_to_rotate_bar_fig.legend(fontsize=6,loc='center',bbox_to_anchor=(0.28,0.92),ncol=len(df.index))
        


        
        #save figures
        opt_react_norm_fig.savefig(root_prefix+'1_bit_opt_react_norm_'+category+'.eps',format='eps')
        ideal_react_norm_fig.savefig(root_prefix+'1_bit_ideal_react_norm_'+category+'.eps',format='eps')
        no_opt_react_norm_fig.savefig(root_prefix+'1_bit_no_opt_react_norm_'+category+'.eps',format='eps')
        opt_to_no_opt_normalized_fig.savefig(root_prefix+'1_bit_opt_to_no_opt_'+category+'.eps',format='eps')
        opt_to_no_opt_bar_fig.savefig(root_prefix+'1_bit_opt_to_no_opt_bar_'+category+'.pdf',format='pdf')
        ideal_to_no_opt_fig.savefig(root_prefix+'1_bit_ideal_to_no_opt_'+category+'.eps',format='eps')
        ideal_to_no_opt_bar_fig.savefig(root_prefix+'1_bit_ideal_to_no_opt_bar_'+category+'.pdf',format='pdf')
        opt_react_raw_fig.savefig(root_prefix+'1_bit_opt_react_raw_'+category+'.eps',format='eps')
        ideal_react_raw_fig.savefig(root_prefix+'1_bit_opt_react_raw_'+category+'.eps',format='eps')
        no_opt_react_raw_fig.savefig(root_prefix+'1_bit_no_opt_react_raw_'+category+'.eps',format='eps')
        ideal_to_opt_bar_fig.savefig(root_prefix+'1_bit_ideal_to_opt_bar_'+category+'.pdf',format='pdf')

        #save transform figs
        transform_to_no_opt_fig.savefig(root_prefix+'1_bit_transform_to_no_opt_'+category+'.eps',format='eps')
        transform_to_no_opt_bar_fig.savefig(root_prefix+'1_bit_transform_to_no_opt_bar_'+category+'.pdf',format='pdf')
        ideal_to_transform_bar_fig.savefig(root_prefix+'1_bit_ideal_to_transform_bar_'+category+'.pdf',format='pdf')
        transform_react_raw_fig.savefig(root_prefix+'1_bit_transform_react_raw_'+category+'.eps',format='eps')
        transform_react_norm_fig.savefig(root_prefix+'1_bit_transform_react_norm_'+category+'.eps',format='eps')

        
        #save rotate figs
        rotate_to_no_opt_fig.savefig(root_prefix+'1_bit_rotate_to_no_opt_'+category+'.eps',format='eps')
        rotate_to_no_opt_bar_fig.savefig(root_prefix+'1_bit_rotate_to_no_opt_bar_'+category+'.pdf',format='pdf')
        rotate_to_opt_bar_fig.savefig(root_prefix+'1_bit_rotate_to_opt_'+category+'.pdf',format='pdf')
        ideal_to_rotate_bar_fig.savefig(root_prefix+'1_bit_ideal_to_rotate_bar_'+category+'.pdf',format='pdf')
        rotate_react_raw_fig.savefig(root_prefix+'1_bit_rotate_react_raw_'+category+'.eps',format='eps')
        rotate_react_norm_fig.savefig(root_prefix+'1_bit_rotate_react_norm_'+category+'.eps',format='eps')
        
