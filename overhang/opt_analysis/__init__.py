import logging
import os
import shutil
import pickle as pi

tlogger=logging.getLogger('dna.overhang.tools.tree_analysis')
tlogger.addHandler(logging.NullHandler())

class tree_analysis: #class for analyzing reaction trees based on input files
    from overhang.opt_analysis._1_bit_analysis import _sweep_overhangs_1_bit,analyze_1_bit,draw_1_bit #import analysis dedicated to 1 bit codewords
    from overhang.opt_analysis._opt_codeword_analysis import _sweep_overhangs_codewordsize, analyze_opt_codewordsize, draw_opt_codewordsize
    #import analysis to find optimal overhang/codeword combination 
    
    def __init__(self,**kwargs):
        self._workloadDict={} #dictionary to keep buffers for workloads
        self._primer3=""
        self._primer5=""
        self._format_ID=0x3000
        self._out_dir={} #this dictionary takes in a category directory and file name and outputs a path to the appropriate output directory
        self._1_bit_results={}#dictionary to hold results of 1 bit blocks  analysis, structure is set up like workloadDict, a hierarchical manner to reflect the category and workload, this should be general results
        self._opt_codeword_results={} #this dictionary is a container of results for analyzing an optimal codewordsize/overhang combination


        #load already existing pickled data to short circuit launching the actual analysis
        if kwargs['pickle'] is not None:
            #assume this is the path to the pickled data for 1 bit
            pickle_file=open(kwargs['pickle'],'rb')
            self._1_bit_results=pi.load(pickle_file)
        if kwargs['pickle'] is not None:
            pickle_file=open(kwargs['pickle'],'rb')
            self._opt_codeword_results=pi.load(pickle_file)
            
        if kwargs.has_key('out_dir'):
            if not os.path.exists(kwargs['out_dir']):
                #shutil.rmtree(kwargs['out_dir'])
                #make the path if it does not exist
                os.mkdir(kwargs['out_dir'])

        if kwargs.has_key('w_dir'):
            w_dir=kwargs['w_dir']
            #directory structure is expected to have a subdirectory for each category of input files: e.g images, text, random should all be in different subdirectories
            for _dir in os.listdir(w_dir):
                category_path=os.path.join(w_dir,os.path.basename(os.path.normpath(_dir)))
                tlogger.debug(category_path + '--seen at workloadDict construction')
                if os.path.isdir(category_path):
                    self._out_dir[_dir]={}
                    if kwargs.has_key('out_dir'):
                        self._out_dir[_dir]["root"]=os.path.join(kwargs['out_dir'],os.path.normpath(_dir))#root result directory for 
                        if not os.path.exists(self._out_dir[_dir]["root"]):
                            #shutil.rmtree(self._out_dir[_dir]["root"])
                            os.mkdir(self._out_dir[_dir]["root"])
                    if _dir not in self._workloadDict:
                        self._workloadDict[_dir]={}
                        for work_file in os.listdir(category_path):
                            work_file_path=os.path.join(category_path,os.path.basename(os.path.normpath(work_file)))
                            tlogger.debug(work_file_path+'--seen at workloadDict construction') 
                            if os.path.isdir(work_file_path):
                                tlogger.warning('found directory in sub-directory of workload directory')
                            else:
                                if kwargs.has_key('out_dir'):
                                    self._out_dir[_dir][work_file]=os.path.join(self._out_dir[_dir]["root"],os.path.normpath(work_file))#file specific result path
                                    if not os.path.exists(self._out_dir[_dir][work_file]):
                                        #shutil.rmtree(self._out_dir[_dir][work_file])
                                        os.mkdir(self._out_dir[_dir][work_file])
                                if work_file not in self._workloadDict[_dir]:
                                    with open(work_file_path) as data: #read in whole file into buffer may need to change in the future
                                        self._workloadDict[_dir][work_file]=data.read()
                                        data.close()
                else:
                    tlogger.warning('found non-directory file at top of workload directory')
        else: tlogger.critical('No workload dictionary specified')
        tlogger.debug('Loaded all workloads')
        if kwargs.has_key('random_data'):
            if kwargs['random_data']==True:
                tlogger.debug('Generating random data buffers')
                #add random data to the workload Dictionary
                self._workloadDict['random']={}
                size_list=[100, 1000, 10000] #hardcode to 100 kB, 1000 kB, 10000 kB for now
                for size in size_list:
                    self._workloadDict['random'][str(size)]=bytearray(os.urandom(size))
