import numpy as np
import matplotlib
import pandas as pd
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import os
import json
import copy
import pandas as pd
from math import sqrt

#Utilities for plotting data



def pandas_barchart(dump_frame,y_list,title,ylabel,y_ticks,y_log,x_tick_labels,yrange):
     #write data into bloks in the heat map
    font_xticks = {'family' : 'serif',
             'weight' : 'normal',
                   'size'   : 8.7}
    font_yticks= {'family' : 'serif',
             'weight' : 'normal',
                  'size'   : 8.7}
    font_ylabel= {'family' : 'serif',
             'weight' : 'normal',
                  'size'   : 8.7}
    #color=['0.2','0.1','0.5','0.6','0.9']
    #if len(color)<len(x):
    #       color = np.linspace(0, 1, len(x))

    bar_axes=dump_frame.plot.bar(rot=0,figsize=(9,2),y=y_list,width=0.7,edgecolor='k',logy=y_log)
    bar_fig=bar_axes.get_figure()
    if ylabel!=False:
        bar_axes.set_ylabel(ylabel,fontsize=font_ylabel['size'])

    if y_ticks==False:
        bar_axes.get_yaxis().set_visible(False)

    if x_tick_labels==False:
        bar_axes.xaxis.set_ticklabels([])

    bar_axes.set_title(title,fontsize=font_ylabel['size'])
    bar_axes.set_ylim(yrange[0],yrange[1])
    matplotlib.rc('hatch', color='k', linewidth=0.5)
    patterns =['xxxxx','//////','......','-----','***']
    bars=bar_axes.patches
    hatches = [p for p in patterns for i in range(len(dump_frame))]
    for bar, hatch in zip(bars, hatches):
        bar.set_hatch(hatch)
        
    bar_axes.get_legend().remove()
    bar_axes.legend(fontsize=6,loc='upper left',ncol=len(dump_frame))

    bar_axes.tick_params(labelsize=font_xticks['size'])
    return bar_fig



#wrapper function to plot multiple components on the same axis
#labelSet provides a label for each component at a corresponding location
def plot_components_wrapper(x,comp_axes,components,title,xlabel,ylabel,labelSet=None,linestyleSet=None,colorSet=None,log=False,dumpFile=None,linewidth_=None, markerSet=None,markeverySet=None):
    if dumpFile!=None:#make data frame then dump it 
        dump_dict={}
        for label_index,label in enumerate(labelSet):
            if isinstance(components[label_index],(np.ndarray,list)):
                dump_dict[label]=components[label_index]
            else:
                dump_dict[label]=list(components[label_index])

        if isinstance(x,(np.ndarray,list)):
            dump_frame=pd.DataFrame(dump_dict,index=x)
        else:
            dump_frame=pd.DataFrame(dump_dict,index=list(x))
        dump_frame.to_csv(dumpFile)
        

    if linewidth_ == None:
        linewidth_ = 1.8
    prop_cycle=plt.rcParams['axes.prop_cycle']
    defaultColors=prop_cycle.by_key()['color']
    #defaultColors=[str(t) for t in np.linspace(0.0,0.6,len(components))]
    #print defaultColors
    linestyleArray=[]
    labelArray=[]
    colorArray=[]
    markerArray=[]
    markeveryArray=[]
    
    if markeverySet==None:
        for comp in components:
            markeveryArray.append(100)
    else:
        for m in markeverySet:
            markeveryArray.append(m)
        while len(markeveryArray)<len(components):
            markeveryArray.append(100)

    
    if markerSet ==None:
        for comp in components:
            markerArray.append(None)
    else:
        for m in markerSet:
            if m == '-' or m=='--':
                markerArray.append(None)
            else:
                markerArray.append(m)
            
        while len(markerArray)<len(components):
            markerArray.append(None)
    

    #Check styling arrays
    if linestyleSet is None:
        for comp in components:
            linestyleArray.append('-') #default to solid lines
    else:
        for l in linestyleSet:
            linestyleArray.append(l)
        while len(linestyleArray) < len(components):
            linestyleArray.append('-') #fill in missing line styles with solid lines

    if colorSet is None:
        colorArray=defaultColors #default to the default color cycle
        colorIndex=0
        while len(colorArray)<len(components):
            colorArray.append(defaultColors[colorIndex%len(colorArray)])
            colorIndex+=1
    else:
        for c in colorSet:
            colorArray.append(c)
        colorIndex=0
        while len(colorArray) < len(components): #fill in missing colors with default colors
            colorArray.append(defaultColors[colorIndex])
            colorIndex+=1
    if labelSet is None:
        for comp in components:
            labelArray.append('')
    else:
        for compID, comp in enumerate(components):
            #print "compID {}".format(compID)
            labelArray.append(labelSet[compID])

    assert len(components)<=len(linestyleArray)
    assert len(components)<=len(colorArray)
    assert len(components)<=len(markerArray)
    assert len(components)<=len(markeveryArray)
    assert len(components)<=len(labelArray)

            
    for compID,comp in enumerate(components):
        comp_axes.plot(x,comp,label=labelArray[compID],color=colorArray[compID],linestyle=linestyleArray[compID],linewidth=linewidth_,marker=markerArray[compID],markevery=markeveryArray[compID],markersize=2)
    comp_axes.set_xlabel(xlabel)

    
    comp_axes.set_ylabel(ylabel)
    comp_axes.set_title(title)
    if labelSet is not None: #only display legend if there are multiple lines to display
        comp_axes.legend()
    if log is True:
        comp_axes.set_yscale('log')


def plot_bargraph(fig, ax, rectSet, labelSet,xlabel,ylabel,xnameSet,log,title,dumpFile=None):
    #use frames to dump to a csv file
    if dumpFile!=None:
        dump_dict={}
        for rectIndex, rect in enumerate(rectSet):
            if isinstance(rect,(np.ndarray,list)):
                dump_dict[labelSet[rectIndex]]=rect
            else:
                dump_dict[labelSet[rectIndex]]=list(rect)
            dump_frame=pd.DataFrame(dump_dict,index=list(xnameSet))
            dump_frame.to_csv(dumpFile)

            
    font = {'family' : 'serif',
        'weight' : 'normal',
            'size'   : 6}
    matplotlib.rc('font',**font)


    ax.grid(zorder=0) # Add grid lines
    ax.grid(which='minor',linestyle=':',zorder=0)
    
    bar_width = 0.7 / len(rectSet)
    index = np.arange(len(rectSet[0]))
    #print index
    color = np.linspace(0, 1, len(rectSet))

    dump_list=[[]]*len(xnameSet) #have x name set empty lists, append to each one the rectSets data 
    
    for rectIndex, rect in enumerate(rectSet):
        rect_plt = ax.bar(index+rectIndex*bar_width, rect, bar_width,
                       color=str(color[rectIndex]),
                       label=labelSet[rectIndex],
                          edgecolor="black",zorder=2)
                
    # Add labels
    plt.xlabel(xlabel,fontsize=6) # Add x lines
    plt.ylabel(ylabel,fontsize=6) # Add y lines
    # Add tick names
    plt.xticks(index + bar_width*len(rectSet)/2, xnameSet,fontsize=6)
    ax.xaxis.set_ticks_position('none') 
    plt.yticks(fontsize=6)


    # Add legend
    plt.legend(fontsize=4, loc="upper left", bbox_to_anchor=[0,1],ncol=2)
    if log==True:
        ax.set_yscale('log')
        ax.set_ylim(top=100)



    

#plots a heatmap for a 2D array 
def plot_heatmap(data,ax,ylabelSet,xlabelSet,figLabel,xLabel,yLabel,cbarLabel,cmap="summer",dumpFile=None, scientific=True,fontsize=None):
    #dump heat map out by using data frames
    if dumpFile!=None:
        dump_frame=pd.DataFrame(data,index=list(ylabelSet),columns=list(xlabelSet))
        dump_frame.to_csv(dumpFile)

    heatmap=ax.imshow(data,cmap,aspect='auto')
    #print data
    #show all ticks
    ax.set_xticks(np.arange(len(xlabelSet)))
    ax.set_yticks(np.arange(len(ylabelSet)))
    ax.set_xticklabels(xlabelSet)
    ax.set_yticklabels(ylabelSet)

    if fontsize==None:
        _fontsize=5.5
    else:
        _fontsize=fontsize
    #this line will rotate the x axis labels to a 45 degree angle 

    #write data into bloks in the heat map
    font = {'family' : 'serif',
             'weight' : 'normal',
             'size'   : _fontsize}
    for i in range(len(ylabelSet)):
        for j in range(len(xlabelSet)):
            if scientific: #plot text in scientific notation
                #convert number of scientific notation
                s = '{x:0.{ndp:d}e}'.format(x=np.float32(data[i,j]), ndp=2)
                m, e = s.split('e')
                sci_string=r'{m:s}\times 10^{{{e:d}}}'.format(m=m, e=int(e))
                text=ax.text(j,i,"${0:s}$".format(sci_string),ha="center",va="center",color="black", fontdict=font)
            else:
                text=ax.text(j,i,data[i,j],ha="center",va="center",color="black",fontdict=font)
    #set axis labels
    ax.set_ylabel(yLabel)
    ax.set_xlabel(xLabel)
    ax.set_title(figLabel)

    #create colorbar
    cbar=ax.figure.colorbar(heatmap,ax=ax)
    cbar.ax.set_ylabel(cbarLabel,rotation=-90,va="bottom")
    
    
    
