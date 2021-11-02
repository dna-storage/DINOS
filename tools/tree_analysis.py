'''
Author: Kevin Volkel

filename: tree_analysis.py

Description: Call analysis functions to explore optimizations on reaction trees for different files

'''


if __name__ == "__main__":
    import argparse
    from overhang.opt_analysis import tree_analysis
    
    parser=argparse.ArgumentParser(description="Analyze optimized and unoptimized trees for different workloads")

    parser.add_argument('--w_dir',dest="w_dir",action="store", default="tree_workloads",help="Directory to use for workloads")
    parser.add_argument('--out_dir',dest='out_dir',action="store",default="",help="Directory to store output figures and results")
    parser.add_argument('--cw_size',dest='cw_size',action="store_true",default=False,help="Do codeword size experiment")
    parser.add_argument('--1_bit',dest='_1_bit',action="store_true",default=False,help="Do 1 bit analysis")
    parser.add_argument('--pickled_results', dest='pickled_results',action="store",default=None,help="path to pickled results that will be plotted by the analysis chosen")
    args = parser.parse_args()

    t_analysis=tree_analysis(w_dir=args.w_dir, random_data=False, out_dir=args.out_dir, pickle=args.pickled_results)

    if args._1_bit:
        #codeword size = 1 bit experiments
        if args.pickled_results is None:
            t_analysis.analyze_1_bit() #if path specified, short circuit to drawing results
        t_analysis.draw_1_bit()
    if args.cw_size:
        #codeword size experiments
        if args.pickled_results is None:
            t_analysis.analyze_opt_codewordsize()
        t_analysis.draw_opt_codewordsize()
        
    
    
