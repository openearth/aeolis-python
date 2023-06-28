"""
Benchmarking the parallellization of the solve_matrix function using the concurrent.futures module.
"""

import scipy.sparse
import numpy as np

def solve_matrix(A,Ct_i):
    # solve the matrix
    result  = scipy.sparse.linalg.spsolve(A,Ct_i)

    return result.reshape(80,301)


def solve_matrix4(A,Ct_i):
    # solve the matrix
    result  = scipy.sparse.linalg.spsolve(A,Ct_i)

    return result

if __name__ == '__main__':
    import concurrent.futures
    import time
    # from itertools import izip

    ######################################
    ## 2 domains
    ######################################

    # Load data
    A = scipy.sparse.load_npz('parallell_solve/A.npz')
    yCt_i = np.load('parallell_solve/yCt_i.npy')

    start_time = time.time()

    # devide matrices in 2 domains
    decs =    [ 
        (A[0:24080,0:24080], yCt_i.flatten()[0:24080]), # first dec
        (A[21371:,21371:], yCt_i.flatten()[21371:]), # second dec
    ]

    results = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for dec, result in zip(decs, executor.map(solve_matrix, *zip(*decs) )):
            results.append(result)
        
    # combine the results
    comb_pool2 = np.vstack((
        results[0][0:75,:],
        results[1][4:,:]
        ))
    
    end_time = time.time()
    multiprocess_time = end_time - start_time

    solve_sp = scipy.sparse.linalg.spsolve(A, yCt_i.flatten()).reshape(151,301)

    diff_comb_pool2 = (comb_pool2-solve_sp).max()

    print(f'decompostion using Futures and 2 domains gives {multiprocess_time} seconds with an accuracy of {diff_comb_pool2}' )

    ######################################
    ## 4 domains
    ######################################                                                                                                  

     # Load data
    A = scipy.sparse.load_npz('parallell_solve/A.npz')
    yCt_i = np.load('parallell_solve/yCt_i.npy')

    start_time = time.time()

    # devide matrices in 4 domains
    decs =    [ 
        (A[0:13545,0:13545], yCt_i.flatten()[0:13545]), # first dec
        (A[12040:25585,12040:25585], yCt_i.flatten()[12040:25585]), # second dec
        (A[24080:37625,24080:37625], yCt_i.flatten()[24080:37625]), # third dec
        (A[36120:45451,36120:45451], yCt_i.flatten()[36120:45451]) # fourth dec
    ]
    
    results = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for dec, result in zip(decs, executor.map(solve_matrix4, *zip(*decs) )):
            results.append(result)
  
    reshaped_results = [result.reshape(45,301) for result in results[:-1]]
    last_result = results[-1].reshape(31,301)

    # combine the results
    comb_pool4 = np.vstack((
        reshaped_results[0][0:41,:],reshaped_results[1][1:-5,:],reshaped_results[2][:-5,:],last_result[:,:]
        ))
    
    end_time = time.time()
    multiprocess_time = end_time - start_time
    
    diff_comb_pool4 = (comb_pool4-solve_sp).max()

    print(f'decompostion using Pool and 4 domains gives {multiprocess_time} seconds with an accuracy of {diff_comb_pool2}' )