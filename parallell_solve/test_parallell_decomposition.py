import scipy.sparse
import numpy as np
import timeit

def solve_matrix(A,Ct_i):
    # solve the matrix
    result  = scipy.sparse.linalg.spsolve(A,Ct_i)

    return result.reshape(80,301)

# similar to the above, but without reshaping
def solve_matrix4(A,Ct_i):
    # solve the matrix
    result  = scipy.sparse.linalg.spsolve(A,Ct_i)

    return result

if __name__ == '__main__':
    from multiprocessing import Pool
    import matplotlib.pyplot as plt

    ######################################
    ## 2 domains
    ######################################

    # Load data
    A = scipy.sparse.load_npz('parallell_solve/A.npz')
    yCt_i = np.load('parallell_solve/yCt_i.npy')

    # a function for the parallelization
    # this is only for compatibility with timeit
    def pool_2():

        # devide matrices in 2 domains
        decs =    [ 
            (A[0:24080,0:24080], yCt_i.flatten()[0:24080]), # first dec
            (A[21371:,21371:], yCt_i.flatten()[21371:]), # second dec
        ]

        with Pool(processes=5) as pool:
            results = pool.starmap(solve_matrix, decs)
    
        # combine the results
        comb_pool2 = np.vstack((
            results[0][0:75,:],
            results[1][4:,:]
            ))
        return comb_pool2
    

    solve_sp = scipy.sparse.linalg.spsolve(A, yCt_i.flatten()).reshape(151,301)

    comb_pool2 = pool_2()
    diff_comb_pool2 = (comb_pool2-solve_sp).max()


    ######################################
    ## 4 domains
    ######################################                                                                                                  

     # Load data
    A = scipy.sparse.load_npz('parallell_solve/A.npz')
    yCt_i = np.load('parallell_solve/yCt_i.npy')

    # a function for the parallelization
    # this is only for compatibility with timeit
    def pool_4():
        # devide matrices in 4 domains
        decs =    [ 
            (A[0:13545,0:13545], yCt_i.flatten()[0:13545]), # first dec
            (A[12040:25585,12040:25585], yCt_i.flatten()[12040:25585]), # second dec
            (A[24080:37625,24080:37625], yCt_i.flatten()[24080:37625]), # third dec
            (A[36120:45451,36120:45451], yCt_i.flatten()[36120:45451]) # fourth dec
        ]
        
        with Pool(processes=5) as pool:
            results = pool.starmap(solve_matrix4, decs)
    
        reshaped_results = [result.reshape(45,301) for result in results[:-1]]
        last_result = results[-1].reshape(31,301)

        # combine the results
        comb_pool4 = np.vstack((
            reshaped_results[0][0:41,:],reshaped_results[1][1:-5,:],reshaped_results[2][:-5,:],last_result[:,:]
            ))

        return comb_pool4

    comb_pool4 = pool_4()
    
    diff_comb_pool4 = (comb_pool4-solve_sp).max()

    print('decompostion using Pool with 2 domains solved sequentially gives %.2f seconds with an accuracy of %.2e ' %  ( timeit.timeit('comb_pool2 = pool_2()', number=1,globals=globals()), diff_comb_pool2 ) )
    print('decompostion using Pool with 4 domains solved sequentially gives %.2f seconds with an accuracy of %.2e ' %  ( timeit.timeit('comb_pool4 = pool_4()', number=1,globals=globals()), diff_comb_pool4 ) )

    # plot the results
    time_spsolve = [[1,2,4, 8, 16, 32, 64, 128, 256, 512, 1024],[0.08, 0.13, 0.33, 0.52, 1.04, 2.05, 4.11, 8.35, 16.33, 32.55, 67.74]]
    time_comb2 = [[1,2,4, 8, 16, 32, 64, 128, 256, 512, 1024],[0.05, 0.11, 0.21, 0.4, 0.79, 1.60, 3.14, 6.25, 12.43, 25.93, 53.02]]
    time_comb4 = [[1,2,4, 8, 16, 32, 64, 128, 256, 512, 1024],[0.05, 0.09, 0.17, 0.33, 0.63, 1.21, 2.40, 4.97, 9.52, 20.72, 42.44]]

    plt.plot(time_spsolve[0], time_spsolve[1], label='current (sp_solve)', marker='+', color='gray')
    plt.plot(time_comb2[0], time_comb2[1], label='2 domains', marker='o')
    plt.plot(time_comb4[0], time_comb4[1], label='4 domains', marker='x', color='r')
    plt.xlabel('Time [s]')
    plt.ylabel('iterations')
    plt.legend()
    plt.show()

