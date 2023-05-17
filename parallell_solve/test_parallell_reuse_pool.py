import scipy.sparse
import numpy as np
import timeit

def solve_matrix(A,Ct_i):
    # solve the matrix
    result  = scipy.sparse.linalg.spsolve(A,Ct_i)

    return result.reshape(80,301)


def solve_matrix4(A,Ct_i):
    # solve the matrix
    result  = scipy.sparse.linalg.spsolve(A,Ct_i)

    return result

if __name__ == '__main__':
    from multiprocessing import Pool, Value, Array
    import time
    import matplotlib.pyplot as plt

    ######################################
    ## 4 domains with Pool Reuse
    ######################################                                                                                                  
     # Load data
    A = scipy.sparse.load_npz('parallell_solve/A.npz')
    yCt_i = np.load('parallell_solve/yCt_i.npy')

    # create a pool of workers
    pool = Pool()

    start_time = time.time()
    # using a look to simulate a loop over time steps
    for i in range(1):
        # devide matrices in 4 domains
        decs =    [ 
            (A[0:13545,0:13545], yCt_i.flatten()[0:13545]), # first dec
            (A[12040:25585,12040:25585], yCt_i.flatten()[12040:25585]), # second dec
            (A[24080:37625,24080:37625], yCt_i.flatten()[24080:37625]), # third dec
            (A[36120:45451,36120:45451], yCt_i.flatten()[36120:45451]) # fourth dec
        ]
        
        # using same pool 
        results = pool.starmap(solve_matrix4, decs)

        reshaped_results = [result.reshape(45,301) for result in results[:-1]]
        last_result = results[-1].reshape(31,301)

        # combine the results
        comb_pool4 = np.vstack((
            reshaped_results[0][0:41,:],reshaped_results[1][1:-5,:],reshaped_results[2][:-5,:],last_result[:,:]
            ))

    # close the pool and wait for the work to finish    
    pool.close()
    pool.join()

    end_time = time.time()

    print('decompostion using Pool with 4 domains solved sequentially gives %.2f seconds ' %  ( end_time-start_time ) )

    # plot the results
    time_spsolve = [[1,2,4, 8, 16, 32, 64, 128, 256, 512, 1024],[0.08, 0.13, 0.33, 0.52, 1.04, 2.05, 4.11, 8.35, 16.33, 32.55, 67.74]]
    time_comb2 = [
        [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024],
        [0.5, 0.11, 0.21, 0.4, 0.79, 1.60, 3.14, 6.25, 12.43, 25.93, 53.02]]
    time_comb4 = [[1,2,4, 8, 16, 32, 64, 128, 256, 512, 1024],[0.5, 0.09, 0.17, 0.33, 0.63, 1.21, 2.40, 4.97, 9.52, 20.72, 42.44]]
    time_comb4_pool = [[1,2,4, 8, 16, 32, 64, 128, 256, 512, 1024],[0.3, 0.06, 0.11, 0.22, 0.41, 0.86, 1.68, 3.45, 6.88, 13.70, 27.55]]

    plt.plot(time_spsolve[0], time_spsolve[1], label='current (sp_solve)', marker='+', color='gray')
    plt.plot(time_comb2[0], time_comb2[1], label='2 domains', marker='o')
    plt.plot(time_comb4[0], time_comb4[1], label='4 domains', marker='x', color='r')
    plt.plot(time_comb4_pool[0], time_comb4_pool[1], label='4 domains + reuse pool', marker='*', color='g')
    plt.xlabel('iterations')
    plt.ylabel('Time [s]')
    plt.legend()
    plt.show()

