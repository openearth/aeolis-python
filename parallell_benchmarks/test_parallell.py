import scipy.sparse
import numpy as np
import timeit
import matplotlib.pyplot as plt





def decomp2(Ct_i,A,yCt_i):
    # devide domain specifying overlap
    # we have a 151 by 301 domain which we devide by 2 with 5 cells ovelap
    # the first domain is 80*301 = 24080 indeces large
    # we perform the sub matrix solve
    Ct_i_dec1 = scipy.sparse.linalg.spsolve(A[0:24080,0:24080], yCt_i.flatten()[0:24080]).reshape(80,301)
    #the second domain starts on index 71*301 = 21371 to the end
    # we perform the sub matrix solve
    Ct_i_dec2 = scipy.sparse.linalg.spsolve(A[21371:,21371:], yCt_i.flatten()[21371:]).reshape(80,301)
    # now we combine the matrices again
    comb = np.vstack((Ct_i_dec1[0:75,:],Ct_i_dec2[4:,:]))
    return comb



def decomp4(Ct_i,A,yCt_i):
    # devide domain specifying overlap
    # we have a 151 by 301 domain which we devide by 4 with about 5 cells ovelap
    # the first domain is 45*301 = 13545 indeces large
    # we perform the sub matrix solve
    Ct_i_dec1 = scipy.sparse.linalg.spsolve(A[0:13545,0:13545], yCt_i.flatten()[0:13545]).reshape(45,301)
    #the second domain starts on index 40*301 = 12040 until 85*301 = 25585
    # we perform the sub matrix solve
    Ct_i_dec2 = scipy.sparse.linalg.spsolve(A[12040:25585,12040:25585], yCt_i.flatten()[12040:25585]).reshape(45,301)
    #the third domain starts on index 80*301 = 24080 until 125*301 = 37625
    # we perform the sub matrix solve
    Ct_i_dec3 = scipy.sparse.linalg.spsolve(A[24080:37625,24080:37625], yCt_i.flatten()[24080:37625]).reshape(45,301)
    #the fourth domain starts on index 120*301 = 36120 until 151*301 = 45451
    #we perform the sub matrix solve
    Ct_i_dec4 = scipy.sparse.linalg.spsolve(A[36120:45451,36120:45451], yCt_i.flatten()[36120:45451]).reshape(31,301)
    # now we combine the matrices again
    comb4 = np.vstack((Ct_i_dec1[0:41,:],Ct_i_dec2[1:-5,:],Ct_i_dec3[:-5,:],Ct_i_dec4[:,:]))
    return comb4


def solve_matrix(A,Ct_i):
    # solve the matrix
    result  = scipy.sparse.linalg.spsolve(A,Ct_i)

    return result.reshape(80,301)


if __name__ == '__main__':
    from multiprocessing import Pool, Value, Array
    import time

    Ct_i = np.load('parallell_solve/Ct_i.npy')
    A = scipy.sparse.load_npz('parallell_solve/A.npz')
    yCt_i = np.load('parallell_solve/yCt_i.npy')

    def pool_2():

        start_time = time.time()
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
        
        # end_time = time.time()
        # print(end_time)
        # multiprocess_time = end_time - start_time

    solve_sp = scipy.sparse.linalg.spsolve(A, yCt_i.flatten()).reshape(151,301)
    comb2 = decomp2(Ct_i,A,yCt_i)
    comb4 = decomp4(Ct_i,A,yCt_i)

    diff_comb2 = (comb2-solve_sp).max()
    diff_comb4 = (comb4-solve_sp).max()

    comb_pool2 = pool_2()
    diff_comb_pool2 = (comb_pool2-solve_sp).max()

    start_time = time.time()
    solve_sp = scipy.sparse.linalg.spsolve(A, yCt_i.flatten()).reshape(151,301)
    end_time = time.time()
    solve_sp_time = end_time - start_time

    print('')
    print('spsolve gives %.2f seconds, this is the benchmark' % timeit.timeit('solve_sp = scipy.sparse.linalg.spsolve(A, yCt_i.flatten()).reshape(151,301)',number=1024,globals=globals()))
    print('decompostion in 2 domains solved sequentially gives %.2f seconds with an accuracy of %.2e ' %  (timeit.timeit('comb2 = decomp2(Ct_i,A,yCt_i)',number=1,globals=globals()), diff_comb2))
    print('decompostion in 4 domains solved sequentially gives %.2f seconds with an accuracy of %.2e ' %  ( timeit.timeit('comb4 = decomp4(Ct_i,A,yCt_i)', number=1,globals=globals()), diff_comb4 ) )

    print('decompostion using Pool with 2 domains solved sequentially gives %.2f seconds with an accuracy of %.2e ' %  ( timeit.timeit('comb_pool2 = pool_2()', number=1,globals=globals()), diff_comb_pool2 ) )

                                                                                                                         
    print('')
    print('next step is to solve for domains in parallell')
    print('')

#plt.imshow(comb4-solve_sp)
