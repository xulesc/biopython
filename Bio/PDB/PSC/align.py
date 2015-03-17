# Copyright (C) 2015, Anuj Sharma (anuj.sharma80@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
"""

from .munkres import run_munkres
import numpy as np
from scipy.spatial import distance
import math

class Munkres:
  def run_munkres(self, cost_matrix, do_max = 0):
    orig_shape = cost_matrix.shape
    if orig_shape[0] != orig_shape[1]:
      self.cost_matrix = np.copy(cost_matrix)
      max_dim = max(orig_shape)
      self.cost_matrix.resize(max_dim, max_dim)
    else:
      self.cost_matrix = cost_matrix
    
    c = np.array(run_munkres(self.cost_matrix, do_max)).reshape(self.cost_matrix.shape)
    return c[0:orig_shape[0], 0:orig_shape[1]]
  
  def self_test(self):
    print 'Test 1 - forces matrix copy'
    a = np.array([[1, 2, 3, 6], [2, 4, 6, 5], [3, 6, 9, 1]]).astype(np.int32)  
    print 'a = %s' %a
    print 'c = %s' %(self.run_munkres(a))
    print 'Test 2 - no matrix copy'
    a = np.array([[1, 2, 3], [2, 4, 6], [3, 6, 9]]).astype(np.int32)  
    print 'a = %s' %a
    print 'c = %s' %(self.run_munkres(a))

class DP:
    #http://www.avatar.se/molbioinfo2001/dynprog/adv_dynamic.html
    def dp(self, gapcost, S):
        w = gapcost; (rows, cols) = S.shape
        # init
        M = np.zeros((rows + 1, cols + 1), dtype=np.float)
        ## fill
        for i in xrange(1, rows + 1):
            for j in xrange(1, cols + 1):
                M[i][j] = max(M[i-1][j-1]+S[i-1][j-1], M[i][j-1]+w, M[i-1][j]+w)
        X = rows; Y = cols
        ## traceback
        alignment = []
        while X > 0 and Y > 0:
            ## diagonal, left, up
            u = [(X-1,Y-1),(X-1,Y),(X,Y-1)]
            v = np.array([M[X-1][Y-1] + S[X-1][Y-1], M[X-1][Y] + w, M[X][Y-1] + w])
            ## which value is equal to the current value
            i = np.where(v==M[X][Y])[0][0]
            (X,Y) = u[i]
            if i == 0:
                alignment.append((X, Y))
        alignment.reverse()
        return alignment                                                                                                                    

    def dp2(self, gapcost, S):
        w = gapcost; (rows, cols) = S.shape
        ## init
        dir = np.zeros((rows+1, cols+1), dtype=np.float)
        val = np.zeros((rows+1, cols+1), dtype=np.float)
        ## decide matrix and path
        for i in xrange(1, rows+1):
            for j in xrange(1, cols+1):
                D = val[i-1][j-1] + S[i-1][j-1]
                H = val[i-1][j]
                if dir[i-1][j] == 1:
                    H = H + gapcost
                V = val[i][j-1]
                if dir[i][j-1] == 1:
                    V = V + gapcost
            
                if D >= H and D >= V:
                    dir[i][j] = 1
                    val[i][j] = D
                else:
                    dir[i][j] = 0
                    val[i][j] = max(V, H)
        #print dir
        #print val
        ## extract the alignment
        (i, j) = (rows, cols); path = []
        while i > 0 and j > 0:
            if dir[i][j] == 1:
                path.append((i - 1, j - 1))
                i -= 1
                j -= 1
            else:
                H = val[i-1][j]
                if dir[i-1][j] == 1:
                    H = H + gapcost
                V = val[i][j-1]
                if dir[i][j-1] == 1:
                    V = V + gapcost
                if V >= H:
                    j -= 1
                else:
                    i -= 1
        path.reverse()
        return path

class Align:
    def __init__(self):
	self.munkres = Munkres()
        self.dp = DP()
      
    def non_sequential(self, coords1, coords2):
	dist_matrix = distance.cdist(coords1, coords2, 'euclidean').astype(np.int32)
	cost_matrix = self.munkres.run_munkres(dist_matrix, 0)
	non_zero = cost_matrix > 0
	return np.column_stack(np.where(non_zero))
	
    def _zeros(self, shape):
        retval = []
        for x in range(shape[0]):
            retval.append([])
        for y in range(shape[1]):
            retval[-1].append(0)
        return retval  

    def sequential(self, coords1, coords2):
        gap_penalty = -0.6
        dist_matrix = distance.cdist(coords1, coords2, 'euclidean').astype(np.float)
        # simlarity like in tmalign
        Lmin = min(coords1.shape[0], coords2.shape[0])
        d0Lmin = (1.24 * ((1 + math.cos(Lmin-1.5)) ** (1 / 3.0))) ** 2
        siml_matrix = 1 / (1 + (dist_matrix * 2)/d0Lmin)
        # inverse of rmsd
        #siml_matrix = 1 - dist_matrix / dist_matrix.max()
        return self.dp.dp2(gap_penalty, siml_matrix)
      
if __name__ == '__main__':
    print 'local test'
    l = Align().sequential(None, None)
    l.reverse()
    print l
    
