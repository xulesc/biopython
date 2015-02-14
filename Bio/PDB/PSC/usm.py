# Copyright (C) 2015, Anuj Sharma (anuj.sharma80@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
This module provides an implementation for calculating similarity between two
protein structures using the universal similarity metric (usm). Description and
theoretical basis for its applicability to protein structure comparison can be 
found in 

"Measuring the similarity of protein structures by means of the universal 
similarity metric", N. Krasnogor and D. A. Pelta, Bioinformatics, Volume 20 
Issue 7, May 2004, Pages 1015-1021.
"""

import numpy as np
from scipy.spatial import distance
import zlib

class USM:
  """
  USM implements a method for pairwise comparison of protein structure. First
  the contact map for each protein structure is calculated. The similarity of
  the protein pair is calculated from the compression ratios of the contact
  maps.
  
  Reference:
  "Measuring the similarity of protein structures by means of the universal 
  similarity metric", N. Krasnogor and D. A. Pelta, Bioinformatics, Volume 20 
  Issue 7, May 2004, Pages 1015-1021.
  """
  def __init__(self, contact_dist = 6.5, max_coords = 1000, metric = 'euclidean'):
    """
    Initialize the USM object. Default values for the contact_dist is set to
    6.5, max_coords is set to 1000 and the metric is set to euclidean. These
    values are set to values described in the original paper.
    """
    self.contact_dist = contact_dist
    self.max_coords = max_coords
    self.metric = metric
    
  def _get_cm_file_string(self, X, Y, n_contacts, n_atoms):
    "Returns the contact map as a string."
    st = '\n'.join(map(lambda x : '%s %s' %(x[0], x[1]), zip(X, Y)))
    return '%d\t# Number of Residues\n%d\t# Number of Contacts\n%s' \
      %(n_atoms, n_contacts, st)
    
  def _compute_distance(self, x, y, xy, yx):
    "Compute the distance given the compression ratios."
    return max(yx - y, xy - x) / max(x, y)
    
  def _in_memory_compress(self, to_zip):
    "Calculate the size of the compressed string in memory."
    return float(len(zlib.compress(to_zip)))
    
  def _make_contact_map(self, coords):
    """
    Returns the contact map of a protein structure given the coordinates of the
    structure. 
    
    First the distances between all residues are calculated. Two residues are
    considered to be in contact if they are less than contact_dist apart and
    the residues are separated by 2 indices. 
    
    The implementation ensures that contact r1-r2 & r2-r1 is reported only once.
    """
    dist = distance.cdist(coords, coords, self.metric)
    idx1, idx2 = np.where(dist <= self.contact_dist)
    fidx = np.where((idx2-idx1) >= 2)
    return [idx1[fidx], idx2[fidx], len(fidx[0]), len(coords)]

  def get_contact_map(self, coords):
    """
    Calculates the contact map for the structure given coordinates of the 
    residues.
    
    Returns 
      - the indices of residues that are in contact
      - string representation for storage to file
    """
    [X, Y, n_contacts, n_atoms] = self._make_contact_map(
      coords[0:min(len(coords), self.max_coords)])
    return (zip(X, Y), self._get_cm_file_string(X, Y, n_contacts, n_atoms))
  
  def dist(self, cm1, cm2):
    """
    Calculates the distance between two proteins given the string representation
    of their contact maps.
    """
    x = self._in_memory_compress(cm1)
    y = self._in_memory_compress(cm2)
    xy = self._in_memory_compress(cm1 + cm2)
    yx = self._in_memory_compress(cm2 + cm1)
    return self._compute_distance(x, y, xy, yx)
