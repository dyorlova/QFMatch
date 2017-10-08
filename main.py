"""Cluster matching.

Algorightm:
There are 2 data files with clustered 2d points: we call them "left" file
and "right" file in the scope of this project.
User can specify which cluster from the left file to find matches for.
Algorithm will take this cluster, mix it with all the clusters from the
right file. Bin the mix, then separate the mix into pairs 
(left cluster, one right cluster) for each right cluster. Then apply the same
"bin grid" for each of the pairs and calculate dissimilarities using quadratic
form based comparison.
"""


import collections
import datetime
import math
import threading
from matplotlib import pyplot
import matplotlib.lines as lines
import matplotlib.patches as mpatches

import binner
import color_generator
import data_loader
import point
from scipy.stats import mstats
import numpy as np
from scipy.spatial import distance


# SETTINGS.
# TODO: consider moving them to sys.args and let user provide them from CLI.
_LEFT_FILENAME = '/Users/user/file1.txt'
_RIGHT_FILENAME = '/Users/user/file2.txt'

# Minimal bin size for binning the mix.
_BIN_SIZE = 1000

# How many first rows in data files contain bogus data (headers, description 
# etc)
_NUM_FIRST_ROWS_TO_SKIP_IN_THE_DATA_FILES = 2
# Which character is used for next line in data files.
_DATA_FILES_LINE_SEPARATOR = '\r'
_COLUMNS_SEPARATOR_REGEX= r'\s+'
# In which column we have X coordinate in data files.
_DATA_FILES_X_COLUMN = 1
# In which column we have Y coordinate in data files.
_DATA_FILES_Y_COLUMN = 2
# In which column we have cluster id in data files.
_DATA_FILES_CLUSTER_ID_COLUMN = 0

# Whether we want not to show clusters with ids < 0 on plot.
_DO_NOT_SHOW_NEGATIVE_CLUSERS_ON_PLOT = True

# Coefficient which will be multiplied with left cluster standart deviation
# to compare with geometric mean of right clusters.
# If right cluster's gmean is laying within this variable multiplied
# with left cluster standart deviation then right cluster will be taken into
# consideration.
_SIGMA_MULTIPLIER_TO_CONSIDER_RIGHT_CLUSTERS_WITH_GMEAN_WITHIN = 3

# OTHER STRING LITERALS.
_DATASET_TYPE_CUSTOM_ATTRIBUTE_NAME = 'dataset_type'
_RIGHT_DATASET = 'right_dataset'
_LEFT_DATASET = 'left_dataset'


class _BinsCollection(object):
  """Stores list of bins + metadata characterizing the whole collection.

  (I.e. total number of points in all bins).
  """

  def __init__(self):
    # List of binner.Bin objects.
    self._bins = []
    # Total number of points in all bins.
    self._total_num_points = 0

  def AddBin(self, bin_to_add):
    """Adds bin to collection.

    Args:
      bin_to_add: binner.Bin object.
    """
    self._bins.append(bin_to_add)
    self._total_num_points += len(bin_to_add.GetPoints())

  def GetBins(self):
    return self._bins

  def GetBin(self, bin_index):
    return self._bins[bin_index]

  def GetTotalNumPoints(self):    
    return self._total_num_points
  


class _Dissimilarity(object):
  """Structure containing information about dissimilarity between 2 clusters."""

  def __init__(self, left_cluster_id, right_cluster_id, dissimilarity_score):
    self.left_cluster_id = left_cluster_id
    self.right_cluster_id = right_cluster_id
    self.dissimilarity_score = dissimilarity_score


class _Matcher(object):
  """Matches clusters.

  Concepts of "left" and "right" entity correspond to first (left) file and
  second (right) file which we match clusters for.
  Concept of "mix" / "mixed" entity means that it somehow includes information 
  from left and right files (i.e. "Bin" contains points from both files).
  """

  def __init__(self):
    # Dict with cluster id as key and all points (points.Point objects) 
    # related to this cluster as value in LEFT (base) file.
    self._all_left_points_by_cluster_id = {}

    # Dict with cluster id as key and all points (points.Point objects) 
    # related to this cluster as value in RIGHT (matching) file.
    self._all_right_points_by_cluster_id = {}

    # List of bins after binning on mixed dataset.
    self._mix_bins = []

    # Dict with cluster id as key and _BinsCollection object with bins from 
    # this right cluster as value.
    self._right_bin_collection_by_cluster_id = collections.defaultdict(
        _BinsCollection)
    self._right_clusters = []

    # Dict with cluster id as key and _BinsCollection object with bins from 
    # this left cluster as value.
    self._left_bin_collection_by_cluster_id = collections.defaultdict(
        _BinsCollection)
    self._left_clusters = []

    self._dissimilarities = []

    self._matched_pairs = []

    # Max distance between means of bins on mixed set of points.
    self._max_distance_between_bins = None

    # Contains information about which axis and which coordinate space
    # was divided by when doing 'kd-tree building'-like binning.
    self._binning_splitting_coordinates = []

  def ProcessAndReturnDissimilarities(self):
    dt = datetime.datetime.now()
    self._LoadLeft()
    self._LoadRight()
    self._MixAndBin()
    self._SeparateMixedBins()
    self._CalculateMaxDistanceBetweenBins()
    self._CalculateDissimilarities()
    self._Match()
    self._MergeUnmatchedRight()
    self._MergeUnmatchedLeft()

    print 'Took %s' % (datetime.datetime.now() - dt).total_seconds()

    self._DrawGraphs()
    return self._dissimilarities
  
  def _LoadLeft(self):
    # Load left points and mark each of it as 'left point' to be able
    # later to separate them after binning.
    cust_attrs_to_set = {_DATASET_TYPE_CUSTOM_ATTRIBUTE_NAME: _LEFT_DATASET}
    self._all_left_points_by_cluster_id = data_loader.DataLoader(
        _LEFT_FILENAME,
         num_first_rows_to_skip=
         _NUM_FIRST_ROWS_TO_SKIP_IN_THE_DATA_FILES,
         line_separator=_DATA_FILES_LINE_SEPARATOR,
         x_column=_DATA_FILES_X_COLUMN,
         y_column=_DATA_FILES_Y_COLUMN,
         cluster_id_column=_DATA_FILES_CLUSTER_ID_COLUMN,
         cluster_ids_to_exclude={0, -1000},
         columns_separator_regex=_COLUMNS_SEPARATOR_REGEX
    ).LoadAndReturnPointsDividedByClusterId(
        point_custom_attributes=cust_attrs_to_set)
    print 'Left points are loaded. Clusters are %s' % (', ').join(
        [str(s) for s in self._all_left_points_by_cluster_id.iterkeys()])    

  def _LoadRight(self):
    # Load right points and mark each of it as 'right point' to be able
    # later to separate them after binning.
    cust_attrs_to_set = {_DATASET_TYPE_CUSTOM_ATTRIBUTE_NAME: _RIGHT_DATASET}
    self._all_right_points_by_cluster_id = data_loader.DataLoader(
         _RIGHT_FILENAME, 
         num_first_rows_to_skip=
         _NUM_FIRST_ROWS_TO_SKIP_IN_THE_DATA_FILES,
         line_separator=_DATA_FILES_LINE_SEPARATOR,
         x_column=_DATA_FILES_X_COLUMN,
         y_column=_DATA_FILES_Y_COLUMN,
         cluster_id_column=_DATA_FILES_CLUSTER_ID_COLUMN,
         cluster_ids_to_exclude={0, -1000},
         columns_separator_regex=_COLUMNS_SEPARATOR_REGEX
    ).LoadAndReturnPointsDividedByClusterId(
        point_custom_attributes=cust_attrs_to_set)
    print 'Right points are loaded. Clusters are %s' % (', ').join(
        [str(s) for s in self._all_right_points_by_cluster_id.iterkeys()])

  def _MixAndBin(self):
    # Mix left and right files.
    the_mix = []
    for dict_of_points_divided_by_cluster_id in [
        self._all_right_points_by_cluster_id, 
        self._all_left_points_by_cluster_id]:
      for points in dict_of_points_divided_by_cluster_id.itervalues():
        for p in points:
          the_mix.append(p)

    # And bin the medley.
    good_binner = binner.SplittingInHalfBinner(
        the_mix, min_points_per_bin=_BIN_SIZE)
    self._mix_bins = good_binner.GetBins()
    self._binning_splitting_coordinates = good_binner.GetSplittingCoordinates()

  def _SeparateMixedBins(self):
    # Separate the medley keeping the same bin borders as were calculated on the
    # medley.
    for mix_bin in self._mix_bins:
      self._SeparateMixedBin(mix_bin)

  def _SeparateMixedBin(self, mix_bin):
    left_bin_by_cluster_id = {}
    right_bin_by_cluster_id = {}

    # Each bin from mixed set can potentially have points related
    # to each cluster from left and right dataset.
    # Here we will create Bin object corresponding to the bin from mixed set
    # for each left cluster and each right cluster. If there are no points
    # related to particular left or right cluster in this mix bin, then this new
    # Bin object will have no points in it.
    for c_id in self._all_left_points_by_cluster_id.iterkeys():
      left_bin = binner.Bin()
      # In the bin containing only left points keep the mean which was
      # calculated on the mixed bin.
      left_bin.SetFixedMean(mix_bin.GetFixedMean())
      left_bin_by_cluster_id[c_id] = left_bin

    for c_id in self._all_right_points_by_cluster_id.iterkeys():
      right_bin = binner.Bin()
      # In the bin containing only right points keep the mean which was
      # calculated on the mixed bin.
      right_bin.SetFixedMean(mix_bin.GetFixedMean())
      right_bin_by_cluster_id[c_id] = right_bin

    # Do the actual separation of bin with mixed points.
    for cur_point in mix_bin.GetPoints():
      if cur_point.GetCustomAttribute(
          _DATASET_TYPE_CUSTOM_ATTRIBUTE_NAME) == _LEFT_DATASET:
        left_bin_by_cluster_id[cur_point.GetClusterId()].AddPoint(cur_point)
      elif cur_point.GetCustomAttribute(
          _DATASET_TYPE_CUSTOM_ATTRIBUTE_NAME) == _RIGHT_DATASET:
        right_bin_by_cluster_id[cur_point.GetClusterId()].AddPoint(cur_point)
      else:
        raise ValueError(
            'Can not define which dataset point %s belongs to' % cur_point)

    for right_cluster_id, right_bin in right_bin_by_cluster_id.iteritems():
      self._right_bin_collection_by_cluster_id[right_cluster_id].AddBin(
          right_bin)

    for left_cluster_id, left_bin in left_bin_by_cluster_id.iteritems():
      self._left_bin_collection_by_cluster_id[left_cluster_id].AddBin(
          left_bin)

  def _CalculateMaxDistanceBetweenBins(self):
    """Calculate max distance between two farthest-apart mixed bins."""
    total_ops = len(self._mix_bins) * len(self._mix_bins)
    print 'Calculating max distance between bins. Total calculations: %s' % (
        total_ops)
    self._max_distance_between_bins = 0
    
    for bin_i in self._mix_bins:
      for bin_j in self._mix_bins:
        d = _Dist(bin_i.GetFixedMean(), bin_j.GetFixedMean())
        if self._max_distance_between_bins < d:
          self._max_distance_between_bins = d
    print 'Max distance is calculated'

  def _CalculateDissimilarities(self):
    """Calculates dissimilarities between each left and right clusters."""
    print 'Calculating dissimilarities'
    num_bins = len(self._mix_bins)

    for left_cluster_id, left_bin_collection in (self
        ._left_bin_collection_by_cluster_id.iteritems()):
      for right_cluster_id, right_bin_collection in (self
          ._right_bin_collection_by_cluster_id.iteritems()):
        # This operation can be easily parallelized via multiprocessing.
        d = self._CalculateDissimilarityBetweenClusters(
            left_cluster_id, left_bin_collection, right_cluster_id,
            right_bin_collection)
        print 'Left cluster: %s, Right cluster: %s, dissimilarity: %s' % (
            d.left_cluster_id, d.right_cluster_id, d.dissimilarity_score)
        self._dissimilarities.append(d)

    print 'Dissimilarities are calculated'

  def _CalculateDissimilarityBetweenClusters(
      self, left_cluster_id, left_bin_collection, 
      right_cluster_id, right_bin_collection):
    num_bins = len(self._mix_bins)
   
    # Sanity check.
    assert num_bins == len(left_bin_collection.GetBins())
    assert num_bins == len(right_bin_collection.GetBins())

    dissimilarity_score = 0 
    mixed = _MixCollections(left_bin_collection, right_bin_collection)
    max_dist = _CalculateMaxDistanceBetweenBinCollections(mixed, mixed)

    for i in xrange(num_bins):
      for j in xrange(num_bins):
        # Weight of the bin in the left cluster.
        h_i = (len(left_bin_collection.GetBin(i).GetPoints())
               / float(left_bin_collection.GetTotalNumPoints()))
        
        # Weight of the bin in the left cluster.
        h_j = (len(left_bin_collection.GetBin(j).GetPoints())
               / float(left_bin_collection.GetTotalNumPoints()))
        
        f_i = (len(right_bin_collection.GetBin(i).GetPoints())
               / float(right_bin_collection.GetTotalNumPoints()))
        f_j = (len(right_bin_collection.GetBin(j).GetPoints())
               / float(right_bin_collection.GetTotalNumPoints()))
     
        # Using mean calculated on mixed bins.
        i_mean = mixed.GetBin(i).GetFixedMean()
        j_mean = mixed.GetBin(j).GetFixedMean()
        importance_coef = _Dist(j_mean, i_mean) / max_dist

        dissimilarity_score += (
           math.pow((1 - importance_coef), 2) * (h_i - f_i) * (h_j - f_j))

    return _Dissimilarity(
        left_cluster_id, right_cluster_id, dissimilarity_score)

  def _Match(self):
    """Match clusters."""
    # Key is left cluster id, value is _Dissimilarity object containing
    # information about dissimilarity between given left cluster
    # and closest right cluster.
    closest_for_left = {}
    # Key is right cluster id, value is _Dissimilarity object containing
    # information about dissimilarity between given right cluster
    # and closest left cluster.
    closest_for_right = {}

    for diss in self._dissimilarities:
      cur_diss = closest_for_left.get(diss.left_cluster_id)
      if (not cur_diss
          or cur_diss.dissimilarity_score > diss.dissimilarity_score):
        closest_for_left[diss.left_cluster_id] = diss

      cur_diss = closest_for_right.get(diss.right_cluster_id)
      if (not cur_diss 
          or cur_diss.dissimilarity_score > diss.dissimilarity_score):
        closest_for_right[diss.right_cluster_id] = diss

    # Now - find trivial matches for left and right clusters. Leave (in
    # "closest" dicts) only clusters which we were not able to find matches for.
    for left_cluster_id, diss in closest_for_left.items():
      # Make sure that right cluster closest to left cluster X and left 
      # cluster closest to right cluster X match. left.closest = right and 
      # right.closest = left.
      if diss.right_cluster_id in closest_for_right:
        if (closest_for_right[diss.right_cluster_id].left_cluster_id 
            == left_cluster_id):
          print ('Left cluster: %s. Closest right cluster: %s. '
                 'Closest left cluster for the right cluster: %s. Matches.') % (
                     left_cluster_id, diss.right_cluster_id,
                     closest_for_right[diss.right_cluster_id].left_cluster_id)
          self._matched_pairs.append((left_cluster_id, diss.right_cluster_id))
          # We found the pairs for these clusters, delete them from closests
          # dicts.
          del closest_for_right[diss.right_cluster_id]
          del closest_for_left[left_cluster_id]
        else:
          print ('Left cluster: %s. Closest right cluster: %s. '
                 'Closest left cluster for the right cluster: %s. '
                 'Does not match.') % (
                     left_cluster_id, diss.right_cluster_id,
                     closest_for_right[diss.right_cluster_id].left_cluster_id)
      else:
        # Right cluster was already matched to another cluster before.
        # It likely means that there were 2 left clusters which dissimilarity
        # was smallest with the same right cluster.
        print ('Left cluster: %s. Closest right cluster: %s. '
               'Already found the match for right cluster.') % (
                   left_cluster_id, diss.right_cluster_id)

    self._unmatched_right = [(k, v) for k, v in closest_for_right.iteritems()]
    self._unmatched_left = [(k, v) for k, v in closest_for_left.iteritems()]

    print 'Non-matched left clusters: %s' % ', '.join(
        [str(c) for c in closest_for_left.iterkeys()])
    print 'Non-matched right clusters: %s' % ', '.join(
        [str(c) for c in closest_for_right.iterkeys()])
    print 'Initially matched: %s' % self._matched_pairs

  def _MergeUnmatchedRight(self):
    for right_cluster_id, diss in self._unmatched_right:
      merge_with = None
      matched_pair_index = 0
      # Find candidate among right clusters which were matched earlier to merge
      # with currently processed unmatched right cluster.
      for l_id, r_id in self._matched_pairs:
        if (l_id == diss.left_cluster_id 
            # l_id is tuple if left cluster was merged earlier and now it
            # contains ids of merged clusters.
            or (isinstance(l_id, tuple) and diss.left_cluster_id in l_id)
            or (isinstance(r_id, tuple) and diss.right_cluster_id in r_id)):
          merge_with = r_id
          break
        matched_pair_index += 1

      if not merge_with:
        continue

      # X is currently processed right cluster and Y is left cluster closest
      # to it identified on previous step. Z is right cluster which we will
      # merge with X.
      # Then this var will contain dissimilarity between Y and Z.
      merge_with_and_left_diss = None
      for d in self._dissimilarities:
        if (d.left_cluster_id == diss.left_cluster_id 
            and d.right_cluster_id == merge_with):
          merge_with_and_left_diss = d
          break

      bin_collection1 = self._right_bin_collection_by_cluster_id[
          right_cluster_id]
      bin_collection2 = self._right_bin_collection_by_cluster_id[
          merge_with]

      num_bins = len(self._mix_bins)
      assert num_bins == len(bin_collection1.GetBins())
      assert num_bins == len(bin_collection2.GetBins())

      print ('Will try to do right merge. Right cluster: %s, '
             'merging candidate: %s. Left cluster: %s') % (
                 merge_with, right_cluster_id, diss.left_cluster_id)

      # Merge right clusters.
      right_bin_collection = _MixCollections(bin_collection1, bin_collection2)

      # If we merge right cluster with already merged clusters (rather than
      # single cluster).
      if isinstance(merge_with, tuple):
        merged_cluster_id = tuple(list(merge_with) + [right_cluster_id])
      else:
        merged_cluster_id = (merge_with, right_cluster_id)

      left_bin_collection = self._left_bin_collection_by_cluster_id[
          diss.left_cluster_id]

      # Dissimilarity between merged clusters and left cluster.
      new_diss = self._CalculateDissimilarityBetweenClusters(
          diss.left_cluster_id, left_bin_collection, 
          merged_cluster_id, right_bin_collection)
      if (new_diss.dissimilarity_score 
          < merge_with_and_left_diss.dissimilarity_score):
        # If merged right clusters result has better dissimilarity score with
        # left clusters then consider the merged clusters as single cluster.
        print (
            'Merged. Right cluster: %s, merging candidate: %s. Left cluster: '
            '%s\nOld score: %s\nNew score: %s') % (
                merge_with, right_cluster_id, l_id, 
                merge_with_and_left_diss.dissimilarity_score, 
               new_diss.dissimilarity_score)
        # Delete right-left cluster pair from result since we will replace it
        # with merged_right-left cluster dissimilarity.
        del self._matched_pairs[matched_pair_index]
        self._matched_pairs.append((l_id, merged_cluster_id))
        self._right_bin_collection_by_cluster_id[
            merged_cluster_id] = right_bin_collection
        self._dissimilarities.append(new_diss)
      else:
        # Keep status quo - merged clusters -> left cluster dissimilarity
        # is higher than single right -> left dissimilarity.
        print ('NOT Merged. Right cluster: %s, merging candidate: %s. '
               'Left cluster: %s\nOld score: %s\nNew score: %s') % (
                   merge_with, right_cluster_id, l_id, 
                   merge_with_and_left_diss.dissimilarity_score, 
                   new_diss.dissimilarity_score)

  def _MergeUnmatchedLeft(self):
    for left_cluster_id, diss in self._unmatched_left:
      merge_with = None
      matched_pair_index = 0
      for l_id, r_id in self._matched_pairs:
        if (r_id == diss.right_cluster_id
            or (isinstance(l_id, tuple) and diss.left_cluster_id in l_id)
            or (isinstance(r_id, tuple) and diss.right_cluster_id in r_id)):
          merge_with = l_id
          break
        matched_pair_index += 1

      if not merge_with:
        continue

      merge_with_and_right_diss = None
      for d in self._dissimilarities:
        if (d.right_cluster_id == diss.right_cluster_id and 
            d.left_cluster_id == merge_with):
          merge_with_and_right_diss = d
          break

      print ('Will try to do left merge. Left cluster: %s, merging candidate: '
             '%s. Right cluster: %s') % (merge_with, left_cluster_id, 
                                         diss.right_cluster_id)
      
      bin_collection1 = self._left_bin_collection_by_cluster_id[left_cluster_id]
      bin_collection2 = self._left_bin_collection_by_cluster_id[merge_with]

      num_bins = len(self._mix_bins)
      assert num_bins == len(bin_collection1.GetBins())
      assert num_bins == len(bin_collection2.GetBins())

      left_bin_collection = _MixCollections(bin_collection1, bin_collection2)

      if isinstance(merge_with, tuple):
        merged_cluster_id = tuple(list(merge_with) + [left_cluster_id])
      else:
        merged_cluster_id = (merge_with, left_cluster_id)

      right_bin_collection = self._right_bin_collection_by_cluster_id[
          diss.right_cluster_id]

      new_diss = self._CalculateDissimilarityBetweenClusters(
          merged_cluster_id, left_bin_collection, 
          diss.right_cluster_id, right_bin_collection)
      if (new_diss.dissimilarity_score 
          < merge_with_and_right_diss.dissimilarity_score):
        print ('Merged. Left cluster: %s, merging candidate: %s. '
               'Right cluster: %s\nOld score: %s\n'
               'New score: %s') % (
                   merge_with, left_cluster_id, diss.right_cluster_id, 
                   merge_with_and_right_diss.dissimilarity_score, 
                   new_diss.dissimilarity_score)
        del self._matched_pairs[matched_pair_index]
        self._left_bin_collection_by_cluster_id[
            merged_cluster_id] = left_bin_collection
        self._matched_pairs.append((merged_cluster_id, r_id))
        self._dissimilarities.append(new_diss)
      else:
        print ('NOT Merged. Left cluster: %s, merging candidate: %s. '
               'Right cluster: %s\nOld score: %s\n'
               'New score: %s') % (
                   merge_with, left_cluster_id, diss.right_cluster_id, 
                   merge_with_and_right_diss.dissimilarity_score, 
                   new_diss.dissimilarity_score)


  def _DrawGraphs(self):
    print 'Matched: %s' % self._matched_pairs
    left_xs = []
    left_ys = []
    left_colors = []
    left_sizes = []

    right_xs = []
    right_ys = []
    right_colors = []
    right_sizes = []
    
    chunk_ids_for_color_generation = [i for i in xrange(len(self._matched_pairs))]
    color_gen = color_generator.ColorGenerator(
        chunk_ids_for_color_generation, 
        exclude_colors=[
            color_generator.KELLY_COLORS_BY_COLOR_NAME[color_generator.STRONG_BLUE]])
    colors_by_left_cluster_id = {}
    colors_by_right_cluster_id = {}

    for i, match in enumerate(self._matched_pairs):
      left_color = color_gen.GetColor(i)
      if isinstance(match[0], tuple):
        for c_id in match[0]:
          colors_by_left_cluster_id[c_id] = left_color
      else:
        colors_by_left_cluster_id[match[0]] = left_color

      right_color = color_gen.GetColor(i)
      if isinstance(match[1], tuple):
        for c_id in match[1]:
          colors_by_right_cluster_id[c_id] = right_color
      else:
        colors_by_right_cluster_id[match[1]] = right_color

    for cluster_id, points in self._all_left_points_by_cluster_id.iteritems():
      for cur_point in points:
        if (_DO_NOT_SHOW_NEGATIVE_CLUSERS_ON_PLOT 
            and cur_point.GetClusterId() < 0):
          continue
        elif cur_point.GetClusterId() == 0:
          continue
        else:
          x, y = cur_point.GetCoordinates()
          left_xs.append(x)
          left_ys.append(y)
          if cur_point.GetClusterId() in colors_by_left_cluster_id:
            left_colors.append(colors_by_left_cluster_id[cur_point.GetClusterId()])
            left_sizes.append(20)
          else:
            left_colors.append(
                color_generator.GetKellyColor(color_generator.STRONG_BLUE))
            left_sizes.append(7)

    for cluster_id, points in self._all_right_points_by_cluster_id.iteritems():
      for cur_point in points:
        if (_DO_NOT_SHOW_NEGATIVE_CLUSERS_ON_PLOT
            and cur_point.GetClusterId() < 0):
          continue
        elif cur_point.GetClusterId() == 0:
          continue          
        else:
          x, y = cur_point.GetCoordinates()
          right_xs.append(x)
          right_ys.append(y)
          if cur_point.GetClusterId() in colors_by_right_cluster_id:
            right_colors.append(colors_by_right_cluster_id[cur_point.GetClusterId()])
            right_sizes.append(20)
          else:
            right_colors.append(
                color_generator.GetKellyColor(color_generator.STRONG_BLUE))
            right_sizes.append(7)

    fig = pyplot.figure()

    left_ax = fig.add_subplot(121)
    
    left_patches = []
    for c_id, color in colors_by_left_cluster_id.iteritems():
      if c_id == 0:
        continue
      elif c_id < 0 and _DO_NOT_SHOW_NEGATIVE_CLUSERS_ON_PLOT:
        continue
      else:
        patch = mpatches.Patch(color=color, label=str(c_id))
        left_patches.append(patch)

    right_patches = []
    for c_id, color in colors_by_right_cluster_id.iteritems():
      if c_id == 0:
        continue
      elif c_id < 0 and _DO_NOT_SHOW_NEGATIVE_CLUSERS_ON_PLOT:
        continue
      else:
        patch = mpatches.Patch(color=color, label=str(c_id))
        right_patches.append(patch)

    left_ax.scatter(left_xs, left_ys, c=left_colors, s=left_sizes)
    left_ax.legend(handles=left_patches, loc=4)

    right_ax = fig.add_subplot(122)
    right_ax.scatter(right_xs, right_ys, c=right_colors, s=right_sizes)
    right_ax.legend(handles=right_patches, loc=4)

    pyplot.show()

  def _CalculateWithinHowManySigmaRightClusterGmeanLays(
      self, right_cluster_gmean):
    assert self._left_cluster_gmean is not None
    assert self._left_cluster_std is not None
    return np.absolute(
        right_cluster_gmean - self._left_cluster_gmean) / self._left_cluster_std    

  def _PlotPoints(self, plot_left=True, plot_right=True, with_bins=False):
    """Can be used to plot points and their binning borders.

    Useless in calculations.
    """
    xs = []
    ys = []
    colors = []

    fig = pyplot.figure()
    ax = fig.add_subplot(111)

    if plot_left:
      for _, ps in self._all_left_points_by_cluster_id.iteritems():
        for p in ps:
          xs.append(p.GetCoordinate(0))
          ys.append(p.GetCoordinate(1))
          colors.append('blue')

    if plot_right:
      for _, ps in self._all_right_points_by_cluster_id.iteritems():
        for p in ps:
          xs.append(p.GetCoordinate(0))
          ys.append(p.GetCoordinate(1))
          colors.append('red')

    if with_bins:
      for splitting_coordinate in self._binning_splitting_coordinates:
        ax.add_line(lines.Line2D(
          [splitting_coordinate[0][0], splitting_coordinate[0][1]], 
          [splitting_coordinate[1][0], splitting_coordinate[1][1]], 
          color='black'))

    ax.scatter(xs, ys, c=colors)
    # Change these values according to points min and max coordinates.
    pyplot.ylim([0,1])  
    pyplot.xlim([0,1])
    pyplot.show()       


def _CalculateMaxDistanceBetweenBinCollections(
    bin_collection1, bin_collection2):  
  """Calculate max distance between means of bins in two collections."""
  max_distance_between_bins = 0

  assert len(bin_collection1.GetBins()) == len(bin_collection2.GetBins())
  if len(bin_collection1.GetBins()) == 1:
    raise ValueError('Can not calculate distance between bin collections.')
  
  for bin_i in bin_collection1.GetBins():
    for bin_j in bin_collection2.GetBins():
      if bin_i.GetPoints() and bin_j.GetPoints():
        d = _Dist(bin_i.GetFixedMean(), bin_j.GetFixedMean())
        if max_distance_between_bins < d:
          max_distance_between_bins = d 
  return max_distance_between_bins


def _MixCollections(bc1, bc2):
  """Mixes two _BinCollection objects."""
  new = _BinsCollection()
  for i, b1 in enumerate(bc1.GetBins()):
    b2 = bc2.GetBin(i)
    new_bin = binner.Bin()
    new_bin.SetFixedMean(b1.GetFixedMean())
    assert all(b1.GetFixedMean() == b2.GetFixedMean())
    for p in b1.GetPoints():
      new_bin.AddPoint(p)
    for p in b2.GetPoints():
      new_bin.AddPoint(p)
    new.AddBin(new_bin)
  return new


def _Dist(coordinates1, coordinates2):
  """Euclidean distance between N-dimensional points."""
  if len(coordinates2) == 1:
    return math.abs(coordinates2[0] - coordinates1[0])
  elif len(coordinates2) == 2:
    return math.sqrt(
        math.pow(
            coordinates2[0] - coordinates1[0], 2) 
        + math.pow(coordinates2[1] - coordinates1[1], 2)
    )
  else:
    return distance.euclidean(coordinates2, coordinates1)


def main(unused_argv): 
  _Matcher().ProcessAndReturnDissimilarities()


if __name__ == '__main__':
  main(None)
