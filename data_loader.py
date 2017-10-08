import collections
import re
import point


class DataLoader(object):
  """Loads data from file to list of 2D point.Point objects."""

  def __init__(self, filename, num_first_rows_to_skip=2, line_separator='\r', 
               x_column=0, y_column=1, cluster_id_column=2, cluster_ids_to_exclude=None,
               columns_separator_regex=r'\s'):
    self._filename = filename
    self._num_first_rows_to_skip = num_first_rows_to_skip
    self._line_separator = line_separator
    self._columns_separator_regex = columns_separator_regex
    self._x_column = x_column
    self._y_column = y_column
    self._cluster_id_column = cluster_id_column
    self._cluster_ids_to_exclude = cluster_ids_to_exclude or set()
    assert self._x_column != self._y_column != self._cluster_id_column
    # We support 2D only - thus file should contain cluster id, x, y.
    assert self._x_column <= 2
    assert self._x_column >= 0
    assert self._y_column <= 2
    assert self._y_column >= 0
    assert self._cluster_id_column <= 2
    assert self._cluster_id_column >= 0

  def LoadAndReturnPoints(self, point_custom_attributes=None):
    return list(self._OpenFileAndYieldPoints(
        point_custom_attributes=point_custom_attributes))

  def LoadAndReturnPointsDividedByClusterId(self, point_custom_attributes=None):
    points_by_cluster_id = collections.defaultdict(list)
    for point in self._OpenFileAndYieldPoints(
        point_custom_attributes=point_custom_attributes):      
      points_by_cluster_id[point.GetClusterId()].append(point)
    return points_by_cluster_id

  def _OpenFileAndYieldPoints(self, point_custom_attributes=None):
    with open(self._filename, 'r') as file_descr:
      for i, row in enumerate(file_descr.read().split(self._line_separator)): 
        row = row.strip()
        if i >= self._num_first_rows_to_skip and row: 
          try:
            data_list = [float(s) for s in re.split(self._columns_separator_regex, row)]
            x = data_list[self._x_column]
            y = data_list[self._y_column]
            cluster_id = data_list[self._cluster_id_column]
            if cluster_id in self._cluster_ids_to_exclude:
              continue
          except (ValueError, TypeError):
            print 'Failed on processing row "%s"' % row
            raise
          else:
            cur_point = point.Point((x, y), cluster_id)
            if point_custom_attributes:
              for key, value in point_custom_attributes.iteritems():
                cur_point.SetCustomAttribute(key, value)
            yield cur_point