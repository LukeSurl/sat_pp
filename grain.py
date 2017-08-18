#!/usr/bin/env python
# encoding: utf-8
"""
Author: Nick Bearson, nickb@ssec.wisc.edu
Copyright (c) 2014 University of Wisconsin SSEC. All rights reserved.

Get new leap second files @ ftp://time.nist.gov/pub/leap-seconds.list
"""

from datetime import datetime, timedelta

# Useful constants...
VIIRS_EPOCH = datetime(1958, 1, 1)
MODIS_EPOCH = datetime(1993, 1, 1)
UNIX_EPOCH  = datetime(1970, 1, 1)

NTP_EPOCH = datetime(1900, 1, 1)

DEFAULT_EPOCH = NTP_EPOCH

# Use the leap-seconds file included with this package if none is specified
from pkg_resources import resource_filename
DEFAULT_LEAP_SECONDS = resource_filename(__name__, 'leap-seconds')

# ----------------------------------------------------------------

class Grain(object):
  """An object for parsing and utilizing the information contained in a leap second file.
  
  Initialize with a file object, ie: the result of open(...)"""
  def __init__(self, leap_second_file=None):
    if leap_second_file is None:
      leap_second_file = open(DEFAULT_LEAP_SECONDS)
    leap_times = []
    offsets = []
    for line in leap_second_file:
      li=line.strip()
      if not li.startswith('#'):
        pieces = li.split()
        leap_time = NTP_EPOCH + timedelta(seconds=int(pieces[0]))
        leap_times.append(leap_time)
        offset = int(pieces[1])
        offsets.append(offset)

    # Convert our offsets from leapseconds-from-beginning to the change at each instant - it'll make
    # things easier to deal with later!
    offsets.insert(0,0) # add 0 to the start so our math works in the next line
    offsets = [j-i for i, j in zip(offsets[:-1], offsets[1:])]
    self.leaps = zip(leap_times, offsets)



  def _leaps_between(self, date1, date2):
    """
    Counts the number of leap seconds that have occurred between two datetimes
    """
    if date1 > date2:
      raise RuntimeError, "date1 > date2"
    between_times = [i for i in self.leaps if date1 <= i[0] and i[0] <= date2] # FIXME: should these be > or >=?
    offset = sum(leap[1] for leap in between_times) # sum all the offsets in self.leaps
    return offset


  def utc2tai(self, utc, epoch=DEFAULT_EPOCH):
    """
    Takes datetime object (utc) and returns TAI seconds since given epoch. 
    """
    offset = self._leaps_between(epoch, utc)
    tai = utc - epoch
    seconds_since_epoch = (tai.days * (24 * 60 * 60)) + tai.seconds + offset
    return seconds_since_epoch


  def tai2utc(self, seconds_since_epoch, epoch=DEFAULT_EPOCH):
    """
    Takes TAI seconds since given epoch and returns a datetime.
    """
    td_sse = timedelta(seconds=seconds_since_epoch)
    utc_unadjusted = td_sse + epoch

    offset = self._leaps_between(epoch, utc_unadjusted)
    td_offset = timedelta(seconds=offset)
    utc = utc_unadjusted - td_offset
    return utc


def test():
  ls_file = open('leap-seconds')
  g = Grain(ls_file)
  now = datetime.utcnow()
  print "test TAI times, seconds since NTP_EPOCH: ", NTP_EPOCH
  print "VIIRS EPOCH: ", g.utc2tai(VIIRS_EPOCH), VIIRS_EPOCH
  print "MODIS EPOCH: ", g.utc2tai(MODIS_EPOCH), MODIS_EPOCH
  print "NOW: ", g.utc2tai(now), now


if __name__ == '__main__':
  test()
