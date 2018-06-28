import sys
import os
import numpy as np
from numpy import sin, cos, exp, pi, e, sqrt

class DataFrame:

    def __init__(self, title, directory=os.curdir):
        self.file_info = os.path.join(directory, '%s.info' % title)
        with open(self.file_info, "r") as fh_info:
            tmp_mx1, tmp_mx2, tmp_dx = fh_info.readline().split()
            tmp_my1, tmp_my2, tmp_dy = fh_info.readline().split()
            tmp_mz1, tmp_mz2, tmp_dz = fh_info.readline().split()
            tmp_mt1, tmp_mt2, tmp_dt = fh_info.readline().split()
        self.mx1, self.mx2 = np.int(tmp_mx1), np.int(tmp_mx2)
        self.my1, self.my2 = np.int(tmp_my1), np.int(tmp_my2)
        self.mz1, self.mz2 = np.int(tmp_mz1), np.int(tmp_mz2)
        self.mt1, self.mt2 = np.int(tmp_mt1), np.int(tmp_mt2)
        self.dx = float(tmp_dx)
        self.dy = float(tmp_dy)
        self.dz = float(tmp_dz)
        self.dt = float(tmp_dt)
        self.mx = self.mx2 - self.mx1 + 1
        self.my = self.my2 - self.my1 + 1
        self.mz = self.mz2 - self.mz1 + 1
        self.mt = self.mt2 - self.mt1 + 1
        self.file_bin = os.path.join(directory, '%s.bin' % title)
        self.fh_bin = open(self.file_bin, "r")

    def close(self):
        return self.fh_bin.close()

    def _seek_position(self, ix, iy, iz, it, check_range):
        if check_range:
            if (ix < self.mx1 or self.mx2 < ix):
                raise IndexError("x-index out of range")
            if (iy < self.my1 or self.my2 < iy):
                raise IndexError("y-index out of range")
            if (iz < self.mz1 or self.mz2 < iz):
                raise IndexError("z-index out of range")
            if (it < self.mt1 or self.mt2 < it):
                raise IndexError("time-index out of range")
        self.fh_bin.seek(
            8 * 3 * (
                (ix - self.mx1) + self.mx * (
                (iy - self.my1) + self.my * (
                (iz - self.mz1) + self.mz * (
                (it - self.mt1)
            ))))
        )

    def read_field_on_xyz(self, it, check_range=True):
        self._seek_position(self.mx1, self.my1, self.mz1, it, check_range)
        temp = np.fromfile(
            self.fh_bin, dtype=np.float64,
            count=3 * self.mx * self.my * self.mz
        )
        print(3 ,self.mx, self.my, self.mz)
        return temp.reshape([self.mz, self.my, self.mx, 3]).T

    def read_field_on_xy(self, iz, it, check_range=True):
        self._seek_position(self.mx1, self.my1, iz, it, check_range)
        temp = np.fromfile(
            self.fh_bin, dtype=np.float64,
            count=3 * self.mx * self.my
        )
        return temp.reshape([self.my, self.mx, 3]).T

    def read_field_on_x(self, iy, iz, it, check_range=True):
        self._seek_position(self.mx1, iy, iz, it, check_range)
        temp = np.fromfile(
            self.fh_bin, dtype=np.float64,
            count=3 * self.mx
        )
        return temp.reshape([self.mx, 3]).T

    def read_field(self, ix, iy, iz, it, check_range=True):
        self._seek_position(ix, iy, iz, it, check_range)
        temp = np.fromfile(self.fh_bin, dtype=np.float64, count=3)
        return temp.reshape([self.mx, 3]).T

    def __enter__(self):
        return self

    def __exit__(self):
        return self.close()
