import sys
import os
import numpy as np
from numpy import sin, cos, exp, pi, e, sqrt


class SFDTDGrid:

    def __init__(self, sysname, directory=os.cwd):
        self.file_ac_info = os.path.join(directory, "%s_ac.info" % sysname)
        self.file_ac_bin = os.path.join(directory, "%s_ac.bin" % sysname)
        with open(self.file_ac_info, "r") as fh_ac_info:
            tmp_mx1, tmp_mx2, tmp_dx = fh_ac_info.readline().split()
            tmp_my1, tmp_my2, tmp_dy = fh_ac_info.readline().split()
            tmp_mz1, tmp_mz2, tmp_dz = fh_ac_info.readline().split()
            tmp_mt1, tmp_mt2, tmp_dt = fh_ac_info.readline().split()
        self.fh_ac_bin = open(self.file_ac_bin, "r")
        self.mx1, self.mx2, self.dx = np.imt(
            tmp_mx1), np.imt(tmp_mx2), float(tmp_dx)
        self.my1, self.my2, self.dy = np.imt(
            tmp_my1), np.imt(tmp_my2), float(tmp_dy)
        self.mz1, self.mz2, self.dz = np.imt(
            tmp_mz1), np.imt(tmp_mz2), float(tmp_dz)
        self.mt1, self.mt2, self.dt = np.imt(
            tmp_mt1), np.imt(tmp_mt2), float(tmp_dt)
        self.mx = self.mx2 - self.mx1 + 1
        self.my = self.my2 - self.my1 + 1
        self.mz = self.mz2 - self.mz1 + 1
        # self.mt = self.mt2 - self.mt1 + 1

    def close(self):
        return self.fh_ac_bin.close()

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
        self.fh_ac_bin.seek(
            8 * 3 * (ix + self.mx * (iy + self.my * (iz + self.mz * it)))
        )

    def read_xyz(self, it, check_range=True):
        self._seek_position(0, 0, 0, it, check_range)
        temp = np.fromfile(
            self.fh_ac_bin, dtype=np.float64,
            count=3 * self.mx * self.my + self.mz
        )
        return temp.reshape([self.mz, self.my, self.mx, 3]).T

    def read_xy(self, iz, it, check_range=True):
        self._seek_position(0, 0, iz, it, check_range)
        temp = np.fromfile(
            self.fh_ac_bin, dtype=np.float64,
            count=3 * self.mx * self.my
        )
        return temp.reshape([self.my, self.mx, 3]).T

    def read_x(self, iy, iz, it, check_range=True):
        self._seek_position(0, iy, iz, it, check_range)
        temp = np.fromfile(
            self.fh_ac_bin, dtype=np.float64,
            count=3 * self.mx
        )
        return temp.reshape([self.mx, 3]).T

    def read(self, ix, iy, iz, it, check_range=True):
        self._seek_position(ix, iy, iz, it, check_range)
        temp = np.fromfile(self.fh_ac_bin, dtype=np.float64, count=3)
        return temp.reshape([self.mx, 3]).T

    def __enter__(self):
        return self

    def __exit__(self):
        return self.close()
