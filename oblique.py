import sys
import os
import numpy as np
from numpy import sin, cos, exp, pi, e, sqrt

import dataframe

c_light = 137.036

class ObliqueImporter:

    def __init__(self, sysname, angle=0.0):
        title = "%s_ac" % sysname
        self.df = dataframe.DataFrame(title)
        self.theta = angle * (pi / 180.0)
        self._r_tx = (c_light * self.df.dt) / (cos(self.theta) * self.df.dx)
        self._r_yt = (sin(self.theta) * self.df.dy) / (c_light * self.df.dt)

    def close(self):
        return self.df.close()

    def interp_field_on_x(self, iy, iz, r_it):
        if r_it < self.df.mt1:
            idx = (r_it - self.df.mt1) * self._r_tx
            f = self.df.read_field_on_x(iy, iz, self.df.mt1)
            g = np.zeros_like(f)
            i = np.arange(self.df.mx, dtype=float)
            for j in range(3):
                g[j] = np.interp(i, i + idx, f[j])
        else:
            it1, it2 = int(r_it), int(r_it) + 1
            wt1, wt2 = it2 - r_it, r_it - it1
            f1 = self.df.read_field_on_x(iy, iz, it1)
            f2 = self.df.read_field_on_x(iy, iz, it2)
            g = f1 * wt1 + f2 * wt2
        return g

    def predict_field_on_x(self, iy, iz, it, iy_piv):
        r_it = it + (iy_piv - iy) * self._r_yt
        if r_it < self.df.mt2:
            return self.interp_field_on_x(iy_piv, iz, r_it)
        else:
            return np.zeros([3, self.df.mx])

    def predict_field_on_x2(self, iy, iz, it, jy1, jy2):
        fs = np.zeros([3, self.mx])
        for jy_piv in range(jy1, jy2):
            fs += self.predict_field_on_x(iy, iz, it)
        return fs * (1.0 / (jy2 - jy1))

    def __enter__(self):
        return self

    def __exit__(self):
        return self.close()


def predict_bc(ob, title, nt, dt, directory=os.curdir):
    # Default jy-smearing range
    jy1 = int(ob.df.my1 * 0.75 + ob.df.my2 * 0.25)
    jy2 = int(ob.df.my1 * 0.25 + ob.df.my2 * 0.75)
    # Export .info file:
    file_info = os.path.join(directory, '%s.info' % title)
    with open(file_info, 'w') as fh_info:
        sys.stderr.write('# Export file=%s\n' % fh_info.name)
        fh_info.write('%d %d %e\n' % (ob.df.mx1, ob.df.mx2, ob.df.dx))
        fh_info.write('%d %d %e\n' % (1, 2, 0.0))
        fh_info.write('%d %d %e\n' % (1, 1, 0.0))
        fh_info.write('%d %d %e\n' % (1, nt + 1, dt))
    # Export .bin file:
    file_bin = os.path.join(directory, '%s.bin' % title)
    with open(file_bin, 'w') as fh_bin:
        iz = int((ob.df.mz1 + ob.df.mz2) / 2)
        sys.stderr.write('# Export file=%s\n' % fh_bin.name)
        for it in range(1, ob.df.mt2 + 1 + 1):
            jt = (it * dt) / ob.df.dt
            f1 = ob.predict_field_on_x2(ob.df.my1, iz, jt, jy1, jy2)
            f2 = ob.predict_field_on_x2(ob.df.my2, iz, jt, jy1, jy2)
            f1.tofile(fh_bin)
            f2.tofile(fh_bin)
            if it % 1000 == 0:
                sys.stderr.write('# Generated it=%d\n' % it)
