import sys
import os
import numpy as np
from numpy import sin, cos, exp, pi, e, sqrt

import dataframe

c_light = 137.036


class ObliqueData:

    def __init__(self, title, angle=0.0):
        self.df = dataframe.DataFrame(title)
        self.theta = angle * (pi / 180.0)
        self._r_tx = (c_light * self.fd.dt) / (cos(self.theta) * self.fd.dx)
        self._r_yt = (sin(self.theta) * self.fd.dy) / (c_light * self.fd.dt)

    def close(self):
        return self.df.close()

    def interp_field_on_x(self, iy, iz, r_it):
        if r_it < self.fd.mt1:
            idx = (r_it - self.fd.mt1) * self._r_tx
            f = self.fd.read_field_on_x(iy, iz, self.fd.mt1)
            g = np.zeros_like(f)
            i = np.arange(self.fd.mx, dtype=float)
            for j in range(3):
                g[j] = np.interp(i, i + idx, f[j])
        else:
            it1, it2 = int(r_it), int(r_it) + 1
            wt1, wt2 = it2 - r_it, r_it - it1
            f1 = self.fd.read_field_on_x(iy, iz, it1)
            f2 = self.fd.read_field_on_x(iy, iz, it2)
            g = f1 * wt1 + f2 * wt2
        return g

    def predict_field_on_x(self, iy, iz, it, iy_piv):
        r_it = it + (iy_piv - iy) * self._r_yt
        if r_it < self.fd.mt2:
            return self.interp_field_on_x(iy_piv, iz, r_it)
        else:
            return np.zeros([3, self.fd.mx])

    def predict_field_on_x2(self, iy, iz, it, jy1, jy2):
        fs = np.zeros([3, self.mx])
        for jy_piv in range(jy1, jy2):
            fs += self.predict_field_on_x(iy, iz, it)
        return fs * (1.0 / (jy2 - jy1))

    def predict_bc(self, title, jylim=None, directory=os.cwd):
        if jylim is None:
            jy1 = int(self.fd.my1 * 0.75 + self.fd.my2 * 0.25)
            jy2 = int(self.fd.my1 * 0.25 + self.fd.my2 * 0.75)
        else:
            jy1, jy2 = jylim
        # Export .info file:
        file_info = os.path.join(directory, '%s.info' % title)
        with open(file_info, 'w') as fh_info:
            sys.stderr.write('# Export file=%s\n' % fh_info.name)
            fh_info.write('%d %d %e' % (self.fd.mx1, self.fd.mx2, self.fd.dx))
            fh_info.write('%d %d %e' % (1, 2, 0.0))
            fh_info.write('%d %d %e' % (1, 1, 0.0))
            fh_info.write('%d %d %e' % (1, self.fd.mt2 + 1, self.fd.dt))
        # Export .bin file:
        file_bin = os.path.join(directory, '%s.bin' % title)
        with open(file_bin, 'w') as fh_bin:
            iz = int((self.fd.mz1 + self.fd.mz2) / 2)
            sys.stderr.write('# Export file=%s\n' % fh_bin.name)
            for it in range(1, self.fd.mt2 + 1 + 1):
                f1 = self.predict_field_on_x2(self.fd.my1, iz, it, jy1, jy2)
                f2 = self.predict_field_on_x2(self.fd.my2, iz, it, jy1, jy2)
                f1.tofile(fh_bin)
                f2.tofile(fh_bin)
                if it % 1000 == 0:
                    sys.stderr.write('# Generated it=%d\n' % it)

    def __enter__(self):
        return self

    def __exit__(self):
        return self.close()
