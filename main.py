'''
    Copyright (C) <2020> <Author: Weikang Tang>
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''
import sys
import reqdsk
import poincareRK4

if __name__ == '__main__':
    gf=reqdsk.Gfile(sys.argv[1])
    gf.g2h5()
#    gf.pltpsi()
#    gf.pltf()
    ms=poincareRK4.Mesh(gf.Raxis, gf.Zaxis, 0.99*max(gf.bdr), int(sys.argv[2]), gf.funcbr, gf.funcbz)
    ms.crt()
    ms.itp(gf.Psi_bound, gf.Psi_axis, gf.funcpsi, gf.qin, gf.pin, gf.fin, gf.ff, gf.fp)
