import numpy as np
import PyDealII.Release as dealii

streth_factor = 2.6

dx, dy, dz = np.array([50, 100, 100, 50]), np.array([50, 100, 100, 50]), np.array([50, 100, 100, 50])

for _ in range(5):
    x, y, z = dx[-1] * streth_factor, dy[-1] * streth_factor, dz[-1] * streth_factor
    dx, dy, dz = np.append(dx, x), np.append(dy, y), np.append(dz, z)

dxt, dyt, dzt = dx, dy, dz

if streth_factor > 1.9:
    dx, dy, dz = np.sort(dx)[::-1], np.sort(dy)[::-1], np.sort(dz)[::-1]
    dx, dy, dz = dx[:-4], dy[:-4], dz[:-4]
else:
    dx, dy, dz = np.sort(dx[3:])[::-1], np.sort(dy[3:])[::-1], np.sort(dz[3:])[::-1]
    dx, dy, dz = dx[:-1], dy[:-1], dz[:-1]

dx, dy, dz = np.append(dx, dxt), np.append(dy, dyt), np.append(dz, dzt)


p_origin = dealii.Point([-np.sum(dx) / 2.0, -np.sum(dy) / 2.0, 0.0])
p_end = dealii.Point([np.sum(dx) / 2.0, np.sum(dy) / 2.0, np.sum(dz)*1.0])

tria = dealii.Triangulation('3D', dealii.MeshSmoothing.limit_level_difference_at_vertices)
tria.generate_subdivided_steps_hyper_rectangle(
    [dx.tolist(), dy.tolist(), dz.tolist()], p_origin, p_end, False)

for cell in tria.active_cells():
    cell.material_id = 0
    c = cell.center().to_list()
    if abs(c[0])<100 and abs(c[1])<100 and abs(c[2])<(np.sum(dz)/2+100) and abs(c[2])>(np.sum(dz)/2-100):
        cell.refine_flag = 'isotropic'

tria.execute_coarsening_and_refinement()

rhos = [100]
with open('halfspace.rho', 'w') as rhof:
    print('%d' % (len(rhos)), file=rhof)
    for rho in rhos:
        print('%g' % (rho), file=rhof)

source = [0, 0, 0]
sites = [0]

tria.save('halfspace.tria')
tria.write('halfspace.vtu', 'vtu')
with open('halfspace.emd', 'w') as emdf:
    print('%g %g %g' % (source[0], source[1], source[2]), file=emdf)
    print('%d' % (len(sites)), file=emdf)
    for x in sites:
        print('%g 0 0' % (x), file=emdf)