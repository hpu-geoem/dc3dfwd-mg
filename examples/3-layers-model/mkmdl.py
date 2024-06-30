import numpy as np
import PyDealII.Release as dealii

dx = np.array([150, 100, 150])
dy = np.array([150, 100, 150])
dz = np.array([100, 150, 200])
p_origin = dealii.Point([-np.sum(dx) / 2.0, -np.sum(dy) / 2.0, 0.0])
p_end = dealii.Point([np.sum(dx) / 2.0, np.sum(dy) / 2.0, np.sum(dz)*1.0])

tria = dealii.Triangulation('3D', dealii.MeshSmoothing.limit_level_difference_at_vertices)
tria.generate_subdivided_steps_hyper_rectangle(
    [dx.tolist(), dy.tolist(), dz.tolist()], p_origin, p_end, False)

for cell in tria.active_cells():
    c = cell.center().to_list()
    if  c[2]> 300:
        cell.material_id = 2
    if  c[2]<300 and c[2]> 100:
        cell.material_id = 1
    if  c[2]<100 and c[2]> 0:
        cell.material_id = 0

rhos = [100, 1, 1000]
with open('3-layers-model.rho', 'w') as rhof:
    print('%d' % (len(rhos)), file=rhof)
    for rho in rhos:
        print('%g' % (rho), file=rhof)

source = [0, 0, 0]
sites = [ 1, 2, 5, 10, 20, 30, 40, 50, 60, 80]

tria.save('3-layers-model.tria')
tria.write('3-layers-model.vtu', 'vtu')
with open('3-layers-model.emd', 'w') as emdf:
    print('%g %g %g' % (source[0], source[1], source[2]), file=emdf)
    print('%d' % (len(sites)), file=emdf)
    for x in sites:
        print('%g 0 0' % (x), file=emdf)