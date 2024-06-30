import numpy as np
import PyDealII.Release as dealii

dx , dy , dz= np.array([80, 70, 60, 60, 60, 60, 70, 80]), np.array([80, 70, 60, 60, 60, 60, 70, 80]), np.array([80, 70, 60, 60, 60, 60, 70, 80])
base_value = 80
for i in range(0, 6):
    if i < 2:
        base_value = base_value * 1.2
    elif 2 <= i< 4:
        base_value = base_value * 1.4
    elif 4 <= i:
        base_value= base_value*1.6
    add_value = np.array([base_value])
    dx ,dy , dz = np.concatenate((add_value,dx, add_value)),np.concatenate((add_value,dy, add_value)),np.concatenate((dz, add_value))

p_origin = dealii.Point([-np.sum(dx) / 2.0, -np.sum(dy) / 2.0, 0.0])
p_end = dealii.Point([np.sum(dx) / 2.0, np.sum(dy) / 2.0, np.sum(dz)*1.0])
tria = dealii.Triangulation('3D', dealii.MeshSmoothing.limit_level_difference_at_vertices)
tria.generate_subdivided_steps_hyper_rectangle(
    [dx.tolist(), dy.tolist(), dz.tolist()], p_origin, p_end, False)

for cell in tria.active_cells():
    cell.material_id = 0

tria.save('halfspace.tria')
tria.write('halfspace.vtu', 'vtu')

rhos = [100]
with open('halfspace.rho', 'w') as rhof:
    print('%d' % (len(rhos)), file=rhof)
    for rho in rhos:
        print('%g' % (rho), file=rhof)

source = [0, 0, 0]
sites = [1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 80, 100, 150, 200, 250, 300]
with open('halfspace.emd', 'w') as emdf:
    print('%g %g %g' % (source[0], source[1], source[2]), file=emdf)
    print('%d' % (len(sites)), file=emdf)
    for x in sites:
        print('%g 0 0' % (x), file=emdf)