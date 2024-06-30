import numpy as np
import PyDealII.Release as dealii

dx = np.array([100, 80, 70, 60, 50, 40, 30, 20, 18, 16, 14, 12, 10, 10, 10, 10, 10, 
               10, 10, 10, 10, 10, 12, 14, 16, 18, 20, 30, 40, 50, 60, 70, 80, 100])
dy = np.array([100, 80, 70, 60, 50, 40, 30, 20, 18, 14, 12, 10, 10, 12, 14, 18, 20, 30, 40, 50, 60, 70, 80, 100])
dz = np.array([10, 12, 14, 16, 18, 20, 30, 40, 50, 60, 70, 80, 100])

mid_position_control = float(np.sum(dx[0:12])) 

p_origin = dealii.Point([-mid_position_control, -mid_position_control, 0.0])
p_end = dealii.Point([mid_position_control, mid_position_control, np.sum(dz)*1.0])

tria = dealii.Triangulation('3D', dealii.MeshSmoothing.limit_level_difference_at_vertices)
tria.generate_subdivided_steps_hyper_rectangle(
    [dx.tolist(), dy.tolist(), dz.tolist()], p_origin, p_end, False)

for cell in tria.active_cells():
    cell.material_id = 0

rhos = [100]
with open('halfspace.rho', 'w') as rhof:
    print('%d' % (len(rhos)), file=rhof)
    for rho in rhos:
        print('%g' % (rho), file=rhof)

source = [0, 0, 0]
sites = [ 1, 2, 4, 6, 8, 10, 20, 40, 60, 80, 100]

tria.save('halfspace.tria')
tria.write('halfspace.vtu', 'vtu') 
with open('halfspace.emd', 'w') as emdf:
    print('%g %g %g' % (source[0], source[1], source[2]), file=emdf)
    print('%d' % (len(sites)), file=emdf)
    for x in sites:
        print('%g 0 0' % (x), file=emdf)