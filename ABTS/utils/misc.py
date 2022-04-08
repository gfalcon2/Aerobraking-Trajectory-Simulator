def find_neighbour(point, list):
    dist_list = [abs(point-p) for p in list]
    temp = sorted(set(dist_list))
    index_1 = dist_list.index(temp[0])
    index_2 = dist_list.index(temp[1])
    # index_sorted = sorted(set(range(len(dist_list)),key=dist_list.__getitem__))

    temp = (point-list[index_1]), (point-list[index_2])
    temp_index = [index_1,index_2]
    index_sorted = sorted(range(len(temp)),key=temp.__getitem__)

    return temp_index[index_sorted[0]],temp_index[index_sorted[1]], (temp[index_sorted[0]]), (temp[index_sorted[1]])


def new_periapsis(m,r,v,args):
    import numpy as np
    Energy = (np.linalg.norm(v) ** 2) * 0.5 - m.planet.mu / (np.linalg.norm(r))
    a = - m.planet.mu / (2 * Energy)
    h = np.cross(r, v)
    h += 0.
    e = (1 + (2.0 * Energy * (np.inner(h, h)) / (m.planet.mu) ** 2)) ** 0.5
    r_p = a * (1.0 - e)
    if args.print_res:
        print('NEW PERIAPSIS VACUUM ALTITUDE:', (r_p-m.planet.Rp_e) * 1e-3, ' km')