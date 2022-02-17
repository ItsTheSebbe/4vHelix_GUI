import re
import numpy as np

def open_rpoly(rpoly_file):

    polyFile = open(rpoly_file, 'r')

    # Matrix 'data' stores helix coordinates
    data = []
    # fwd_helix_connections and rev_helix_connections: number rows is amount of helices,
    # every row stores oligo connections
    fwd_helix_connections = []
    rev_helix_connections = []
    count = 0

    try:
        for line in polyFile:
            if line.startswith('hb'):
                data.insert(count, line.split(' '))
                count += 1
            elif line.startswith('c'):
                if 'f3' not in line:
                    # rev_helix_connections.append(line.split(' '))
                    rev_helix_connections.append([int(re.search('c helix_(.+?) ',
                                                                line).group(1)),
                                                  int(re.search('\' helix_(.+?) ', line).group(1))])
                else:
                    fwd_helix_connections.append([int(re.search('c helix_(.+?) ',
                                                                line).group(1)),
                                                  int(re.search('\' helix_(.+?) ', line).group(1))])
    except Exception:
        print('Failed to read the file' )

    return data, fwd_helix_connections, rev_helix_connections

def open_ntrail(ntrail_file):

    with open(ntrail_file, 'r') as f:
        content= f.read()
        n_trail_list=content.split()
    n_trail_list = [int(n) for n in n_trail_list]
    return n_trail_list

def open_ply(ply_file):

    number_vertices = 0
    i=0

    with open(ply_file, 'r') as f:
        content = f.read()
        content=re.split('\n',content)
    for line in content:
        if "element vertex" in line:
            number_vertices = [int(s) for s in line.split() if s.isdigit()]
        if "element face" in line:
            number_face = [int(s) for s in line.split() if s.isdigit()]
        if "end_header" in line:
            i+=1
            break
        else:
            i+=1

    faces_list = []
    vertices_list = []

    for line in content[i: i+number_vertices[0]]:
        line = line.split()
        line = [float(x) for x in line]
        vertices_list.append(line)

    for line in content[(i+number_vertices[0]):-1]:
        line = [int(s) for s in line.split() if s.isdigit()]

        faces_list.append((line))
        if int(line[0]) != 3:
            #print(line)
            print("Presence of non triangulated faces detected. The script might not work properly.")

    return number_vertices[0], np.array(vertices_list), number_face[0], np.array(faces_list)

def move_along_vector(point, vector,length):  # rpoly file contains center coordinate of helix, "generate" needs end coordiates of helix:
    move_distance = float(
        length) * 0.4 / 2.0  # 0.4 is the length of a base pair in oxDNA units, move half the helixlength down
    return [point[0] - move_distance * vector[0], point[1] - move_distance * vector[1],
            point[2] - move_distance * vector[2]]