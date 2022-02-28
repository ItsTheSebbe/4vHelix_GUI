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

        loaded_correctly = True
        return data, fwd_helix_connections, rev_helix_connections, loaded_correctly

    except Exception:
        loaded_correctly = False
        print('Failed to read the file')
        return data, fwd_helix_connections, rev_helix_connections, loaded_correctly


def open_ntrail(ntrail_file):

    with open(ntrail_file, 'r') as f:
        content = f.read()
        n_trail_list = content.split()
    n_trail_list = [int(n) for n in n_trail_list]
    return n_trail_list


def load_ply_data(ply_file):
    vertices_list = []
    faces_list = []
    with open(ply_file, 'r') as f:
        content = f.read()
        content = re.split('\n', content)

        i = 0
        for line in content:
            if "element vertex" in line:
                number_vertices = [int(s) for s in line.split() if s.isdigit()]
            if "element face" in line:
                number_face = [int(s) for s in line.split() if s.isdigit()]
            if "end_header" in line:
                i += 1
                break
            else:
                i += 1

        for line in content[i: i+number_vertices[0]]:
            line = line.split()
            line = [float(x) for x in line]
            vertices_list.append(line)

        for line in content[(i+number_vertices[0]):-1]:
            line = [int(s) for s in line.split() if s.isdigit()]

            faces_list.append((line))
            if int(line[0]) != 3:
                # print(line)
                print(
                    "Presence of non triangulated faces detected. The script might not work properly.")
        
        # Remove number of points for edge
        faces_list = np.array(faces_list)
        faces_list = faces_list[:,1:]
    return number_vertices[0], np.array(vertices_list), number_face[0], faces_list


def open_ply(ply_file):
    try:
        vertNum, vertices, faceNum, faces = load_ply_data(ply_file)
        loaded_correctly = True
        return vertNum, vertices, faceNum, faces, loaded_correctly

    except Exception:
        loaded_correctly = False
        number_vertices = 0
        number_face = 0
        vertices_list = []
        faces_list = []
        print('Failed to read the file')
        return number_vertices, np.array(vertices_list), number_face, np.array(faces_list), loaded_correctly


# rpoly file contains center coordinate of helix, "generate" needs end coordiates of helix:
def move_along_vector(point, vector, length):
    move_distance = float(
        length) * 0.4 / 2.0  # 0.4 is the length of a base pair in oxDNA units, move half the helixlength down
    return [point[0] - move_distance * vector[0], point[1] - move_distance * vector[1],
            point[2] - move_distance * vector[2]]
