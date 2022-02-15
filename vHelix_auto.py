import sys
import re
from pyquaternion import Quaternion
import math

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from tacoxDNA.src.libs import cadnano_utils as cu


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

    # print (data, fwd_helix_connections, rev_helix_connections)
    return data, fwd_helix_connections, rev_helix_connections


def unit_vector(vector):  #from here https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):          #from here https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def angle(vector1, vector2):            #from here https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
    """ Returns the angle in radians between given vectors"""
    v1_u = unit_vector(vector1)
    v2_u = unit_vector(vector2)
    minor = np.linalg.det(
        np.stack((v1_u[-2:], v2_u[-2:]))
    )
    if minor == 0:
        sign = 1
    else:
        sign = -np.sign(minor)
    dot_p = np.dot(v1_u, v2_u)
    dot_p = min(max(dot_p, -1.0), 1.0)
    return sign * np.arccos(dot_p)

def move_along_vector(point, vector,
                      length):  # rpoly file contains center coordinate of helix, "generate" needs end coordiates of helix:
    move_distance = float(
        length) * 0.4 / 2.0  # 0.4 is the length of a base pair in oxDNA units, move half the helixlength down
    return [point[0] - move_distance * vector[0], point[1] - move_distance * vector[1],
            point[2] - move_distance * vector[2]]

def rotation_matrix_from_vectors(vec1, vec2):   #https://stackoverflow.com/questions/63525482/finding-the-rotation-matrix-between-two-vectors-in-python align to a parallel vector

    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix

def to_cadnano(data):
    tot_n_bp = 0
    # Reads orientation from the "data" and produces rotations from the Quaternian coordinates
    allowed_bases = {}
    allowed_for_vectors_base = {}

    generator = cu.StrandGenerator()

    start_bases = []
    vec_list = []
    vec2_list = []

    position_list = []
    new_position_list = []
    x_vector = np.array([1.0, 0.0, 0.0])

    for n, i in enumerate(data):

        position = [float(i[3]) / 0.84, float(i[4]) / 0.84,
                    float(i[5]) / 0.84]  # 0.84 scaling is ad hoc solution to get good looking models
        position_list.append(position)

        q = Quaternion(w=float(i[9]), x=float(i[6]), y=float(i[7]),
                       z=float(i[8]))  # find the helix rotation Info from file
        vec = q.rotate(np.array([0.0, 0.0, 1.0]))  # use it to figure out direction vec = q.rotate(np.array([0.0, 0.0, 1.0]))
        vec2 = q.rotate([0.65, -0.76, 0.0])   # ad hoc conversion between rpoly rotation and cadnano utils --> for oxDNA. vec2 = q.rotate([0.65, -0.76, 0.0]) vec2 = q.rotate([1.0, 0.0, 0.0])
                                             # vec2 is going to be used as the perpendicular vector to the direction,
                                             # for the first base
        n_bp = int(i[2])
        vec=unit_vector(vec)
        new_position = move_along_vector(position, vec, n_bp)       #calculate the position of start base
        new_position_list.append(new_position)
        vec_list.append(vec)
        vec2_list.append(vec2)
        #print(vec)

    #print(new_position_updated)
    v2_reference_list = []
    vector_from_origin = []
    for n, i in enumerate(data):

        vec = vec_list[n]
        vec2 = vec2_list[n]

        new_strand = generator.generate_or_sq(bp=n_bp, start_pos=new_position_list[n], direction=vec_list[n], perp=vec2_list[n]) #generate strands
        tot_n_bp += n_bp

        v2_reference = np.array([1.0,0.0,0.0]) #create a vector reference, to check the angle of the base direction, copied from cadnano_utilsnp.array([1.0, 0.0, 0.0]) np.random.randn(3)
        v2_reference -= v2_reference.dot(vec) * vec/np.linalg.norm(vec)**2 #the vector has to be perpendicular to vec, copied from cadnano_utils and https://stackoverflow.com/questions/33658620/generating-two-orthogonal-vectors-that-are-orthogonal-to-a-parti
        v2_reference /= np.linalg.norm(v2_reference)

        u = np.array([0.0, 0.0, 1.0]) #help from https://math.stackexchange.com/questions/2906314/how-to-calculate-angle-between-two-vectors-in-3d-with-clockwise-or-counter-clock

        #list of vectors
        sequence = new_strand[0]._get_Marco_output()  # get the sequence only for scaffold (like the one in .conf file)
        sequence_list = sequence.split('\n')
        base_coord_list = []
        #print(vec)

        for sequence in sequence_list:         #get coordinates
            base_coord = sequence.split(' ')
            base_coord_list.append(base_coord)

        vector = unit_vector(np.array([float(base_coord_list[0][0]), float(base_coord_list[0][1]), float(base_coord_list[0][2])]))
        #print(vector, n+1)
        vector_from_origin.append(vector)

        v2_reference = np.array(vector)  # create a vector reference, to check the angle of the base direction, copied from cadnano_utilsnp.array([1.0, 0.0, 0.0]) np.random.randn(3)
        v2_reference -= v2_reference.dot(vec) * vec / np.linalg.norm(vec) ** 2  # the vector has to be perpendicular to vec, copied from cadnano_utils and https://stackoverflow.com/questions/33658620/generating-two-orthogonal-vectors-that-are-orthogonal-to-a-parti
        v2_reference /= np.linalg.norm(v2_reference)

        v2_reference_list.append(v2_reference)

        allowed_for_helix = []
        allowed_for_helix_vectors = []
        vec_on_x_list = []
        R = (rotation_matrix_from_vectors(vec, x_vector))   #
        vec_on_x = np.dot((R), vec)
        vec_on_x_list.append(vec_on_x)

        for base in range(len(base_coord_list)-1):
            vector_from_coord = np.array(([float(base_coord_list[base][3]), float(base_coord_list[base][4]), float(base_coord_list[base][5])])) #obtain vector dir for each base
            vector_from_coord = unit_vector(vector_from_coord)
            vectors_angle_x = ((angle(vector_from_coord, v2_reference)))
            #print(vectors_angle_x, n + 1)
            #print((vectors_angle), n+1, len(base_coord_list))
            if vectors_angle_x > 0.52 and vectors_angle_x < 1.5708:
                allowed_for_helix.append(base)
                #print((np.cross(vec_on_x, vz_reference)), n+1)
                allowed_for_helix_vectors.append(
                    [float(base_coord_list[base][0]), float(base_coord_list[base][1]), float(base_coord_list[base][2]),
                     #float(vector_from_coord_rotated[0]), float(vector_from_coord_rotated[1]), float(vector_from_coord_rotated[2])])
                     float(base_coord_list[base][3]), float(base_coord_list[base][4]), float(base_coord_list[base][5])])
                #print(math.degrees(angle_rotated), n+1)
            if base == 0:
                start_bases.append([float(base_coord_list[base][0]), float(base_coord_list[base][1]), float(base_coord_list[base][2])])

        allowed_bases[n+1] = allowed_for_helix
        allowed_for_vectors_base[n+1] = allowed_for_helix_vectors

    fig = plt.figure()
    ax=fig.add_subplot(111, projection = '3d')

    for v1,coord in zip(vector_from_origin,start_bases):
        x,y,z = (coord[0], coord[1], coord[2])
        u,v,w = v1[0],v1[1],v1[2]
        #ax.quiver(x, y, z, u, v, w, length=10, color = 'blue')
    for vec,coord in zip(vec_list,start_bases):
        x,y,z = (coord[0], coord[1], coord[2])
        u,v,w = vec[0], vec[1], vec[2]
        ax.quiver(x, y, z, u, v, w, length=60, color = 'red')
    for vec2, coord in zip(v2_reference_list, start_bases):
        x, y, z = (coord[0], coord[1], coord[2])
        u, v, w = vec2[0], vec2[1], vec2[2]
        ax.quiver(x, y, z, u, v, w, length=10, color='green')

    for v1,coord in zip(vec_on_x_list,start_bases):
        x,y,z = (coord[0], coord[1], coord[2])
        u,v,w = v1[0],v1[1],v1[2]
        #ax.quiver(x, y, z, u, v, w, length=20, color = 'blue')


    for helix in allowed_for_vectors_base:
        i = 0
        for base in allowed_for_vectors_base[helix]:
            x,y,z = float(base[0]), float(base[1]), float(base[2])
            u,v,w = float(base[3]), float(base[4]), float(base[5])
            if i==len(allowed_for_vectors_base[helix])//2:
                ax.text(x,y,z,s=str(helix))
            i+=1
            #ax.quiver(x, y, z, u, v, w, length = 10, color='orange')
            ax.set_xlim(60, -60)
            ax.set_ylim(60, -60)
            ax.set_zlim(60, -60)

    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    plt.show()

    return allowed_bases


def create_cadnano(data, fwd_helix_connections, rev_helix_connections):
    print(data, fwd_helix_connections, rev_helix_connections)


if __name__ == "__main__":
    data, fwd_helix_connections, rev_helix_connections = open_rpoly("octahedron.rpoly")

    allowed_bases = to_cadnano(data)

    # print(len(position_list))
