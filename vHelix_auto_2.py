import sys
import re
from pyquaternion import Quaternion
import math
import json
import io
from pathlib import Path
import random
from virtual_scaffold import save_workbook
import numpy as np
from tacoxDNA.src.libs import cadnano_utils as cu
from collections import Counter


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
        if "end_header" in line:
            i+=1
            break
        else:
            i+=1

    faces_list = []
    for line in content[(i+number_vertices[0]):-1]:
        line = [int(s) for s in line.split() if s.isdigit()]
        faces_list.append((line))
        if int(line[0]) != 3:
            #print(line)
            print("Presence of non triangulated faces detected. The script might not work properly.")

    return faces_list

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

def move_along_vector(point, vector,length):  # rpoly file contains center coordinate of helix, "generate" needs end coordiates of helix:
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

def rotation_matrix(axis, theta):   #https://stackoverflow.com/questions/6802577/rotation-of-3d-vector
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

def vertices_from_ntrail(n_trail_list): #return a dictionary of the edges and the corresponding connections between the vertices

    I=1
    vert_to_helix_number = {}

    for current, next in zip(n_trail_list,n_trail_list[1:]):
        vert_to_helix_number[I]=(current, next)
        I+=1

    return vert_to_helix_number

def face_finder(edge,vert_to_helix_number,faces_list): #find to which faces the selected edge belongs


    for edge_number, vertices in vert_to_helix_number.items():  #find which vertices are connected by the selected edge
        if edge_number == edge:
            vertices_found=vertices

    faces_found = []
    for ele in faces_list:
        if vertices_found[0] in ele and vertices_found[1] in ele:
            if ele[0] == 3: #looking for triangular faces
                faces_found.append(ele[1:])
            else:
                print("The faces are not triangualarized. The algorithm might not work.")
                faces_found.append(ele[1:])


    return faces_found

def adjacent_edges(selected_edge, faces_found, vert_to_helix_number): #find edges adjacent to the selected edge, on the adjacent faces

    adjacent_edges_dict = {}

    for face in faces_found:
        mod_face = face.copy()
        mod_face.append(mod_face[0])
        for first, second in zip(mod_face, mod_face[1:]):   #these are the faces that have the same direction in the ply as the rpoly
            for edge, vertices in vert_to_helix_number.items():
                if vertices == (first, second):
                    adjacent_edges_dict[edge] = face
        for first, second in zip(mod_face, mod_face[1:]):   #these are the faces with the opposite direction
            for edge, vertices in vert_to_helix_number.items():
                if vertices == (second, first):
                    adjacent_edges_dict[edge] = face
    try:
        del adjacent_edges_dict[selected_edge]
    except KeyError:
        pass

    return adjacent_edges_dict  #the first face will have the same direction as the ply, the second one will have the opposite
                                # compare order of edges btw ply and rpoly. From ply the normal should be always outside: if the order not the same,
                                # I should take the negative vector from the rpoly
                                #for the double vertices, these should be the faces that should be connected (empirycally tried)

def double_vertices_faces(faces_list, vert_to_helix_number, to_reinforce):

    double_vertices_faces_list = []

    for face in faces_list:     #get a list of faces with the starting value (number of edges of the face) and a repetition of the first edge value
        mod_face = face.copy()
        mod_face.append(mod_face[1])
        for first, second in zip(mod_face[1:], mod_face[2:]):
            for edge, vertices in vert_to_helix_number.items():
                if vertices == (first, second):
                    double_vertices_faces_list.append(mod_face)

    faces_to_edges_list_old = []

    for a in (double_vertices_faces_list):
        edges_list = []
        for i in range(0, a[0]):
            try:
                edge_face = list(vert_to_helix_number.keys())[list(vert_to_helix_number.values()).index((a[i+1], a[i+2]))]
                edges_list.append(edge_face)
            except ValueError:
                continue
        faces_to_edges_list_old.append(edges_list)

    faces_to_edges_list = []

    for ele in faces_to_edges_list_old:
        if ele not in faces_to_edges_list:
            faces_to_edges_list.append(ele)

    return faces_to_edges_list

def double_edges_counter(vert_to_helix_number):

    double_edges = []
    double_edges_coupled = []


    for x in vert_to_helix_number.values():
        if x and x[::-1] in vert_to_helix_number.values():
            double_edges.append(list(vert_to_helix_number.keys())[list(vert_to_helix_number.values()).index(x)])
            double_edges_coupled.append((list(vert_to_helix_number.keys())[list(vert_to_helix_number.values()).index(x)], list(vert_to_helix_number.keys())[list(vert_to_helix_number.values()).index(x[::-1])]))

    for x in double_edges_coupled:
        if x[::-1] in double_edges_coupled:
            double_edges_coupled.remove(x[::-1])

    return double_edges, double_edges_coupled

def vectors_from_rpoly(data): #find vectors from rpoly

    position_list = []
    new_position_list = []
    vec_list = []
    vec2_list = []
    lengths_list = []
    tot_edges_vectors_dict = {}

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
        lengths_list.append(n_bp)
        #vec=unit_vector(vec)
        new_position = move_along_vector(position, vec, n_bp)       #calculate the position of start base
        new_position_list.append(new_position)
        vec_list.append(vec)
        vec2_list.append(vec2)
        tot_edges_vectors_dict[n+1] = vec   #dict with all edges and direction vectors

    return tot_edges_vectors_dict, position_list, lengths_list, vec_list, new_position_list, vec2_list

def adjacent_vectors(adjacent_edges_dict,selected_edge, vec_list):

    adj_edges_vectors_dict = {}

    adj_edges_vectors_dict[selected_edge] = vec_list[selected_edge - 1]

    for edge in adjacent_edges_dict:
        adj_edges_vectors_dict[edge] = vec_list[edge - 1]

    return adj_edges_vectors_dict

def normal_from_2vec(adj_edges_vectors_dict, selected_edge, adjacent_edges_dict): #find the normal btw 2 vectors, that will be the normal to the surface

    adjacent_edges_list = []

    for key in adjacent_edges_dict:
        adjacent_edges_list.append(key)

    normal_face_1 = np.cross(adj_edges_vectors_dict[adjacent_edges_list[1]], adj_edges_vectors_dict[selected_edge])

    normal_face_2 = np.cross( adj_edges_vectors_dict[selected_edge],adj_edges_vectors_dict[adjacent_edges_list[-1]]) #if i take the vectors in the opposite order, I should have the normal pointing "out",
                                                                                                             # because the two faces have different directions in the ply. Do they always have different directions?
                                                                                                             #Probably not

    return normal_face_1, normal_face_2

def angle_from_normals (normal_face_1, normal_face_2):

    angle = angle_between(normal_face_1, normal_face_2)

    return angle

def vector_sum (normal_face_1, normal_face_2):

    sum = normal_face_1 + normal_face_2

    return sum

def cross_vec(edge, tot_edges_vectors_dict, sum):

    cross = np.cross(tot_edges_vectors_dict[edge], sum)

    return  cross

def define_helices(cross, vec, position_list, edge, lengths_list): #https://math.stackexchange.com/questions/149434/how-do-you-find-a-vector-in-the-form-a-b-when-only-the-angle-and-magnitude-are/149435#149435
                        # https://stackoverflow.com/questions/11773889/how-to-calculate-a-vector-from-an-angle-with-another-vector-in-2d
                        #https://scipython.com/book/chapter-6-numpy/examples/creating-a-rotation-matrix-in-numpy/
    helix_0 = {}
    helix_1 = {}
    helix_2 = {}
    helix_3 = {}

    center_0 = np.array(position_list[edge-1])
    helix_0["center"] = center_0
    center_1 = center_0 + ((np.dot(cross,(rotation_matrix(vec, math.radians(135))))))
    helix_1["center"] = center_1
    center_2 = center_1 + (np.dot(cross, (rotation_matrix(vec, math.radians(45))))) #using center_0 + (np.dot(cross,
                                                                                    # (rotation_matrix(vec, math.radians(90))))) was giving problems (not the right position)
    helix_2["center"] = center_2
    center_3 = center_0 + (np.dot(cross, (rotation_matrix(vec, math.radians(45)))))
    helix_3["center"] = center_3

    helix_0["direction"] = vec
    helix_1["direction"] = - vec
    helix_2["direction"] = vec
    helix_3["direction"] = -vec

    helix_0["length"] = lengths_list[edge-1]
    helix_1["length"] = round(round(lengths_list[edge - 1]/(32/3)-2)*(32/3)) #shorter, so there is some space to avoid collisions. Leaves 21 bp in total(10 on one side and 11 on the other)
    helix_2["length"]= round(round(lengths_list[edge - 1]/(32/3)-2)*(32/3))
    helix_3["length"] = lengths_list[edge - 1]


    helices_dict = {}
    #helices_dict["Reinforced"] = "yes"
    helices_dict["helix_0"] = helix_0
    helices_dict["helix_1"] = helix_1
    helices_dict["helix_2"] = helix_2
    helices_dict["helix_3"] = helix_3

    return helices_dict

def scaffold_crossovers(helices_dict):


    possible_cross = (np.arange(0, 40+int(helices_dict['helix_1']['length']/2), (32/3))).tolist()

    possible_cross= np.round(possible_cross,0)


    central_cross = max(possible_cross)

    cross_zero_one = ([0, int(central_cross-1),1, int(central_cross-1)], [1, int(central_cross), 0, int(central_cross)])
    cross_one_two = ([1, 40, 2, 40], [2, 40 + int(helices_dict['helix_1']['length']) - 1, 1, 40 + int(helices_dict['helix_1']['length']) - 1])

    scaffold_cross = (cross_zero_one + cross_one_two)

    return scaffold_cross

def reinforced_bundle_creation(edge, vert_to_helix_number, faces_list, vec_list, tot_edges_vectors_dict, position_list, lengths_list):

    faces_found = face_finder(edge, vert_to_helix_number, faces_list)

    adjacent_edges_dict = adjacent_edges(edge, faces_found, vert_to_helix_number)

    adj_edges_vectors_dict = adjacent_vectors(adjacent_edges_dict, edge, vec_list)

    normal_face_1, normal_face_2 = normal_from_2vec(adj_edges_vectors_dict, edge, adjacent_edges_dict)

    angle_from_normals(normal_face_1, normal_face_2)

    sum = vector_sum(normal_face_1, normal_face_2)

    cross = cross_vec(edge, tot_edges_vectors_dict, sum)
    vec = tot_edges_vectors_dict[edge]

    helices_dict = define_helices(cross, vec, position_list, edge, lengths_list)

    scaffold_cross = scaffold_crossovers(helices_dict)

    return helices_dict, scaffold_cross

def bundles_colletion(bundle, bundles_dictionary, edge):

    bundles_dictionary["bundle_"+str(edge)] = bundle

    return bundles_dictionary

def helix_3d(vec_list, new_position_list, vec2_list, lengths_list): #build a 3d model of strands to calculate the autofill strand gaps distance

    generator = cu.StrandGenerator()
    new_strands_list = []
    helices_3d_dict = {}
    for n in range(len(vec_list)):

        new_strand = generator.generate_or_sq(bp=lengths_list[n], start_pos=new_position_list[n], direction=vec_list[n], perp=vec2_list[n])
        new_strands_list.append(new_strand)
        scaffold = new_strand[0]._get_Marco_output() #0 should be scaffold and 1 staples. Might be wrong
        scaffold_list = scaffold.split('\n')
        scaffold_list.remove('')
        staple = new_strand[1]._get_Marco_output()
        staple_list = staple.split('\n')
        staple_list.remove('')
        helices_3d_dict["helix"+str(n+1)] = scaffold_list, staple_list

    return helices_3d_dict

def auto_fill_strand_gaps(fwd_helix_connections, rev_helix_connections, helices_3d_dict):    #make it go over the entire structure. If the gaps are big, add bases

    length_base = 1.07 #"length" of a single base in ssdna 0.42, 0.55, 0.63, --> 0.694 in vhelix, 1.07 ad hoc to replicate vHelix here

    bases_to_add_scaffold = {}
    bases_to_add_staples = {}

    for n in range(len(fwd_helix_connections)): #fill gaps on scaffold. Is oxdna in 5'->3' or the opposite? Yes, opposite
        #base3_prime_scaffold = helices_3d_dict["helix"+str(n+1)][0][-1]     #take last base from the strand connected at the 3' and first of the one at 5'
        #base5_prime_scaffold = helices_3d_dict["helix"+str(n+2)][0][0]      #actually inverted, because oxdna is in the opposite order
        base3_prime_scaffold = helices_3d_dict["helix" + str(fwd_helix_connections[n][1])][0][0]    # inverted
        base5_prime_scaffold = helices_3d_dict["helix" + str(fwd_helix_connections[n][0])][0][-1]
        base3_prime_scaffold = base3_prime_scaffold.split(' ')
        base5_prime_scaffold = base5_prime_scaffold.split(' ')
        base5_prime_scaffold = np.array(
            (float(base5_prime_scaffold[0]), float(base5_prime_scaffold[1]), float(base5_prime_scaffold[2])))
        base3_prime_scaffold = np.array(
            (float(base3_prime_scaffold[0]), float(base3_prime_scaffold[1]), float(base3_prime_scaffold[2])))
        dist = (np.linalg.norm(base5_prime_scaffold - base3_prime_scaffold)-0.8)*0.8518  #from oxdna length to nm
        bases_to_add = round((dist) / length_base)
        bases_to_add_scaffold[str(fwd_helix_connections[n][0]) + ", " + str(fwd_helix_connections[n][1])] = bases_to_add

    for n in range(len(rev_helix_connections)): #fill gaps on staples

        #base3_prime_staple = helices_3d_dict["helix"+str(rev_helix_connections[n][0])][1][-1]    #not inverted
        #base5_prime_staple = helices_3d_dict["helix"+str(rev_helix_connections[n][1])][1][0]
        base3_prime_staple = helices_3d_dict["helix" + str(rev_helix_connections[n][1])][1][0] #inverted
        base5_prime_staple = helices_3d_dict["helix" + str(rev_helix_connections[n][0])][1][-1]
        base3_prime_staple = base3_prime_staple.split(' ')
        base5_prime_staple = base5_prime_staple.split(' ')
        base5_prime_staple = np.array(
            (float(base5_prime_staple[0]), float(base5_prime_staple[1]), float(base5_prime_staple[2])))
        base3_prime_staple = np.array(
            (float(base3_prime_staple[0]), float(base3_prime_staple[1]), float(base3_prime_staple[2])))
        dist = (np.linalg.norm(base5_prime_staple - base3_prime_staple)-0.8) * 0.8518  # from oxdna length to nm
        bases_to_add = round((dist) / length_base)
        bases_to_add_staples[str(rev_helix_connections[n][0]) + ", " + str(rev_helix_connections[n][1])] = bases_to_add

    return bases_to_add_scaffold, bases_to_add_staples

def takeFirst(elem):    #copied from below
    return elem[0]

def takeSecond(elem):   #https://www.programiz.com/python-programming/methods/list/sort
    return elem[1]

def add_staples(bundles_collection, to_reinforce):  #add staples. To all edges? Only to reinforced edges?

    base_crossovers = [[3, 39, 0, 39 ], [0,40, 3, 40],\
                         [1, 47, 0, 47], [0,48, 1, 48],\
                        [1, 55, 2, 55], [2, 56, 1, 56],\
                        [3, 63, 2, 63], [2, 64, 3, 64]]
    repeating_length = 32

    #guessed a bunch of things in the two loops, but it seems to work pretty okay.

    staples_crossovers = {}
    staples_nicks = {}


    for bundle, i in zip(bundles_collection.values(), bundles_collection.keys()):
        allowed_crossovers = base_crossovers.copy()
        nicks = []

        length_loop = (bundle[1][3][1]) #last crossover scaffold
        length_helix_0 = (bundle[0]['helix_0']['length']+30)  #using the length of the main vhelix instead of the length of the loop

        for cross in base_crossovers:
            moved_cross = cross[1]

            while moved_cross < length_helix_0 -(32) - 5: #to avoid taking unwanted crossovers
                moved_cross = moved_cross + repeating_length
                new_crossover = [cross[0], moved_cross, cross[2], moved_cross]
                allowed_crossovers.append(new_crossover)
        allowed_crossovers.sort(key=takeSecond)
        if len(allowed_crossovers)%2 !=0:
           allowed_crossovers =  allowed_crossovers[:-1]
        staples_crossovers[i]= allowed_crossovers   #save crossovers for each bundle

              #length_loop = (bundle[1][3][1]) #last crossover scaffold
               #check to remove crossovers after the loop crossover
        #if len(staples_crossovers)%2 !=0:
        #    staples_crossovers = staples_crossovers[:-1]
        #print(staples_crossovers)
        #if len(allowed_crossovers)%2 !=0:
        #    allowed_crossovers =  allowed_crossovers[:-1]

        for cross in list(allowed_crossovers):
            if cross[0] == 1 or cross[0] == 2:
                if cross[1] >= bundle[1][3][1] or cross[3] >= bundle[1][3][1]:
                    allowed_crossovers.remove(cross)
            if cross[2] == 1 or cross[2] == 2:
                 if cross[1] >= bundle[1][3][1] or cross[3] >= bundle[1][3][1]:
                     try:
                         allowed_crossovers.remove(cross)
                     except ValueError:
                         continue

        allowed_crossovers_1 = allowed_crossovers.copy()
        allowed_crossovers_1.remove([3, 39, 0, 39])     #remove to make things cleaner
        length = 0
        starting_point = 40
        starting_helix_1 = 1
        nick_point = starting_point


        should_restart = True
        while should_restart == True:       #oligo 1, starting at 1, 40
            should_restart = False

            if starting_point < (length_loop-23):
                if starting_helix_1 == 1 or starting_helix_1 == 3:  #odd helices
                    should_restart = True

                    for cross in allowed_crossovers_1:
                        if cross[0] == starting_helix_1 and cross[1] != starting_point + 1 and cross[1] != starting_point - 1:  #remove crossovers that are next to the one checked
                            length = length + cross[3] - nick_point + 1
                            starting_helix_1 = cross[2]
                            if length >= 32:
                                nick_point = nick_point + 16
                                length = 0
                                #print(nick_point, 'staple_1', cross[0])     #the save point will be the last base of the staple (arrow)
                                nicks.append([cross[0], nick_point-1])
                            starting_point = cross[3]
                            cross_index = (allowed_crossovers_1.index(cross))
                            try:
                                del allowed_crossovers_1[cross_index+1]
                            except IndexError:
                                continue
                            allowed_crossovers_1.remove(cross)
                            break

                if starting_helix_1 == 0 or starting_helix_1 == 2:
                    should_restart = True
                    for cross in allowed_crossovers_1:
                        if cross[0] == starting_helix_1 and cross[1] != starting_point + 1 and cross[1] != starting_point - 1:
                            length = length +  starting_point - cross[3] +1 #not nick_point,  because it's going "back"
                            starting_helix_1 = cross[2]
                            if length >= 32:
                                nick_point = nick_point + 16
                                length = 0
                                nicks.append([cross[0], nick_point - 1])
                            starting_point = cross[3]
                            cross_index = (allowed_crossovers_1.index(cross))
                            try:
                                del allowed_crossovers_1[cross_index - 1]
                            except IndexError:
                                continue
                            allowed_crossovers_1.remove(cross)
                            break
            elif starting_point >= (length_loop-23):
                length = length + (length_helix_0 - cross[3]) + 1

                should_restart = False

        allowed_crossovers_2 = allowed_crossovers.copy()
        allowed_crossovers_2.remove([0, 40, 3, 40])
        allowed_crossovers_2.remove([3, 39, 0, 39])
        starting_helix_2 = 2
        length = 0
        starting_point = 40
        nick_point = starting_point
        should_restart = True

        while should_restart == True:       #oligo 2, starting at 2, 40 (I think)
            should_restart = False
            if starting_point < (length_loop-12):
                if starting_helix_2 == 1 or starting_helix_2 == 3:  #odd helices
                    should_restart = True
                    for cross in allowed_crossovers_2:
                        if cross[2] == starting_helix_2 and cross[1] != starting_point + 1 and cross[1] != starting_point - 1:  #remove crossovers that are next to the one checked

                            length =  (length + starting_point -  cross[3]) + 1 #negative because it'g going on the opposite direction than the other one
                            starting_helix_2 = cross[0]
                            if length >= 32:
                                nick_point = nick_point + 16
                                length = 0
                                nicks.append([cross[2], nick_point])
                            starting_point = cross[3]
                            cross_index = (allowed_crossovers_2.index(cross))
                            try:
                                del allowed_crossovers_2[cross_index-1]
                            except IndexError:
                                continue
                            allowed_crossovers_2.remove(cross)
                            break


                if starting_helix_2 == 0 or starting_helix_2 == 2 :
                    should_restart = True
                    for cross in allowed_crossovers_2:
                        if cross[2] == starting_helix_2 and cross[1] != starting_point + 1 and cross[1] != starting_point - 1:
                            length = (length +  cross[3] -  nick_point ) +1 #not nick_point,  because it's going "back"
                            starting_helix_2 = cross[0]
                            if length >= 32:
                                nick_point = nick_point + 16
                                nicks.append([cross[2], nick_point])
                                length = 0
                            starting_point = cross[3]
                            cross_index = (allowed_crossovers_2.index(cross))
                            try:
                                del allowed_crossovers_2[cross_index + 1]
                            except IndexError:
                                continue
                            allowed_crossovers_2.remove(cross)
                            break

            elif starting_point >= (length_loop-12):
                length = length + (length_helix_0 - cross[3]) + 1
                should_restart = False
        staples_nicks[i] = nicks

    return staples_crossovers, staples_nicks

def scaff_vhelix_crossovers(lengths_list, gaps_dimensions, bundles, double_edges_coupled):

    gaps_filled_scaffold_lengths = []
    scaff_cross_vhelix = []

    for n, element in zip((lengths_list), gaps_dimensions[0].keys()):   #add the gaps filling to the edges
        gaps_filled_scaffold_lengths.append(n+int(gaps_dimensions[0][element]))

    for n, length in zip(gaps_dimensions[0].keys(), gaps_filled_scaffold_lengths):   #update bundles with the new lengths
        try:
            n = n.split(', ')
            bundles['bundle_'+str(n[0])][0]['helix_0'].update({'length' : length})
            scaff_cross_vhelix.append([int(n[0]), 30 +length, int(n[1]), 30])
        except KeyError:
            scaff_cross_vhelix.append([int(n[0]), 30+length, int(n[1]), 30])
            continue

    index = 0

    for x in scaff_cross_vhelix:
        for edge in double_edges_coupled:
            if x[2] == edge[1]:
                length = 30 + bundles['bundle_'+str(edge[1])][0]['helix_0']['length']
                scaff_cross_vhelix[index][3] = length - 2

            if x[0] == edge[1]:
                scaff_cross_vhelix[index][1] = 30 + 1
        index+=1


    return bundles, gaps_filled_scaffold_lengths, scaff_cross_vhelix   #return the updated bundles and all the updated lengths and crossovers, to use for the normal edges

def staple_vhelix_crossovers(gaps_dimensions, bundles, lengths_list, double_edges_coupled): #using old length list so that scaffold stands ssdna

    staple_cross_vhelix = []

    for n, length, gap in zip(gaps_dimensions[1].keys(), lengths_list, gaps_dimensions[1].values()):
        n = n.split(', ')
        staple_cross_vhelix.append([int(n[0]), 30 - gap, int(n[1]), length +30])

    index = 0

    for x in (staple_cross_vhelix):
        for edge in double_edges_coupled:
            if x[0] ==edge[1]:
                gap = 30-staple_cross_vhelix[index][1]
                staple_cross_vhelix[index][1] = lengths_list[edge[1]-1]+30
            if x[2] == edge[1]:
                staple_cross_vhelix[index][3] = 30 - gap
        index += 1

    #print(staple_cross_vhelix)

    return staple_cross_vhelix

def output_cadnano_2(filename, data, to_reinforce, gaps_filled_scaffold_lengths, bundles, scaff_cross_vhelix, double_vertices, double_edges,double_edges_coupled, staples_crossovers, staples_nicks, staple_cross_vhelix, lengths_list, faces_list, vert_to_helix_number):

    vh_list = []
    vh_order = {}

    max_length = (int(max(gaps_filled_scaffold_lengths)/32)*32)+64

    start_row = 1
    start_col = 1


    double_edges_firsts = []
    double_edges_seconds = []

    for x in double_edges_coupled:
        double_edges_firsts.append(takeFirst(x))
        double_edges_seconds.append(takeSecond(x))

    i = 0


    for n in (range(len(data))):
        if n+1 in to_reinforce:
            if n+1 in double_edges_firsts:
                if i %2 != 0:
                    i+=1
                    vh_order[n + 1] = (i-1,)
                else:
                    vh_order[n+1] = i, i+1, i+2
                    i = i+3
            elif n+1 in double_edges_seconds:
                if i %2 == 0:
                    vh_order[n + 1] = (i+1,)
                    i = i + 1
                else:
                    vh_order[n+1] = (i, )
                    i = i+1
            else:
                if i %2 != 0:
                    #i+=1
                    vh_order[n + 1] = i+1, i+2, i + 3, i + 4
                    i= i + 5
                else:
                    vh_order[n+1] = i, i+1, i+2, i+3
                    i = i+4
        else:
            # vh_order[n+1] = i,
            # i = i+2
            vh_order[n + 1] = i, i + 1, i + 2, i + 3
            i = i + 4
    print((vh_order))
    #stap_colors = [[55,13369344],[87,29184],[119,12060012],[151,16225054],[190,7536862]]    #[55,13369344],[87,29184],[119,12060012],[151,16225054],[190,7536862]

    row = start_row
    col = start_col
    insts = [0 for i in range(max_length)]
    skips = [0 for i in range(max_length)]
    deletions_list = []


    for edge in range(1, len(data)+1):

        id_num = vh_order[edge][0]

        if edge%12 == 0:
            row = row+4
            col = start_col

        if edge in to_reinforce and edge not in double_edges:

            #if edge + 1 in double_edges:
            #    print(edge+1)

            vh_dict = {}
            vh_dict = {"row": row,           #for the main helix
                 "col": col,
                 "num": id_num,
                 "scaf": get_scaf_vhelix_reinforced(max_length, bundles, scaff_cross_vhelix, edge, id_num, vh_order,
                                                    double_edges),
                 "stap": get_stap_vhelix_reinforced(max_length,staples_crossovers, staples_nicks, staple_cross_vhelix, edge, id_num, vh_order),
                 "loop": insts,
                 "skip": insert_deletions(bundles, edge, skips),
                 "scafLoop": [],
                 "stapLoop": [],
                 "stap_colors": get_colors(staples_nicks, edge, helix_num=0)}
            vh_list.append(vh_dict)

            vh_dict = {}

            vh_dict = {"row": row,   #for the first loop
                 "col": col-1,
                 "num": vh_order[edge][1],
                 "scaf": get_scaf_loop_1(max_length, bundles, edge, id_num),
                 "stap": get_stap_loop1(max_length, staples_crossovers, staples_nicks, edge, id_num, bundles),
                 "loop": insts,
                 "skip": insert_deletions(bundles, edge, skips),
                 "scafLoop": [],
                 "stapLoop": [],
                 "stap_colors": get_colors(staples_nicks, edge, helix_num=1)}
            vh_list.append(vh_dict)

            vh_dict = {}

            vh_dict = {"row": row-1,   #for the second loop
                 "col": col-1,
                 "num": vh_order[edge][2],
                 "scaf": get_scaf_loop_2(max_length, bundles, edge, id_num),
                 "stap": get_stap_loop2(max_length,staples_crossovers, staples_nicks, edge, id_num, bundles),
                 "loop": insts,
                 "skip": insert_deletions(bundles, edge, skips),
                 "scafLoop": [],
                 "stapLoop": [],
                 "stap_colors": get_colors(staples_nicks, edge, helix_num=2)}
            vh_list.append(vh_dict)

            vh_dict = {}

            vh_dict = {"row": row-1,   #for the virtual scaffold
                 "col": col,
                 "num": vh_order[edge][3],
                 "scaf": get_virtual_scaff(max_length, edge, id_num, double_vertices, vh_order, scaff_cross_vhelix,
                                           to_reinforce, bundles, faces_list, vert_to_helix_number),
                 "stap": get_stap_virtual(max_length,staples_crossovers, staples_nicks, edge, id_num,  bundles),
                 "loop": insts,
                 "skip": insert_deletions(bundles, edge, skips),
                 "scafLoop": [],
                 "stapLoop": [],
                 "stap_colors": get_colors(staples_nicks, edge, helix_num=3)}

            vh_list.append(vh_dict)
            row = row
            col = col+4
            deletions_list.append(deletions_number(bundles, edge))

        elif edge in to_reinforce and edge in double_edges_firsts:

            vh_dict = {}
            vh_dict = {"row": row,  # for non-reinforced edges
                       "col": col,
                       "num": vh_order[edge][0],
                       "scaf": get_scaff_reinforced_vhelix_BH(max_length, vh_order, id_num, edge, scaff_cross_vhelix, bundles),
                       "stap": get_stap_vhelix_reinforced_BH(max_length, staples_crossovers, staples_nicks, staple_cross_vhelix, edge, id_num, vh_order,double_edges_coupled),
                       "loop": insts,
                       "skip": skips,
                       "scafLoop": [],
                       "stapLoop": [],
                       "stap_colors": [[55, 13369344]]}
            vh_list.append(vh_dict)

            vh_dict = {}

            vh_dict = {"row": row,  # for the first loop
                       "col": col - 1,
                       "num": vh_order[edge][1],
                       "scaf": get_scaf_loop_1(max_length, bundles, edge, id_num),
                       "stap": get_stap_loop1_BH(max_length, staples_crossovers, staples_nicks, edge, id_num, bundles),
                       "loop": insts,
                       "skip": insert_deletions(bundles, edge, skips),
                       "scafLoop": [],
                       "stapLoop": [],
                       "stap_colors": get_colors(staples_nicks, edge, helix_num=1)}
            vh_list.append(vh_dict)

            vh_dict = {}

            vh_dict = {"row": row - 1,  # for the second loop
                       "col": col - 1,
                       "num": vh_order[edge][2],
                       "scaf": get_scaf_loop_2(max_length, bundles, edge, id_num),
                       "stap": get_stap_loop2_BH(max_length, staples_crossovers, staples_nicks, edge, id_num, bundles,
                                                 double_edges_coupled, vh_order),
                       "loop": insts,
                       "skip": insert_deletions(bundles, edge, skips),
                       "scafLoop": [],
                       "stapLoop": [],
                       "stap_colors": get_colors(staples_nicks, edge, helix_num=2)}
            vh_list.append(vh_dict)

            vh_dict = {}
            row = row
            col = col + 4
            deletions_list.append(deletions_number(bundles, edge))

        elif edge in to_reinforce and edge in double_edges_seconds:     #the second of the couple acts as a virtual scaffold in the double edges


            vh_dict = {}
            vh_dict = {"row": row-1,  # for non-reinforced edges
                       "col": col,
                       "num": vh_order[edge][0],
                       "scaf": get_scaff_vhelix_double_edge_BH(max_length, vh_order, id_num, edge, scaff_cross_vhelix),
                       "stap": get_stap_double_edge_BH(max_length, staple_cross_vhelix, edge, id_num, vh_order,
                                                       lengths_list, double_edges_coupled, staples_crossovers, staples_nicks),
                       "loop": insts,
                       "skip": skips,
                       "scafLoop": [],
                       "stapLoop": [],
                       "stap_colors": [[55, 13369344]]}
            vh_list.append(vh_dict)


            row = row
            col = col + 4
            deletions_list.append(deletions_number(bundles, edge))


        else:

            vh_dict = {}

            vh_dict = {"row": row,   #for non-reinforced edges
                 "col": col,
                 "num": vh_order[edge][0],
                 "scaf": get_scaff_vhelix(max_length, vh_order, vh_order[edge][0], edge, scaff_cross_vhelix),
                 "stap": get_stap_vhelix(max_length, staple_cross_vhelix, edge, vh_order[edge][0], vh_order),
                 "loop": insts,
                 "skip": skips,
                 "scafLoop": [],
                 "stapLoop": [],
                 "stap_colors": [[55,13369344]]}
            vh_list.append(vh_dict)
            row=row
            col = col +4
            deletions_list.append(deletions_number(bundles, edge))


    obj = {"name": filename, "vstrands": vh_list}

    json_string = json.dumps(obj, separators=(',', ':'))

    #os.makedirs(os.path.dirname('/' + filename[:-5]), exist_ok=True)
    Path("./" + str(filename[:-5])).mkdir(parents = True, exist_ok=True)
    #filename = os.path.join(folder, filename)
    file_save = "./" + str(filename[:-5]) + "/" + filename

    with io.open(file_save, 'w') as f:
        f.write(json_string)

    (scaffold_length(bundles))
    return deletions_list #, double_edges

def get_scaf_vhelix_reinforced(max_length, bundles,  scaff_cross_vhelix, edge, id_num, vh_order, double_edges):   #create scaffold for helix_0 reinforced

    helix_number = 0
    mesh = [[-1,-1,-1,-1] for x in range(max_length)]


    vhelix_cross_finder = list(item[2] for item in scaff_cross_vhelix)  #list to find the start of the crossover
    vhelix_end_finder = list(item[0] for item in scaff_cross_vhelix)    #list to find the end of the crossover


    for n in vhelix_cross_finder:   #find connections starting from particular edge
        if n == edge:
            index = (vhelix_cross_finder.index(n))
            helix_start = (scaff_cross_vhelix[index])

    for n in vhelix_end_finder: #find connections ending in particular edge
        if n == edge:
            index_1 = (vhelix_end_finder.index(n))
            helix_end = (scaff_cross_vhelix[index_1])

    helix_start_2 = [vh_order[helix_start[0]][0],helix_start[1] - 1 ,id_num,helix_start[3] + 1]

    mesh[30] = helix_start_2

    for x in range(31,helix_end[1]):
        mesh[x]= [id_num,x-1,id_num,x+1]
        #helix_start[3] = helix_start[3] + 1


    helix_end_2 = [id_num,helix_end[1]-1,vh_order[helix_end[2]][0],helix_end[3]]
    mesh[helix_end[1]] = helix_end_2

    loop_cross_1 = (bundles['bundle_' + str(edge)][1][0])
    loop_cross_2 = (bundles['bundle_' + str(edge)][1][1])

    loop_cross_1_2 = [vh_order[edge][0], loop_cross_1[1]-1, vh_order[edge][1], loop_cross_1[3]]
    loop_cross_2_2 = [vh_order[edge][1], loop_cross_2[1], vh_order[edge][0], loop_cross_2[3]+1]
    mesh[loop_cross_1[1]] = loop_cross_1_2
    mesh[loop_cross_2[1]] = loop_cross_2_2

    #print(mesh[29:])

    return mesh

def get_scaf_loop_1(max_length, bundles, edge, id_num):

    mesh = [[-1,-1,-1,-1] for x in range(max_length)]

    cross_start = bundles['bundle_'+str(edge)][1][2]
    cross_end = bundles['bundle_'+str(edge)][1][3]
    length_loop  = bundles['bundle_' + str(edge)][0]['helix_1']['length']

    id_num_loop_1 = id_num+1
    id_num_loop_2 = id_num+2

    mid_cross_1 = bundles['bundle_'+str(edge)][1][0]    #central cross, from helix_0 to helix_1
    mid_cross_2 = bundles['bundle_'+str(edge)][1][0]

    cross_start_rigth = [id_num_loop_1, cross_start[1]+1, id_num_loop_2, cross_start[3]]
    cross_end_rigth = [id_num_loop_2, cross_end[3], id_num_loop_1, cross_end[1]-1]

    mesh[40] = cross_start_rigth
    for x in range(cross_start[1]+1, cross_end[1]):
        mesh[x]= [id_num_loop_1,x+1,id_num_loop_1,x-1]  #order reverse because it'in the reverse direction

    mesh[40+length_loop-1] = cross_end_rigth

    cross_1 = bundles['bundle_'+str(edge)][1][0][1] #position from 0 to 1
    cross_2 = bundles['bundle_'+str(edge)][1][1][1] #position from 1 to 0

    mid_cross_1_2 = [id_num, cross_1, id_num_loop_1, cross_1 -1 ]
    mid_cross_2_2 = [id_num_loop_1, cross_2+1, id_num, cross_2 ]

    mesh[cross_1] = mid_cross_1_2
    mesh[cross_2] = mid_cross_2_2

    return mesh

def get_scaf_loop_2(max_length, bundles, edge, id_num):

    mesh = [[-1,-1,-1,-1] for x in range(max_length)]

    cross_start = bundles['bundle_'+str(edge)][1][2]
    cross_end = bundles['bundle_'+str(edge)][1][3]
    length_loop  = bundles['bundle_' + str(edge)][0]['helix_2']['length']
    id_num_loop_1 = id_num+1
    id_num_loop_2 = id_num+2

    cross_start_rigth = [id_num_loop_1, cross_start[1], id_num_loop_2, cross_start[3]+1]
    cross_end_rigth = [id_num_loop_2, cross_end[3]-1, id_num_loop_1, cross_end[1]]
    mesh[40] = cross_start_rigth

    for x in range(cross_start[1]+1, cross_end[1]):
        mesh[x]= [id_num_loop_2,x-1,id_num_loop_2,x+1]  #order reverse because it'in the reverse direction

    mesh[40+length_loop-1] = cross_end_rigth

    return mesh

def get_virtual_scaff(max_length, edge, id_num, double_vertices, vh_order, scaff_cross_vhelix,to_reinforce,bundles, faces_list, vert_to_helix_number):

    mesh = [[-1,-1,-1,-1] for x in range(max_length)]

    id_num_virt = id_num + 3
    length_virt = bundles['bundle_' + str(edge)][0]['helix_3']['length']  #shorter than helix_0, because no autofill strand gaps
    mesh[30]= [id_num_virt, 31, -1, -1] #not connected

    for x in range(31, 30+length_virt-1):   #-1 otherwise they are too long by 1
        mesh[x] = [id_num_virt, x + 1, id_num_virt, x - 1]
        # helix_start[3] = helix_start[3] + 1

    mesh[30+length_virt-1] = [-1, -1,id_num_virt, x]

    if double_vertices == True:
        #virt_scaff_cross_previous_edge = 1
        #virt_scaff_cross_next_edge = 1
        faces_to_edges_list = double_vertices_faces(faces_list, vert_to_helix_number, to_reinforce)

        for face in faces_to_edges_list:

            if edge in face:
                virt_scaff_cross_previous_edge = face[(face.index(edge))-1]
                try:
                    virt_scaff_cross_next_edge = face[(face.index(edge))+1]
                except IndexError:
                    virt_scaff_cross_next_edge = face[0]

        length_virt = bundles['bundle_' + str(edge)][0]['helix_0'][
            'length']  # as long as helix_0, because I want autofill strand gaps
        try:

            #length_next_edge = bundles['bundle_' + str(virt_scaff_cross_next_edge)][0]['helix_0']['length']
            length_previous_edge = bundles['bundle_' + str(virt_scaff_cross_previous_edge)][0]['helix_0']['length']
            mesh[30] = [id_num_virt, 30 + 1, vh_order[virt_scaff_cross_previous_edge][3], 30 + length_previous_edge - 1]
            mesh[30 + length_virt - 1] = [vh_order[virt_scaff_cross_next_edge][3], 30, id_num_virt,
                                          30 + length_virt - 2]

        except UnboundLocalError:

            length_virt = bundles['bundle_' + str(edge)][0]['helix_3'][
                'length']  # as long as helix_0, because I do not want autofill strand gaps
            mesh[30] = [id_num_virt, 31, -1, -1]
            mesh[30 + length_virt - 1] = [-1, -1, id_num_virt, x]

        except IndexError:

            length_virt = bundles['bundle_' + str(edge)][0]['helix_3'][
                'length']  # as long as helix_0, because I do not want autofill strand gaps
            mesh[30] = [id_num_virt, 31, -1, -1]
            mesh[30 + length_virt - 1] = [-1, -1, id_num_virt, x]



        for x in range(31, 30 + length_virt - 1):  # -1 otherwise they are too long by 1
            mesh[x] = [id_num_virt, x + 1, id_num_virt, x - 1]

        for face in faces_to_edges_list:    #avoid connections to stuff that does not have a virtual scaffold
            if edge in face:
                if len(face)>2:
                    if face[0] not in to_reinforce or face[2] not in to_reinforce or face[2] not in to_reinforce:
                        mesh = [[-1, -1, -1, -1] for x in range(max_length)]
                        length_virt = bundles['bundle_' + str(edge)][0]['helix_3'][
                            'length']  # as long as helix_3, because I do not want autofill strand gaps
                        for x in range(31, 30 + length_virt - 1):  # -1 otherwise they are too long by 1
                            mesh[x] = [id_num_virt, x + 1, id_num_virt, x - 1]
                        mesh[30] = [id_num_virt, 31, -1, -1]
                        mesh[30 + length_virt - 1] = [-1, -1, id_num_virt, x]

                if len(face) < 2:
                    mesh = [[-1, -1, -1, -1] for x in range(max_length)]
                    length_virt = bundles['bundle_' + str(edge)][0]['helix_3'][
                        'length']  # as long as helix_3, because I do not want autofill strand gaps
                    for x in range(31, 30 + length_virt - 1):  # -1 otherwise they are too long by 1
                        mesh[x] = [id_num_virt, x + 1, id_num_virt, x - 1]
                    mesh[30] = [id_num_virt, 31, -1, -1]
                    mesh[30 + length_virt - 1] = [-1, -1, id_num_virt, x]


        #for face in faces_to_edges_list:
            #if len(face)>2:
                #print(face)
    #print(mesh[29:])
    return mesh

def get_stap_vhelix_reinforced(max_length, staples_crossovers, staples_nicks, staple_cross_vhelix, edge, id_num, vh_order): #staple for helix_0

    mesh = [[-1,-1,-1,-1] for x in range(max_length)]

    helix_num = 0

    for x in staple_cross_vhelix:
        if edge == x[0] or edge == x[2]:
            if x[2] == edge:
                max_length_staples = x[3] #+30
                end_vhel_cross = x
            else:
                start_vhel_cross = x

    helix_0_crosses = []


    for x in range(start_vhel_cross[1], max_length_staples):    #starts from when the first staple starts
        mesh[x] = [id_num, x+1, id_num, x-1]

    for n in staples_crossovers['bundle_'+str(edge)]:
        if helix_num in n:
            helix_0_crosses.append(n)

    for x in helix_0_crosses:

        if helix_0_crosses.index(x)%2 != 0:
            mesh[x[1]-1] = [x[2]+id_num, x[1]-1, x[0] +id_num , x[3]-2]
        else:
            mesh[x[1]+1] = [x[2]+id_num, x[1]+2, x[0]+id_num, x[1]+1]

    mesh[start_vhel_cross[1]] = [vh_order[start_vhel_cross[0]][0], start_vhel_cross[1] +1,    #vhelix cross at the 3'
                                 vh_order[start_vhel_cross[2]][0], start_vhel_cross[3]]

    mesh[end_vhel_cross[3]] = [vh_order[end_vhel_cross[0]][0], end_vhel_cross[1],       #vhelix cross at 5'
                                 vh_order[end_vhel_cross[2]][0], end_vhel_cross[3] -1]



    helix_0_nicks = []
    for x in staples_nicks['bundle_'+str(edge)]:
        if 0 in x:                                      #we are in helix_0
            helix_0_nicks.append(x)

    for x in helix_0_nicks:
        mesh[x[1]-1] = [-1, -1, id_num, x[1]-2]
        mesh[x[1]] = [id_num, x[1] +1, -1, -1]

    #print(mesh[27:])
    return mesh

def get_stap_loop1(max_length,staples_crossovers, staples_nicks, edge, id_num, bundles):

    mesh = [[-1, -1, -1, -1] for x in range(max_length)]
    id_num_loop_1 = id_num + 1
    helix_num = 1
    max_staple_length = (bundles['bundle_'+str(edge)][0]['helix_'+str(helix_num)]['length'])+39

    for x in range (41, max_staple_length):
        mesh[x] = [id_num_loop_1, x-1, id_num_loop_1, x+1]

    helix_1_crosses = []

    for n in staples_crossovers['bundle_' + str(edge)]:
        if helix_num in n:
            helix_1_crosses.append(n)

    for x in helix_1_crosses:
        if helix_1_crosses.index(x) % 2 != 0:
            mesh[x[1] - 1] = [id_num_loop_1, x[1] - 2, x[0]+id_num, x[3] - 1]
        else:
            mesh[x[1] + 1] = [x[2]+id_num, x[1] + 1, id_num_loop_1, x[1] + 2]

    helix_1_nicks = []
    for x in staples_nicks['bundle_' + str(edge)]:
        if 1 in x:  # we are in helix_1
            helix_1_nicks.append(x)

    for x in helix_1_nicks:
        mesh[x[1]+1] = [-1, -1, id_num_loop_1, x[1]+2]
        mesh[x[1]] = [id_num_loop_1, x[1]-1, -1, -1]

    mesh[40] = [-1, -1, id_num_loop_1, 41]
    mesh[max_staple_length] = [id_num_loop_1, max_staple_length-1, -1, -1]

    return mesh

def get_stap_loop2(max_length,staples_crossovers, staples_nicks, edge, id_num, bundles):

    mesh = [[-1, -1, -1, -1] for x in range(max_length)]
    id_num_loop_2 = id_num + 2
    helix_num = 2
    max_staple_length = (bundles['bundle_'+str(edge)][0]['helix_'+str(helix_num)]['length'])+39

    for x in range (41, max_staple_length):
        mesh[x] = [id_num_loop_2, x+1, id_num_loop_2, x-1]

    helix_2_crosses = []

    for n in staples_crossovers['bundle_' + str(edge)]:
        if helix_num in n:
            helix_2_crosses.append(n)

    for x in helix_2_crosses:
        if helix_2_crosses.index(x) % 2 != 0:
            mesh[x[1] - 1] = [x[2]+id_num, x[1] - 1, id_num_loop_2, x[3] - 2]
        else:
            mesh[x[1] + 1] = [id_num_loop_2, x[1] + 2, x[0]+id_num, x[1] + 1]

    helix_2_nicks = []
    for x in staples_nicks['bundle_' + str(edge)]:
        if 2 in x:  # we are in helix_2
            helix_2_nicks.append(x)

    for x in helix_2_nicks:
        mesh[x[1] - 1] = [-1, -1, id_num_loop_2, x[1] - 1]
        mesh[x[1]] = [id_num_loop_2, x[1] +2, -1, -1]

    mesh[40] = [id_num_loop_2, 41, -1, -1]
    mesh[max_staple_length] = [-1, -1, id_num_loop_2, max_staple_length-1]

    return mesh

def get_stap_virtual(max_length,staples_crossovers, staples_nicks, edge, id_num,  bundles):

    mesh = [[-1, -1, -1, -1] for x in range(max_length)]

    id_num_virt = id_num + 3
    helix_num = 3
    max_staple_length = (bundles['bundle_' + str(edge)][0]['helix_' + str(helix_num)]['length']) + 30

    for x in range(31, max_staple_length-1):
        mesh[x] = [id_num_virt, x - 1, id_num_virt, x + 1]

    helix_3_crosses = []

    for n in staples_crossovers['bundle_' + str(edge)]:
        if helix_num in n:
            helix_3_crosses.append(n)

    for x in helix_3_crosses:
        if helix_3_crosses.index(x) % 2 != 0:
            mesh[x[1] - 1] = [id_num_virt, x[1] - 2, x[0] + id_num, x[3] - 1]
        else:
            mesh[x[1] + 1] = [x[2]+id_num, x[1] + 1, id_num_virt, x[1] + 2]

    helix_3_nicks = []

    for x in staples_nicks['bundle_' + str(edge)]:
        if 3 in x:  # we are in helix_3
            helix_3_nicks.append(x)

    for x in helix_3_nicks:
        mesh[x[1]+1] = [-1, -1, id_num_virt, x[1] + 2]
        mesh[x[1]] = [id_num_virt, x[1] - 1, -1, -1]

    mesh[30]  = [-1, -1, id_num_virt, 31]
    mesh[max_staple_length-1] = [id_num_virt, max_staple_length - 1, -1, -1]

    return mesh

def get_stap_vhelix(max_length, staple_cross_vhelix, edge, id_num, vh_order):

    mesh = [[-1, -1, -1, -1] for x in range(max_length)]

    for x in staple_cross_vhelix:
        if edge == x[0] or edge == x[2]:
            if x[2] == edge:
                max_length_staples = x[3]#+30
                end_vhel_cross = x
            else:
                start_vhel_cross = x


    for x in range(start_vhel_cross[1], max_length_staples):    #starts from when the first staple starts
        mesh[x] = [id_num, x+1, id_num, x-1]

    mesh[start_vhel_cross[1]] = [vh_order[start_vhel_cross[0]][0], start_vhel_cross[1]+1,  # vhelix cross at the 3'
                                 vh_order[start_vhel_cross[2]][0], start_vhel_cross[3]  -1]

    mesh[end_vhel_cross[3] -1] = [vh_order[end_vhel_cross[0]][0], end_vhel_cross[1],  # vhelix cross at 5'
                                    vh_order[end_vhel_cross[2]][0], end_vhel_cross[3] -2]

    mesh[round((max_length_staples+30)/2)] = [id_num, [round((max_length_staples+30)/2)], -1, -1]
    mesh[round((max_length_staples+30)/2)+1] = [-1, -1, id_num, round((max_length_staples+30)/2)+1]

    return mesh

def get_scaff_vhelix(max_length, vh_order, id_num, edge,scaff_cross_vhelix):   #non reinforced scaffold

    mesh = [[-1,-1,-1,-1] for x in range(max_length)]

    vhelix_cross_finder = list(item[2] for item in scaff_cross_vhelix)  # list to find the start of the crossover
    vhelix_end_finder = list(item[0] for item in scaff_cross_vhelix)  # list to find the end of the crossover

    for n in vhelix_cross_finder:  # find connections starting from particular edge
        if n == edge:
            index = (vhelix_cross_finder.index(n))
            helix_start = (scaff_cross_vhelix[index])

    for n in vhelix_end_finder:  # find connections ending in particular edge
        if n == edge:
            index_1 = (vhelix_end_finder.index(n))
            helix_end = (scaff_cross_vhelix[index_1])

    helix_start_2 = [vh_order[helix_start[0]][0], helix_start[1], id_num,
                     helix_start[3] + 1]  # +1 because cadnano seems to be abit weird, not sure it will work

    mesh[helix_start_2[3]-1] = helix_start_2

    for x in range(31, helix_end[1]-1):
        mesh[x] = [id_num, x - 1, id_num, x + 1]
        # helix_start[3] = helix_start[3] + 1

    helix_end_2 = [id_num, helix_end[1] - 2, vh_order[helix_end[2]][0], helix_end[3]]
    mesh[helix_end[1]-1] = helix_end_2

    #print(scaff_cross_vhelix)
    #print(mesh[29:])
    return mesh

# def get_stap_double_edge(max_length, staple_cross_vhelix, edge, id_num, vh_order):
#
#     mesh = [[-1, -1, -1, -1] for x in range(max_length)]
#
#     helix_num = 0
#
#     for x in staple_cross_vhelix:
#         if edge == x[0] or edge == x[2]:
#             if x[2] == edge:
#                 max_length_staples = x[3] #+ 30
#                 end_vhel_cross = x
#             else:
#                 start_vhel_cross = x
#
#     helix_0_crosses = []
#
#     for x in range(start_vhel_cross[1], max_length_staples ):  # starts from when the first staple starts
#         mesh[x] = [id_num, x + 1, id_num, x - 1]
#
#     for n in staples_crossovers['bundle_' + str(edge)]:
#         if helix_num in n and 3 in n:
#             helix_0_crosses.append(n)
#
#     # for x in helix_0_crosses:
#     #
#     #     if helix_0_crosses.index(x) % 2 != 0:
#     #         mesh[x[1] - 1] = [x[2] + id_num, x[1] - 1, x[0] + id_num, x[3] - 2]
#     #     else:
#     #         mesh[x[1] + 1] = [x[2] + id_num, x[1] + 2, x[0] + id_num, x[1] + 1]
#
#     mesh[start_vhel_cross[1]] = [vh_order[start_vhel_cross[0]][0], start_vhel_cross[1] + 1,  # vhelix cross at the 3'
#                                  vh_order[start_vhel_cross[2]][0], start_vhel_cross[3] ]
#
#     mesh[end_vhel_cross[3] ] = [vh_order[end_vhel_cross[0]][0], end_vhel_cross[1],  # vhelix cross at 5'
#                                     vh_order[end_vhel_cross[2]][0], end_vhel_cross[3] - 1]
#
#     helix_0_nicks = []
#     # for x in staples_nicks['bundle_' + str(edge)]:
#     #     if 0 in x:  # we are in helix_0
#     #         helix_0_nicks.append(x)
#     #
#     # for x in helix_0_nicks:
#     #     mesh[x[1] - 1] = [-1, -1, id_num, x[1] - 2]
#     #     mesh[x[1]] = [id_num, x[1] + 1, -1, -1]
#
#     #print(mesh[28:])
#
#     return mesh

# def get_virt_stap_double_edge(max_length, staples_crossovers,staples_nicks, id_num, edge, bundles):
#
#     mesh = [[-1, -1, -1, -1] for x in range(max_length)]
#
#     id_num_virt = id_num + 3
#     helix_num = 3
#     max_staple_length = (bundles['bundle_' + str(edge)][0]['helix_' + str(helix_num)]['length']) + 30
#
#     for x in range(31, max_staple_length - 1):
#         mesh[x] = [id_num_virt, x - 1, id_num_virt, x + 1]
#
#     helix_3_crosses = []
#
#     for n in staples_crossovers['bundle_' + str(edge)]:
#         if helix_num in n and 0 in n:
#             helix_3_crosses.append(n)
#
#     for x in helix_3_crosses:
#         if helix_3_crosses.index(x) % 2 != 0:
#             mesh[x[1] - 1] = [id_num_virt, x[1] - 2, x[0] + id_num, x[3] - 1]
#         else:
#             mesh[x[1] + 1] = [x[2] + id_num, x[1] + 1, id_num_virt, x[1] + 2]
#
#     helix_3_nicks = []
#
#     for x in staples_nicks['bundle_' + str(edge)]:
#         if 3 in x:  # we are in helix_3
#             helix_3_nicks.append(x)
#
#     for x in helix_3_nicks:
#         mesh[x[1] + 1] = [-1, -1, id_num_virt, x[1] + 2]
#         mesh[x[1]] = [id_num_virt, x[1] - 1, -1, -1]
#
#     mesh[30] = [-1, -1, id_num_virt, 31]
#     mesh[max_staple_length - 1] = [id_num_virt, max_staple_length - 1, -1, -1]
#
#     return mesh

def get_scaff_reinforced_vhelix_BH(max_length, vh_order, id_num, edge, scaff_cross_vhelix,bundles):   #reinforced scaffold for the first of the double edges

    mesh = [[-1,-1,-1,-1] for x in range(max_length)]

    vhelix_cross_finder = list(item[2] for item in scaff_cross_vhelix)  # list to find the start of the crossover
    vhelix_end_finder = list(item[0] for item in scaff_cross_vhelix)  # list to find the end of the crossover

    for n in vhelix_cross_finder:  # find connections starting from particular edge
        if n == edge:
            index = (vhelix_cross_finder.index(n))
            helix_start = (scaff_cross_vhelix[index])

    for n in vhelix_end_finder:  # find connections ending in particular edge
        if n == edge:
            index_1 = (vhelix_end_finder.index(n))
            helix_end = (scaff_cross_vhelix[index_1])

    helix_start_2 = [vh_order[helix_start[0]][0], helix_start[1], id_num,
                     helix_start[3] + 1]  # +1 because cadnano seems to be abit weird, not sure it will work

    mesh[helix_start_2[3]-1] = helix_start_2

    for x in range(31, helix_end[1]-1):
        mesh[x] = [id_num, x - 1, id_num, x + 1]
        # helix_start[3] = helix_start[3] + 1

    helix_end_2 = [id_num, helix_end[1] - 2, vh_order[helix_end[2]][0], helix_end[3]]
    mesh[helix_end[1]-1] = helix_end_2

    loop_cross_1 = (bundles['bundle_' + str(edge)][1][0])
    loop_cross_2 = (bundles['bundle_' + str(edge)][1][1])

    loop_cross_1_2 = [vh_order[edge][0], loop_cross_1[1] - 1, vh_order[edge][1], loop_cross_1[3]]
    loop_cross_2_2 = [vh_order[edge][1], loop_cross_2[1], vh_order[edge][0], loop_cross_2[3] + 1]
    mesh[loop_cross_1[1]] = loop_cross_1_2
    mesh[loop_cross_2[1]] = loop_cross_2_2

    #print(scaff_cross_vhelix)
    #print(mesh[29:])
    return mesh

def get_stap_vhelix_reinforced_BH(max_length, staples_crossovers, staples_nicks, staple_cross_vhelix, edge, id_num, vh_order, double_edges_coupled): #staple for helix_0

    mesh = [[-1,-1,-1,-1] for x in range(max_length)]

    helix_num = 0

    for x in staple_cross_vhelix:
        if edge == x[0] or edge == x[2]:
            if x[2] == edge:
                max_length_staples = x[3] #+30
                end_vhel_cross = x
            else:
                start_vhel_cross = x

    helix_0_crosses = []


    for x in range(start_vhel_cross[1], max_length_staples):    #starts from when the first staple starts
        mesh[x] = [id_num, x+1, id_num, x-1]

    for x in double_edges_coupled:
        if edge in x:
            double_edge_v_scaff = (vh_order[x[1]][0])

    for n in staples_crossovers['bundle_'+str(edge)]:
        if helix_num in n:
            helix_0_crosses.append(n)

    for x in helix_0_crosses:

        if helix_0_crosses.index(x)%2 != 0:
            if x[2] == 3:
                mesh[x[1]-1] = [double_edge_v_scaff, x[1]-1, x[0] +id_num , x[3]-2]
            else:
                mesh[x[1]-1] = [x[2]+id_num, x[1]-1, x[0] +id_num , x[3]-2]
            #mesh[x[1]-1] = [x[2]+id_num, x[1]-1, x[0] +id_num , x[3]-2]
        else:
            if x[0] ==3:
                mesh[x[1]+1] = [x[2]+id_num, x[1]+2, double_edge_v_scaff, x[1]+1]
            else:
                mesh[x[1]+1] = [x[2]+id_num, x[1]+2, x[0]+id_num, x[1]+1]

    mesh[start_vhel_cross[1]] = [vh_order[start_vhel_cross[0]][0], start_vhel_cross[1] +1,    #vhelix cross at the 3'
                                 vh_order[start_vhel_cross[2]][0], start_vhel_cross[3]]

    mesh[end_vhel_cross[3]] = [vh_order[end_vhel_cross[0]][0], end_vhel_cross[1],       #vhelix cross at 5'
                                 vh_order[end_vhel_cross[2]][0], end_vhel_cross[3] -1]



    helix_0_nicks = []
    for x in staples_nicks['bundle_'+str(edge)]:
        if 0 in x:                                      #we are in helix_0
            helix_0_nicks.append(x)

    for x in helix_0_nicks:
        mesh[x[1]-1] = [-1, -1, id_num, x[1]-2]
        mesh[x[1]] = [id_num, x[1] +1, -1, -1]

    #print(mesh[27:])
    return mesh

def get_scaff_vhelix_double_edge_BH(max_length, vh_order, id_num, edge, scaff_cross_vhelix):   #non reinforced scaffold for double edges, second one in the couple

    mesh = [[-1,-1,-1,-1] for x in range(max_length)]

    vhelix_cross_finder = list(item[2] for item in scaff_cross_vhelix)  # list to find the start of the crossover
    vhelix_end_finder = list(item[0] for item in scaff_cross_vhelix)  # list to find the end of the crossover

    for n in vhelix_cross_finder:  # find connections starting from particular edge
        if n == edge:
            index = (vhelix_cross_finder.index(n))
            helix_start = (scaff_cross_vhelix[index])

    for n in vhelix_end_finder:  # find connections ending in particular edge
        if n == edge:
            index_1 = (vhelix_end_finder.index(n))
            helix_end = (scaff_cross_vhelix[index_1])

    helix_start_2 = [vh_order[helix_start[0]][0], helix_start[1] -1, id_num,
                     helix_start[3] + 1]  # +1 because cadnano seems to be abit weird, not sure it will work

    helix_start_2 =[vh_order[helix_start[0]][0], helix_start[1]-1, id_num, helix_start[3]-1]

    mesh[helix_start_2[3]+1] = helix_start_2

    for x in range(30, helix_start_2[3]+1):
        mesh[x] = [id_num, x + 1, id_num, x - 1]
        # helix_start[3] = helix_start[3] + 1

    helix_end_2 = [id_num, 31, vh_order[helix_end[2]][0], helix_end[3]]
    #helix_end_2 = [vh_order[helix_end[2]][0], helix_end[3]+1, id_num, 30]
    mesh[helix_end_2[1]-1] = helix_end_2

    return mesh

def get_stap_double_edge_BH(max_length, staple_cross_vhelix, edge, id_num, vh_order, lengths_list, double_edges_coupled, staples_crossovers, staples_nicks):

    mesh = [[-1, -1, -1, -1] for x in range(max_length)]

    helix_num = 3

    for x in staple_cross_vhelix:
        if edge == x[0] or edge == x[2]:
            if x[2] == edge:
                max_length_staples = x[3] #+ 30
                start_vhel_cross = x
            else:
                end_vhel_cross = x

    for x in range(1, len(lengths_list)+1):
        if edge == x:
            max_length_staples = lengths_list[x-1] + 30
    helix_3_crosses = []

    for x in range(start_vhel_cross[1], max_length_staples):  # starts from when the first staple starts
        mesh[x] = [id_num, x - 1, id_num, x + 1]


    for x in double_edges_coupled:
        if edge in x:
            double_edge_first = (vh_order[x[0]][0])
            double_edge_third = (vh_order[x[0]][2])

    for n in staples_crossovers['bundle_' + str(edge)]:
        if helix_num in n:
            helix_3_crosses.append(n)
    for x in helix_3_crosses:
        if helix_3_crosses.index(x) % 2 != 0:
            if x[0] == 0:
                mesh[x[1] - 1] = [id_num, x[3] - 2, double_edge_first, x[1] - 1]
            else:
                mesh[x[1] - 1] = [id_num, x[3] - 2, double_edge_third, x[1] - 1]
        else:
            if x[2] == 0:
                mesh[x[1] + 1] = [double_edge_first, x[1] + 1, id_num, x[1] + 2]
            else:
                mesh[x[1] + 1] = [double_edge_third, x[1] + 1,  id_num, x[1] + 2]

    mesh[start_vhel_cross[3]] = [vh_order[start_vhel_cross[0]][0], start_vhel_cross[1],  # vhelix cross at the 5', inverted compared to before
                                 vh_order[start_vhel_cross[2]][0], start_vhel_cross[3]+1 ]

    mesh[end_vhel_cross[1]] = [vh_order[end_vhel_cross[0]][0], end_vhel_cross[1]-1,  # vhelix cross at 3', inverted compared to before
                                    vh_order[end_vhel_cross[2]][0], end_vhel_cross[3]]


    helix_0_nicks = []
    for x in staples_nicks['bundle_' + str(edge)]:
        if 0 in x:  # we are in helix_0
            helix_0_nicks.append(x)

    for x in helix_0_nicks:
        mesh[x[1] + 1] = [-1, -1, id_num, x[1] + 2]
        mesh[x[1]] = [id_num, x[1] - 1, -1, -1]

    #print(str(edge)+'st' + str(mesh[28:]))

    return mesh

def get_stap_loop1_BH(max_length,staples_crossovers, staples_nicks, edge, id_num, bundles):

    mesh = [[-1, -1, -1, -1] for x in range(max_length)]
    id_num_loop_1 = id_num + 1
    helix_num = 1
    max_staple_length = (bundles['bundle_'+str(edge)][0]['helix_'+str(helix_num)]['length'])+39

    for x in range (41, max_staple_length):
        mesh[x] = [id_num_loop_1, x-1, id_num_loop_1, x+1]

    helix_1_crosses = []

    for n in staples_crossovers['bundle_' + str(edge)]:
        if helix_num in n:
            #print(n)
            helix_1_crosses.append(n)


    for x in helix_1_crosses:
        if helix_1_crosses.index(x) % 2 != 0:
            mesh[x[1] - 1] = [id_num_loop_1, x[1] - 2, x[0]+id_num, x[3] - 1]
        else:
            mesh[x[1] + 1] = [x[2]+id_num, x[1] + 1, id_num_loop_1, x[1] + 2]


    helix_1_nicks = []
    for x in staples_nicks['bundle_' + str(edge)]:
        if 1 in x:  # we are in helix_1
            helix_1_nicks.append(x)

    for x in helix_1_nicks:
        mesh[x[1]+1] = [-1, -1, id_num_loop_1, x[1]+2]
        mesh[x[1]] = [id_num_loop_1, x[1]-1, -1, -1]

    mesh[40] = [-1, -1, id_num_loop_1, 41]
    mesh[max_staple_length] = [id_num_loop_1, max_staple_length-1, -1, -1]

    return mesh

def get_stap_loop2_BH(max_length,staples_crossovers, staples_nicks, edge, id_num, bundles, double_edges_coupled, vh_order):

    mesh = [[-1, -1, -1, -1] for x in range(max_length)]
    id_num_loop_2 = id_num + 2
    helix_num = 2
    max_staple_length = (bundles['bundle_'+str(edge)][0]['helix_'+str(helix_num)]['length'])+39

    for x in range (41, max_staple_length):
        mesh[x] = [id_num_loop_2, x+1, id_num_loop_2, x-1]

    helix_2_crosses = []

    for x in double_edges_coupled:
        if edge in x:
            double_edge_v_scaff = (vh_order[x[1]][0])

    for n in staples_crossovers['bundle_' + str(edge)]:
        if helix_num in n:
            helix_2_crosses.append(n)

    for x in helix_2_crosses:
         if helix_2_crosses.index(x) % 2 != 0:
            if x[2] == 3:
                mesh[x[1] - 1] = [double_edge_v_scaff, x[1] -1, id_num_loop_2, x[3] -2]
            else:
                mesh[x[1] - 1] = [x[2] + id_num, x[1] - 1, id_num_loop_2, x[3] - 2]
         else:
            if x[0] ==3:
                mesh[x[1] + 1] = [id_num_loop_2, x[1] + 2, double_edge_v_scaff, x[1] + 1]
            else:
                mesh[x[1] + 1] = [id_num_loop_2, x[1] + 2, x[0]+id_num, x[1] + 1]


    helix_2_nicks = []

    for x in staples_nicks['bundle_' + str(edge)]:
        if 2 in x:  # we are in helix_2
            helix_2_nicks.append(x)

    for x in helix_2_nicks:
        mesh[x[1] - 1] = [-1, -1, id_num_loop_2, x[1] - 1]
        mesh[x[1]] = [id_num_loop_2, x[1] +2, -1, -1]

    mesh[40] = [id_num_loop_2, 41, -1, -1]
    mesh[max_staple_length] = [-1, -1, id_num_loop_2, max_staple_length-1]

    #print(mesh[39:])
    return mesh

def get_colors(staples_nicks, edge, helix_num):    #get colour list for staples. singular for each edge and id_num

    #if helix_num != 0:
    #    i=i+1

    helix_nicks = []

    for x in staples_nicks['bundle_' + str(edge)]:
        if helix_num in x:
            helix_nicks.append(x[1])

    colors_list = [13369344, 12060012, 29184, 8947848, 12060012,
                   243362, 16225054, 12060012,7536862]

    colors = []

    for staple in range(len(helix_nicks)):
        colors.append([helix_nicks[staple]-1, random.choice(colors_list)])

    return colors

def scaffold_length(bundles):

    scaffold_length = 0

    for bundle in bundles.values():
        length_0 = (bundle[0]['helix_0']['length'])
        length_1 = (bundle[0]['helix_1']['length'])
        length_2 = (bundle[0]['helix_2']['length'])
        scaffold_length = scaffold_length + length_0 + length_1 + length_2

    print('The scaffold length is ', scaffold_length)

def insert_deletions(bundles, edge, skips):

    length_loop = bundles['bundle_' + str(edge)][0]['helix_1']['length']
    add_del = False

    if length_loop >= 136:
        add_del = True

    if add_del == True:
        del_number = round(length_loop/48)
        for number in range(0, del_number):
            skips[51 + number*48] = -1

    return skips

def deletions_number(bundles, edge):

    length_loop = bundles['bundle_' + str(edge)][0]['helix_1']['length']
    add_del = False

    del_number = 0

    if length_loop >= 136:
        add_del = True


    if add_del == True:
        del_number = round(length_loop / 48)

    return del_number

def GenerateJson(filenameInput, to_reinforce):
    data, fwd_helix_connections, rev_helix_connections = open_rpoly(filenameInput + ".rpoly")
    n_trail_list = open_ntrail(filenameInput + ".ntrail")
    faces_list = open_ply(filenameInput + ".ply")

    double_vertices = False

    filename = filenameInput + ".json"
    #folder = os.path.join('/', filename)

    #data, fwd_helix_connections, rev_helix_connections = open_rpoly("Triangular_BH_4_shorter_ply_3.rpoly")
    #n_trail_list = open_ntrail("Triangular_BH_4_shorter_ply.ntrail")
    #faces_list = open_ply("Triangular_BH_4_shorter_ply_2.ply")

    # data, fwd_helix_connections, rev_helix_connections = open_rpoly("triple_hex_flat_2.rpoly")
    # n_trail_list = open_ntrail("triple_hex_flat_1.ntrail")
    # faces_list = open_ply("triple_hex_flat_1.ply")

    example_edge = 6

    #edge=example_edge

    vert_to_helix_number = vertices_from_ntrail(n_trail_list)

    tot_edges_vectors_dict, position_list, lengths_list, vec_list, new_position_list, vec2_list = vectors_from_rpoly(data)

    #not_reinforce = [3, 9, 19, 30, 40, 43, 44]
    not_reinforce = []
    #to_reinforce = [edge for edge in to_reinforce if edge not in not_reinforce]
    #print(to_reinforce)
    bundles_dictionary = {}

    for element in range(1, len(lengths_list)+1):
        edge=element
        bundle = reinforced_bundle_creation(edge, vert_to_helix_number, faces_list, vec_list, tot_edges_vectors_dict, position_list, lengths_list)
        bundles = bundles_colletion(bundle, bundles_dictionary, edge)


    helices_3d_dict = helix_3d(vec_list, new_position_list, vec2_list, lengths_list)
    gaps_dimensions = (auto_fill_strand_gaps(fwd_helix_connections, rev_helix_connections, helices_3d_dict))

    staples_crossovers, staples_nicks = (add_staples(bundles, to_reinforce))

    double_edges, double_edges_coupled = double_edges_counter(vert_to_helix_number)

    bundles,gaps_filled_scaffold_lengths, scaff_cross_vhelix = scaff_vhelix_crossovers(lengths_list, gaps_dimensions,
                                                                                       bundles, double_edges_coupled)

    staple_cross_vhelix = staple_vhelix_crossovers(gaps_dimensions, bundles, lengths_list, double_edges_coupled)

    filename_excel = str(filename[:-5]) + '_virt_scaff.xlsx'

    deletions_list = output_cadnano_2(filename, data, to_reinforce, gaps_filled_scaffold_lengths, bundles, scaff_cross_vhelix, double_vertices, double_edges, double_edges_coupled, staples_crossovers, staples_nicks, staple_cross_vhelix,lengths_list, faces_list, vert_to_helix_number)

    edges_number = len(bundles)

    faces_to_edges_list = double_vertices_faces(faces_list, vert_to_helix_number, to_reinforce)

    save_workbook(edges_number, bundles, filename_excel, deletions_list, double_vertices, faces_to_edges_list, to_reinforce, double_edges)

if __name__ == "__main__":

    filenameInput = "triple_hex_flat_1"
    to_reinforce = [x for x in range (1, 38)]
    GenerateJson(filenameInput, to_reinforce)
