import numpy as np
from load_files import open_ply
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

file = "import_files/Cuboctahedron/Cuboctahedron.ply"
# file = "import_files/17_triangular_bipyramid.ply"
number_vertices, vertices_list, number_face, faces_list = open_ply(file)
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

for i in range(number_face):
    for j in range(1,faces_list[i][0]+1):
        if j == len(faces_list[0])-1:
            point1 = faces_list[i,j]
            point2 = faces_list[i,1]
        else:            
            point1 = faces_list[i,j]
            point2 = faces_list[i,j+1]

        xx = [vertices_list[point1,0], vertices_list[point2,0]]
        yy = [vertices_list[point1,1], vertices_list[point2,1]]
        zz = [vertices_list[point1,2], vertices_list[point2,2]]
        ax.plot(xx,yy,zz,'r')


ax.scatter(vertices_list[:,0],vertices_list[:,1],vertices_list[:,2])

for i in range(number_vertices):
    ax.text(vertices_list[i,0],vertices_list[i,1],vertices_list[i,2], s=i)
plt.show()