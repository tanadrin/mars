from math import sqrt
from math import cos
from math import sin
from math import atan
from math import acos
from math import log
from math import atan2
import random
import scipy.spatial
import numpy

pi = 3.14159

# Converts a cartesian coordinate to a spherical one
def cart_to_sphere(x,y,z):
    rho = sqrt(x**2+y**2+z**2)
    theta = acos(z/rho)
    phi = atan(y/x)
    return rho,theta,phi

# Converts a spherical coordinate to a cartesian one
def sphere_to_cart(rho,theta,phi):
    x = rho*sin(theta)*cos(phi)
    y = rho*sin(theta)*sin(phi)
    z = rho*cos(theta)
    return x,y,z

# Polar to cartesian coordinates    
def polar_to_cart(radial,angular):
    x=radial*cos(angular)
    y=radial*sin(angular)
    return x,y

# Returns coordinates of n points arranged
# roughly evenly spaced on the surface of a sphere.
def spiral_points(n, radius = 1):
    points_list = []
    s = 3.6/sqrt(n)
    dz = 2.0/n
    long = 0
    z = 1-dz/2
    for k in range(0,n-1):
        r = sqrt(1-z*z)
        node_k = (cos(long)*r,sin(long)*r,z)
        z = z-dz
        long = long+s/r
        points_list.append(node_k)
    return points_list

# As the previous function, but returns
# spherical coordinates.
def spiral_points2(n,radius=1):
    points_list = []
    phi_k = 0
    phi_k0 = 0
    for k in range(1,n):
        h_k = -1+(2*(k-1))/(n-1)
        theta_k = acos(h_k)
        if k == 1 or k == n:
            phi_k = 0
        else:
            phi_k = ( phi_k0 + 3.6/sqrt(n) * 1/sqrt(1-h_k**2) ) % (2*pi)
        phi_k0 = phi_k
        points_list.append([radius,theta_k,phi_k])
    return points_list
    
# Takes a set of points (x,y,z) on the 
# sphere and returns the stereographic projection, in 
# cartesian coordinates.
def stereo_points(points):
    new_points = []
    for point in points:
        x,y,z=point[0],point[1],point[2]
        X = x/(1-z)
        Y = y/(1-z)
        new_point = [X,Y]
        new_points.append(new_point)
    return new_points

# Takes points x,y on the plane and returns points
# x,y,z on the sphere.
def stereo_to_sphere(points):
    new_points = []
    for point in points:
        X,Y=point[0],point[1]
        x=X/(1+X**2+Y**2)
        y=Y/(1+X**2+Y**2)
        z=(-1+X**2+Y**2)/(2+2*X**2+2*Y**2)
        new_point=[x,y,z]
        new_points.append(new_point)
    return (new_points)

# Returns the edges of a 2D Voronoi diagram, in Cartesian
# coordinates; see 
# https://docs.scipy.org/doc/scipy-0.18.1/reference/generated/scipy.spatial.Voronoi.html
def stereo_voronoi_edges(points):
    voronoi = scipy.spatial.Voronoi(points)
    vertices = voronoi.vertices # The vertices of the Voronoi edges
    edges = voronoi.ridge_vertices # Pairs of indices of vertices 
    # in the previous array forming each edge
    return vertices, edges

# Returns the simplices of a 2D Delaunay triangulation of
# a set of points; see
# https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.spatial.Delaunay.html    
def stereo_delaunay(points):
    delaunay = scipy.spatial.Delaunay(points)
    simplices = delaunay.simplices # "Indices of the points forming the
    # simplices in the triangulation... oriented counterclockwise."
    return points, simplices

# Apparently the convex hull of points on a sphere is equivalent
# to their Delaunay triangulation, and means there's no gap in the mesh at the
# south pole. See section 5 at:
# https://www.redblobgames.com/x/1842-delaunay-voronoi-sphere/
def convex_hull(points):
    chull = scipy.spatial.ConvexHull(points)
    return chull.points,chull.simplices

# Rotates a point x,y,z in the plane x,y around the origin 
# by an angle given in radians
def rotatexy(point,angle):
    px,py,pz=point[0],point[1],point[2]
    qx=cos(angle)*px-sin(angle)*py
    qy=sin(angle)*px+cos(angle)*py
    qz=pz
    return qx,qy,qz

# Converts radian coordinates to degree coordinates
def latlong(theta,phi):
    return theta*180/pi,phi*180/pi

# Gets XY coordinates on an equirectangular projection
# from spherical coordinates  
def latlong_to_xy(theta,phi,max_x,max_y):
    x = max_x-1-(phi*(max_x-1)/(2*pi))
    y = theta*(max_y-1)/pi
    return int(x),int(y)

# Finds the centroid of a triangular simplex stored as
# indices of points on a list.
def get_centroid(points,simplex):
    p1,p2,p3=points[simplex[0]],points[simplex[1]],points[simplex[2]]
    cx = (p1[0]+p2[0]+p3[0])/3
    cy = (p1[1]+p2[1]+p3[1])/3
    cz = (p1[2]+p2[2]+p3[2])/3
    return (cx,cy,cz)
    
# Finds common centroid of an arbitrary number of points.
def get_centroid2(points):
    cx,cy,cz = 0,0,0
    for point in points:
        cx += point[0]
        cy += point[1]
        cz += point[2]
    cx = cx/len(points)
    cy = cy/len(points)
    cz = cz/len(points)
    return (cx,cy,cz)
    
# Finds common centroid of an arbitrary number of points in 2d.
def get_centroid3(points):
    cx,cy = 0,0
    for point in points:
        cx += point[0]
        cy += point[1]
    cx = cx/len(points)
    cy = cy/len(points)
    return (cx,cy)

# From Stackoverflow    
def order_points_2d(points):
    c = get_centroid3(points)
    points.sort(key=lambda p: atan2(p[1]-c[1],p[0]-c[0]))
    return points

# Adds Cartesian points of arbitrary dimensionality    
def addp(point1,point2):
    sum_components = []
    for i in range(0,len(point1)):
        sum.component.append(point1[i]+point2[i])
    return tsum_components

# Subtracts Cartesian points of arbitrary dimensionality    
def subp(point1,point2):
    sum_components = []
    for i in range(0,len(point1)):
        sum_components.append(point1[i]-point2[i])
    return sum_components

def order_points_3d(points):
    '''
    Per Tumblr user @jadagul:
    Find the cross product of two vectors (here,
    from the centroid to the first two pts in the
    list), to get the normal vector of the plane
    all the points lie on. Write down the matrix R,
    whose columns are these two points and the
    normal vector, respectively. Multiply each point
    by the multiplicative inverse of R (r_inv),
    using numpy.linalg.inv(R). The result should be
    a list of points in the xy plane; you can order
    them, then reverse the process by multiplying
     back by R.
    '''
    c = numpy.array(get_centroid2(points))
    v1,v2 = numpy.array(subp(points[0],c)),numpy.array(subp(points[1],c)) # v1,v2 are
    # vectors from the centroid c to the first two vtces
    n = numpy.cross(v1,v2) # Normal of the plane the pts lie on
    R = numpy.matrix([[points[0][0],points[0][1],points[0][2]],[points[1][0],points[1][1],points[1][2]],[n[0],n[1],n[2]]])
    r_inv = numpy.linalg.inv(R) # Multiplicative inverse of R
    points.sort(key=lambda p: atan2((p*r_inv).getA1()[1]-(c*r_inv).getA1()[1],(p*r_inv).getA1()[0]-(c*r_inv).getA1()[0]))
    return points