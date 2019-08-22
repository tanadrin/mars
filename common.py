import pyglet
import geometry
from math import sin
from math import cos
from geometry import pi
from PIL import Image
import itertools
import numpy

hmap_x = 2048
hmap_y = 1024
heightmap = Image.open('mars_mola_bumpmap_2048x1024.jpg')
colormap = Image.open('mars_colormap.jpg')
colormap.convert('RGB')

class Location:
    def __init__(self):
        self.name = None
        self.lat=None
        self.long=None
        self.theta = None
        self.phi = None
        self.rho = None
        self.x = None
        self.y = None
        self.z = None
        self.stereox = None
        self.stereoy = None
        self.v = None
        self.c = None
        self.index = None
        self.neighbors = []
        self.vertices = []
        self.pd = None # Flattened array of vertex values for use w Pyglet
        self.coldata = None # Tuple used to render color of loc w Pyglet
    # Generates location from cartesian point data
    def loc_from_point(self,name,x,y,z):
        self.name = name
        self.x = x
        self.y = y
        self.z = z
        self.rho,self.theta,self.phi=geometry.cart_to_sphere(x,y,z)
        self.stereox,self.stereoy = self.stereo()
        self.lat,self.long=geometry.latlong(self.theta,self.phi)
        self.get_rho()
        self.get_rgb()
    # Generates a location from spherical point data
    def loc_from_point2(self,name,rho,theta,phi):
        self.name = name
        self.rho,self.theta,self.phi=rho,theta,phi
        self.x,self.y,self.z=geometry.sphere_to_cart(rho,theta,phi)
        self.stereox,self.stereoy = self.stereo()
        self.lat,self.long=geometry.latlong(self.theta,self.phi)
        self.get_rho()
        self.get_rgb()
    # Generates a location from a line specified in a file
    # using lat and long
    def loc_from_file(self,name,latitude,longitude):
        self.name = name
        self.lat,self.long(latitude,longitude)
        self.theta = float(latitude)*pi/180 # deg -> rad
        self.phi = float(longitude)*pi/180 # deg -> rad
        self.rho = 1
        self.x,self.y,self.z=geometry.sphere_to_cart(self.rho,self.theta,self.phi)
        self.stereox,self.stereoy = self.stereo()
        self.get_rho()
        self.get_rgb()
    # Stores stereographic projection of the location's cartesian coordinates
    def stereo(self):
        return geometry.stereo_points([[self.x,self.y,self.z]])[0]
    # Adjusts the radius of the point based on the elevation of a corresponding
    # coordinate on the heightmap image provided.
    def get_rho(self):
        img_x,img_y = geometry.latlong_to_xy(self.theta,self.phi,hmap_x,hmap_y)
        self.v=heightmap.getpixel((img_x,img_y))
        self.rho = (self.v/255)*0.1+1 # magic number, fix later
        self.x,self.y,self.z=geometry.sphere_to_cart(self.rho,self.theta,self.phi)
    # Gets the color we want to use for the location based on the RGB value
    # of the pixel at the corresponding coordinate of the colormap provided.
    def get_rgb(self):
        img_x,img_y = geometry.latlong_to_xy(self.theta,self.phi,hmap_x,hmap_y)
        self.c=colormap.getpixel((img_x,img_y))
    # Stores color data in a way that's easy to pass to OpenGL.
    def set_coldata(self):
        coldata = []
        i = len(self.vertices)
        for j in range(0,i):
            for k in self.c:
                coldata.append(k)
        self.coldata = tuple(coldata)

# Class for storing and getting information about the whole set of locations a
# little more easily.        
class LocationList:
    def __init__(self,locationfile = None, pointlist = None):
        self.locations = []
        if locationfile != None and pointlist == None:
            for line in locationfile:
                name,diameter,latitude,longitude = line.split(",")
                loc = Location()
                loc.loc_from_file(name,latitude,longitude)
                self.locations.append(loc)
        elif locationfile == None and pointlist != None:
            for point in pointlist:
                loc = Location()
                name = 'loc'
                loc.loc_from_point2(name,point[0],point[1],point[2])
                self.locations.append(loc)
        else:
            raise ValueError('Specify either a list of locations via file OR a list of points, but not both')
    # Returns a list of points corresponding to the 3D coordinates of each
    # location
    def sphere_points(self):
        sphere_points = []
        for location in self.locations:
            sphere_points.append([location.x,location.y,location.z])
        return sphere_points
    # Same but for points on the stereographic projection.
    def stereo_points(self):
        stereo_points = []
        for location in self.locations:
            stereo_points.append([location.stereox,location.stereoy])
        return stereo_points
        
    # Find the neighbor of each point and the Voronoi region associated w/ each point,
    # and order the vertices of the Voronoi region so that they render correctly w/ OpenGL.
    def build_graph(self,points,delaunay):
        
        indices,indptr=delaunay.vertex_neighbor_vertices # Scipy docs: Neighboring vertices of vertices.
        # Tuple of two ndarrays of int: (indices, indptr). The indices of neighboring vertices of vertex 
        # k are indptr[indices[k]:indices[k+1]].

        for location in self.locations:
        
            location.index = self.locations.index(location)
            
            # Each neighbor point for the location
            location.neighbors = indptr[indices[location.index]:indices[location.index+1]]
            
            # Get the actual values of the neighbor points, then sort them by rotation around
            # the reference point; this will mean the vertexes we calculate in the next step
            # will already be in order.
            neighbor_points = []
            for i in location.neighbors:
                neighbor_points.append(points[i])
            sorted_neighbors = geometry.order_points_3d(neighbor_points)
            
            # The Voronoi region vertices are the same as the centroids of the set of triangles
            # touching each vertex.
            centroids = []
            for j in range(0,len(sorted_neighbors)):
                p1 = sorted_neighbors[j]
                p2 = sorted_neighbors[(j+1)%len(sorted_neighbors)]
                centroid = geometry.get_centroid2([points[location.index],p1,p2])
                centroids.append(centroid)
                
            location.vertices = centroids
                
            # Convert the components of the vertices into a flat tuple for use with
            # OpenGL. Will speed up rendering later.
            pd = []
            for vertex in location.vertices:
                for value in vertex:
                    pd.append(value)
            location.pd = tuple(pd)
            location.set_coldata()