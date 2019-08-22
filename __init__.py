import pyglet
import math
import geometry
import random
import numpy
import common
from pyglet.gl import *
import scipy

global camera_set # used for camera initialization
n=3000 # Number of locations to generate
w=500 # 1/2 viewport width
camera_set = None

game_window = pyglet.window.Window(2*w,2*w)

def draw_loc(loc,loc_batch):
    pyglet.graphics.draw(loc.pdl,pyglet.gl.GL_POLYGON,('v3f',loc.pd),('c3B',loc.coldata))

class ProvinceGroup(pyglet.graphics.Group):
    def set_state(self):
        pass
    
@game_window.event
def on_draw():
    global camera_set
    if camera_set == None:
        glTranslatef(0,0,-5)
        glRotatef(90,1,0,0)
        glPushMatrix()
        camera_set = 1
    game_window.clear()
    loc_batch = pyglet.graphics.Batch()
    g=pyglet.graphics.Group()
    
    for loc in location_list.locations:
        loc_batch.add(loc.pdl,pyglet.gl.GL_TRIANGLE_FAN, g,('v3f',loc.pd),('c3B',loc.coldata))
    loc_batch.draw()

@game_window.event    
def on_key_press(key,modifiers):
    if key == pyglet.window.key.C:
        glPopMatrix()
        glPushMatrix()
    if key == pyglet.window.key.SPACE:
        glPushMatrix()
    
@game_window.event
def on_text_motion(motion):
    if motion == pyglet.window.key.MOTION_LEFT:
        glRotatef(15,0,0,1)
    if motion == pyglet.window.key.MOTION_RIGHT:
        glRotatef(-15,0,0,1)
    if motion == pyglet.window.key.MOTION_UP:
        glRotatef(15,0,1,0)
    if motion == pyglet.window.key.MOTION_DOWN:
        glRotatef(-15,0,1,0)
        
@game_window.event
def on_resize(width,height):
    glViewport(0,0,width,height)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(90, width / float(height), 0.0001, 9000)
    glMatrixMode(GL_MODELVIEW)
    return pyglet.event.EVENT_HANDLED
        
if __name__ == '__main__':
    pyglet.graphics.glEnable(pyglet.graphics.GL_DEPTH_TEST)
    location_list = common.LocationList(pointlist=geometry.spiral_points2(n))
    points = location_list.sphere_points()
    stereo_points = location_list.stereo_points()
    delaunay = scipy.spatial.Delaunay(stereo_points)
    location_list.build_graph(points,delaunay)   
    pyglet.app.run()
    