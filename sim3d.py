#!/usr/bin/python2
# python2 movie5.py | mencoder /dev/stdin -demuxer rawvideo -rawvideo format='rgba':w=600:h=600:fps=10 -flip -vf pp=al -ovc lavc -o test.avi
# python2 movie5.py | ffmpeg -f rawvideo -pix_fmt rgb24 -s 500x500 -y -an -vf vflip -r 20 -i /dev/stdin -b 2M -bt 1M -vcodec libx264 test.mp4

from __future__ import print_function
from sys import stdout, stderr
import pyglet
from pyglet.gl import *
import ctypes
import time, math, random
from pyglet import font
from pyglet.window import key
from itertools import islice
import math
from colorsys import hsv_to_rgb as hsv2rgb

class KeyStateHandler(dict):
    mods = (key.MOD_SHIFT, key.MOD_CTRL, key.MOD_ALT, key.MOD_CAPSLOCK, 
    key.MOD_NUMLOCK, key.MOD_WINDOWS, key.MOD_COMMAND, key.MOD_OPTION, 
    key.MOD_SCROLLLOCK)
    
    def on_key_press(self, symbol, modifiers):
        self[symbol] = True
        for m in self.mods:
            if modifiers & m:
                self[m] = True
        #~ print("handled", symbol)
        
    def on_key_release(self, symbol, modifiers):
        if symbol in self:
            self[symbol] = False
        for m in self.mods:
            if modifiers & m:
                self[m] = False

    def __getitem__(self, key):
        return self.get(key, False)

def GLvec(*args):
    return (GLfloat * len(args))(*args) 

def randcol():
    return hsv2rgb(360*random.random(),random.triangular(0,1,1),
                                                    random.triangular())

def towardcol((r,g,b), weight=10):
    (r2,g2,b2) = hsv2rgb(360*random.random(),random.triangular(0,1,1),
                                                    random.triangular())
    r = (r*weight + r2) / (weight + 1)
    g = (g*weight + g2) / (weight + 1)
    b = (b*weight + b2) / (weight + 1)

class Window:
    def __init__(self, xsize=300, ysize=0, rot_init=(30,55), 
                output=False, dAngle=30, fps=20, grid=True, bg=(0,0,0),
                rotate=False):
        if ysize <= 0: ysize = xsize
        self.sizex = xsize
        self.sizey = ysize
        try:
            # Try and create a window with multisampling (antialiasing)
            self.config = Config(sample_buffers=1, samples=4, 
                            depth_size=16, double_buffer=True)
            self.window = pyglet.window.Window(xsize,ysize, config=self.config)
        except pyglet.window.NoSuchConfigException:
            # Fall back to no multisampling for old hardware
            self.window = pyglet.window.Window(xsize,ysize)
            #~ print("fallback", file=stderr)
        
        self.keys = KeyStateHandler()
        self.window.push_handlers(self.keys)
        self.on_resize = self.window.event(self.on_resize)
        self.rot = rot_init
        self.output = output
        self.dAngle = dAngle
        self.fps = fps
        self.grid=grid
        self.bg=bg
        self.rotate = rotate
        if output:
            self.buf = ctypes.create_string_buffer(3*(w*h))
        self.setup()
        

    def on_resize(self, width, height):
        # Override the default on_resize handler to create a 3D projection
        self.sizex, self.sizey = width, height
        glViewport(0, 0, width, height)
        glMatrixMode(GL_PROJECTION)
        gluPerspective(40., width / float(height), .1, 1000.)
        glMatrixMode(GL_MODELVIEW)
        self.set_rot(*self.rot)
        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)
        posL = .8
        glLightfv(GL_LIGHT0, GL_POSITION, GLvec(posL,posL,posL,0))
        #~ glLightfv(GL_LIGHT0, GL_DIFFUSE, GLvec(1, 1, 1, 1))
        ambL = .5 
        glLightfv(GL_LIGHT0, GL_AMBIENT, GLvec(ambL,ambL,ambL,ambL))
        glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, GLvec(0.5, 0.5, 0.5, 1)) 
        #~ glLineWidth(2)
        return pyglet.event.EVENT_HANDLED
    
    def setup(self):
        # One-time GL setup
        glClearColor(*(self.bg+(0,)))
        glColor3f(1, 0, 0)
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_COLOR_MATERIAL)
        glEnable(GL_BLEND)
        glAlphaFunc(GL_EQUAL, .1)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        #~ glEnable(GL_CULL_FACE)
        self.set_rot(*self.rot)

    def set_rot(self, vrot, hrot):
        glLoadIdentity()
        glTranslatef(-.7, -.4, -2.5)
        glRotatef(vrot, 1, 0, 0)
        glRotatef(hrot, 0, 1, 0)

    def keyhandler(self, dt):
        keys = self.keys
        if key.Q in keys: # if it has ever been pressed
            self.window.close()
            pyglet.app.exit()
        ang = self.dAngle*dt
        if self.output:
            ang = self.dAngle / self.fps
        vvec = (1,0,1)
        hvec = (0,-1,0)
        mvec = (.5,.5,.5)
        imvec = tuple(-x for x in mvec)
        glTranslatef(*mvec)
        if keys[key.MOD_SHIFT]:
            ang *= 3
            #~ print('shift')
        if keys[key.MOD_CTRL]:
            ang *= 5
            #~ print('ctrl')
        if keys[key.LEFT]:
            glRotatef(ang, *hvec)
        if keys[key.RIGHT] or self.rotate:
            glRotatef(-ang, *hvec)
        if keys[key.UP]:
            glRotatef(ang, *vvec)
        if keys[key.DOWN]:
            glRotatef(-ang, *vvec)
        glTranslatef(*imvec)
        if key.SPACE in keys:
            set_rot(30,55)
            opts.rotate = False
            del keys[key.SPACE]
        if key.R in keys:
            opts.rotate = not self.rotate
            del keys[key.R]

    def drawsphere(self,mvec, radius, slices=16,stacks=16):
        sphere = gluNewQuadric()
        gluQuadricDrawStyle( sphere, GLU_FILL)
        gluQuadricNormals( sphere, GLU_SMOOTH)
        gluQuadricOrientation( sphere, GLU_OUTSIDE)
        gluQuadricTexture( sphere, GL_TRUE)
        imvec = tuple(-x for x in mvec)
        glTranslatef(*mvec)
        gluSphere(sphere, radius, slices, stacks)
        glTranslatef(*imvec)
        
    def drawline(self, lst,thickness=None):
        if thickness is not None:
            glLineWidth(thickness)
        pts = [tuple(x) for x in lst]
        N = len(pts)
        pts = [coord for pt in pts for coord in pt]
        pyglet.graphics.draw(N, GL_LINE_STRIP,
            ('v3f', pts))
    
    def drawrod(self, lst,thickness):
        if thickness is not None:
            glLineWidth(thickness)
        pts = [tuple(x) for x in lst]
        N = len(pts)
        pts = [coord for pt in pts for coord in pt]
        pyglet.graphics.draw(N, GL_LINE_STRIP,
            ('v3f', pts))
    
    def drawgrid(self):
        glColor3f(1, 0, 0)
        glLineWidth(1   )
        pyglet.graphics.draw(16, GL_LINE_STRIP,
            ('v3f', (0,0,0,
                     0,1,0,
                     1,1,0,
                     1,0,0,
                     0,0,0,
                     0,0,1,
                     0,1,1,  0,1,0,  0,1,1,
                     1,1,1,  1,1,0,  1,1,1,
                     1,0,1,  1,0,0,  1,0,1,
                     0,0,1)))
    
    def setcol(self,col=0,alpha=1):
        try:
            col = tuple(col) + (alpha,)
        except:
            col = randcol() + (alpha,)
        glColor4f(*col)
    
    def clear(self):
        self.window.clear()
        glColor3f(1, 0, 0)
        if self.grid: self.drawgrid()
        
    def caption(self, s):
        self.window.set_caption(s)
    
    def output(self, f=stdout):
        glReadPixels(0,0,w,h, GL_RGB, GL_UNSIGNED_BYTE, buf)
        f.write(buf.raw)
    
    def run(self, updatefunc, fps=0):
        dt = 0
        if fps>0: dt = 1.0/float(fps)
        def update(dt):
            self.clear()
            try:
                updatefunc(dt)
            except StopIteration:
                pyglet.app.exit()
            self.keyhandler(dt)
        pyglet.clock.schedule_interval(update, dt)
        pyglet.app.run()


#~ w = Window()
#~ col = randcol()
#~ 
#~ def updateb(dt):
    #~ w.setcol(alpha=1)
    #~ w.drawsphere((random.random(),random.random(),random.random()),random.random()*.3)

#~ w.run(updateb, 3)
