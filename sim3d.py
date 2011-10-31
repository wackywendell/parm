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
import subprocess

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
    """Create a new window for 3D videos.
    
    A window can then be turned into an application using the run() command.
    
    Use the "run" method to create an application with event loop, with 
    keypresses stored in the 'keys' dict.
    """
    def __init__(self, xsize=300, ysize=0, outfile=None, rot_init=(30,55), 
                dAngle=30, fps=20, grid=True, bg=(0,0,0),
                rotate=False):
        """__init__ args:
        xsize, ysize: size of the window, in pixels
        outfile: an open file descriptor in which to dump raw video (see Vidwriter)
        rot_init: initial rotation of the cube
        dAngle: amount to rotate cube on key press
        fps: expected frames per second (to get rotation looking right)
        grid: draw the red box around simulation
        bg: background color, in rgb floats (e.g. (1,1,1) is white, balck default)"""
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
        self.outfile = outfile
        self.dAngle = dAngle
        self.fps = fps
        self.grid=grid
        self.bg=bg
        self.rotate = rotate
        if outfile is not None:
            self.buf = ctypes.create_string_buffer(3*(ysize*xsize))
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
        """Draw a sphere at location mvec, with a radius.
        
        slices and stacks are GL things; a sphere is drawn as triangles,
        and the more slices and stacks the smoother the sphere, but
        the longer to render.
        
        Call setcol() before for colors."""
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
        """Draw a line from points in lst.
        
        Thickness in pixels."""
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
        """Set color for next element. col should be an RGB tuple of floats,
        e.g. (1,1,1) for white. Alpha is transparency.
        
        Use no arguments for a random color."""
        try:
            col = tuple(col) + (alpha,)
        except:
            col = randcol() + (alpha,)
        glColor4f(*col)
    
    def clear(self):
        """Clear the window."""
        self.window.clear()
        glColor3f(1, 0, 0)
        if self.grid: self.drawgrid()
        
    def caption(self, s):
        """Set the window title. Note that this will not show up in any
        output videos."""
        self.window.set_caption(s)
    
    def output(self):
        glReadPixels(0,0,self.sizex,self.sizey, GL_RGB, GL_UNSIGNED_BYTE, self.buf)
        self.outfile.write(self.buf.raw)
    
    def run(self, updatefunc):
        """Run as an application.
        
        Opens a window and calls updatefunc(timestep), where timestep
        should be roughly 1/fps. Updatefunc should raise a StopIteration
        exception when the simulation is over."""
        def update(dt):
            self.clear()
            try:
                updatefunc(dt)
            except StopIteration:
                pyglet.app.exit()
            if self.outfile is not None:
                self.output()
            self.keyhandler(dt)
        pyglet.clock.schedule_interval(update, 1/float(self.fps))
        pyglet.app.run()

class Vidwriter:
    """
    Run a subprocess in which to dump video.
    
    Use a Window instance with this.inputfile as the Window's outfile
    to get output.
    
    Call 'close' when finished to close the inputfile and wait for the
    video file to finish writing.
    """
    
    args = "/usr/bin/ffmpeg -f rawvideo -pix_fmt rgb24 -s %dx%d -y -an -vf vflip -r %d -i /dev/stdin -b 2M -bt 1M -vcodec libx264 %s"
    # could also be
    # args = "mencoder /dev/stdin -demuxer rawvideo -rawvideo format='rgba':w=%d:h=%d:fps=%d -flip -vf pp=al -ovc lavc -o %s"
    
    
    def __init__(self, fname, sz, fps=20):
        args = self.args % (sz,sz,fps,fname)
        print("Starting", args)
        self.proc = proc = subprocess.Popen(args, stdin=subprocess.PIPE, 
                    stderr = open('/dev/null','wb'), shell=True)
        self.inputfile = proc.stdin
    
    def close(self):
        """Finish writing video, and wait for subprocess to close."""
        self.inputfile.close()
        return self.proc.wait()
    
    def __del__(self):
        self.close()
