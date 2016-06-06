# -*- coding: utf-8 -*-
"""
This example demonstrates many of the 2D plotting capabilities
in pyqtgraph. All of the plots may be panned/scaled by dragging with 
the left/right mouse buttons. Right click on any plot to show a context menu.
"""

from pyqtgraph.Qt import QtGui, QtCore
import numpy as np
import pyqtgraph as pg
import sys
import cPickle

app = QtGui.QApplication([])

win = pg.GraphicsWindow(title="VAF Real Time Viewer")
win.resize(1000,600)
win.setWindowTitle('VAF Real Time Viewer')

# Enable antialiasing for prettier plots
pg.setConfigOptions(antialias=True)

p1 = win.addPlot(title="X centroid")
curve1 = p1.plot(pen=(255,255,0))
bpm1 = p1.plot(pen=None, symbol='s',symbolSize=5,symbolBrush=(100, 100, 255))
cor1 = p1.plot(pen=None, symbol='s',symbolSize=5,symbolBrush=(255, 100, 100))
p1.showGrid(x=True, y=True)
p1.setLabel('left', "X orbit [mm]")
p1.setLabel('bottom', "Z [m]")

p2 = win.addPlot(title="Y centroid")
curve2 = p2.plot(pen=(255,255,0))
bpm2 = p2.plot(pen=None, symbol='s',symbolSize=5,symbolBrush=(100, 100, 255))
cor2 = p2.plot(pen=None, symbol='s',symbolSize=5,symbolBrush=(255, 100, 100))
p2.showGrid(x=True, y=True)
p2.setLabel('left', "Y orbit [mm]")
p2.setLabel('bottom', "Z [m]")

win.nextRow()

p3 = win.addPlot(title="X RMS")
curve3 = p3.plot(pen=(255,255,0))
bpm3 = p3.plot(pen=None, symbol='s',symbolSize=5,symbolBrush=(100, 100, 255))
p3.showGrid(x=True, y=True)
p3.setLabel('left', "X RMS [mm]")
p3.setLabel('bottom', "Z [m]")

p4 = win.addPlot(title="Y RMS")
curve4 = p4.plot(pen=(255,255,0))
bpm4 = p4.plot(pen=None, symbol='s',symbolSize=5,symbolBrush=(100, 100, 255))
p4.showGrid(x=True, y=True)
p4.setLabel('left', "Y RMS [mm]")
p4.setLabel('bottom', "Z [m]")

#ldata = np.zeros(shape=(823,10))
with open('plot.info','r') as f:
    info=cPickle.load(f)

y1bl = []
ptr = 0

def update():
    global curve1, curve2, curve3, curve4
    global curve1o 
    global p1,p2,p3,p4
    global data, ptr
    global x, y1, y2, y3, y4
    global xb, y1b, y2b, y3b, y4b
    global xc, y1c, y2c

    try: 
        ldata = np.loadtxt('ldata.txt')
        if len(ldata) == info[0]:
            x = ldata[:,0]
            y1 = ldata[:,1]
            y2 = ldata[:,2]
            y3 = ldata[:,4]
            y4 = ldata[:,5]
            xb = ldata[info[1],0]
            y1b = ldata[info[1],1]
            y2b = ldata[info[1],2]
            y3b = ldata[info[1],4]
            y4b = ldata[info[1],5]
            xc  = ldata[info[2][::2],0]
            y1c = ldata[info[2][::2],1]
            y2c = ldata[info[2][::2],2]


    except:
        pass


    curve1.setData(x,y1)
    bpm1.setData(xb,y1b)
    cor1.setData(xc,y1c)
    curve2.setData(x,y2)
    bpm2.setData(xb,y2b)
    cor2.setData(xc,y2c)
    curve3.setData(x,y3)
    bpm3.setData(xb,y3b)
    curve4.setData(x,y4)
    bpm4.setData(xb,y4b)
    if ptr == 0 or x[-1] == 0.0:
        p1.enableAutoRange('xy', False)
        p2.enableAutoRange('xy', False)
        p3.enableAutoRange('xy', False)
        p4.enableAutoRange('xy', False)
    ptr += 1

timer = QtCore.QTimer()
timer.timeout.connect(update)
timer.start(50)



'''
p4 = win.addPlot(title="Parametric, grid enabled")
x = ldata[:,0]
y = ldata[:,1]
p4.plot(x, y,pen=(255,255,255))
p4.setLabel('left', "X orbit [mm]")
p4.setLabel('bottom', "Z [m]")
x = ldata[range(1,800,10),0]
y = ldata[range(1,800,10),1]
p4.plot(x, y, pen=None,symbolBrush=(255,0,0), symbolPen='w', symbolSize=3)
p4.showGrid(x=True, y=True)
'''


## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
