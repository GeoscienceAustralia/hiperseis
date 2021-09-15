#!/usr/bin/env python
import os, sys
import cv2
import numpy as np
import time
from scipy.interpolate import interp1d
from PIL.PngImagePlugin import PngImageFile, PngInfo
import json
from pyproj import Geod
import click

class Translator:
    def __init__(self, state):
        self._state = state
        self._meta = None
        self._geod = Geod(ellps="WGS84")
        self._az = None
        self._baz = None

        # load metadata
        try:
            img = PngImageFile(self._state._image_fn)
            self._meta = json.loads(img.text['profile_meta'])
            img.close()
        except Exception as e:
            print(str(e))
            assert 0, 'Failed to load meta-data from {}. Aborting..'.format(self._state._image_fn)
        # end try

        # Compute azimuth and backazimuth
        self._az, self._baz, _ = self._geod.inv(self._meta['lon1'],
                                                self._meta['lat1'],
                                                self._meta['lon2'],
                                                self._meta['lat2'])

    # end func

    def output(self, file=sys.stdout):
        d = '{{"px1":{}, "py1":{},' \
            '"px2":{}, "py2":{},' \
            '"px3":{}, "py3":{},' \
            '"px4":{}, "py4":{},' \
            '"x1":{}, "x2":{},' \
            '"y1":{}, "y2":{}}}'.format(self._state._calibration_coords[0,0],
                                        self._state._calibration_coords[0,1],
                                        self._state._calibration_coords[1,0],
                                        self._state._calibration_coords[1,1],
                                        self._state._calibration_coords[2,0],
                                        self._state._calibration_coords[2,1],
                                        self._state._calibration_coords[3,0],
                                        self._state._calibration_coords[3,1],
                                        self._state.x1, self._state.x2,
                                        self._state.y1, self._state.y2)
        ds = '#' + json.dumps(json.loads(d))
        print(ds, file=file)
        print('#pixel-x pixel-y longitude    latitude    depth', file=file)
        for dc in self._state._digitization_coords:
            px = dc[0]
            py = dc[1]
            dist = dc[2]
            depth = dc[3]
            elon = None
            elat = None

            if(dist > 0):
                elon, elat, _ = self._geod.fwd(self._meta['lon1'], self._meta['lat1'],
                                               self._az, dist * 1e3) # dist in m
            else:
                elon, elat, _ = self._geod.fwd(self._meta['lon1'], self._meta['lat1'],
                                               self._baz, -dist * 1e3) # dist in m
            # end if

            print('{:7d} {:7d} {:3.7f} {:3.7f} {:5.7f}'.format(px, py, elon, elat, depth), file=file)
        # end for
    # end func
# end class

class State:
    def __init__(self, image_fn, profile_fn=None):
        self._states = {'c':'calibrate', 'd':'digitize'}
        self._state = self._states['c']
        self._image_fn = image_fn
        self._profile_fn = profile_fn
        self._master_img = cv2.imread(self._image_fn, 1)
        self._work_img = self._master_img.copy()
        self._window_name = 'image'
        self._calibration_prompts = ['Axes calibration: click on X1',
                                     'Axes calibration: click on X2',
                                     'Axes calibration: click on Y1',
                                     'Axes calibration: click on Y2',
                                     'Enter axes values for X1 X2 Y1 Y2 at the terminal']
        self._calibrated = False
        self._calibration_coords = [] # x1, x2, y1, y2
        self.x1 = 0
        self.x2 = 0
        self.y1 = 0
        self.y2 = 0

        self._xio = None
        self._yio = None
        self._digitization_coords = []
        self._translator = Translator(self)

        cv2.namedWindow(self._window_name, cv2.WINDOW_NORMAL)
        cv2.setWindowProperty(self._window_name, cv2.WND_PROP_FULLSCREEN, cv2.WINDOW_FULLSCREEN)
        cv2.setMouseCallback(self._window_name, State.mouse_callback, self)

        self._gui_message(self._calibration_prompt())
    # end func

    def _calibration_prompt(self):
        if (len(self._calibration_coords) < len(self._calibration_prompts)):
            return self._calibration_prompts[len(self._calibration_coords)]
        else:
            return None
        # end if
    # end func

    def _add_calibration_coord(self, x, y):
        if(len(self._calibration_coords) < 4):
            self._calibration_coords.append((x, y))
        # end if
    # end func

    def _complete_calibration(self):
        if(not self._calibrated and len(self._calibration_coords) == 4):
            self._calibration_coords = np.array(self._calibration_coords)
            time.sleep(0.25) # allow time for gui-thread to catch up

            while(1):
                vals = input('Enter space-separated values X1 X2 Y1 Y2 : ')
                try:
                    vals = [val.strip() for val in vals.strip().split(' ')]

                    if(len(vals) != 4): raise ValueError('Invalid number of values..')

                    self.x1 = float(vals[0])
                    self.x2 = float(vals[1])
                    self.y1 = float(vals[2])
                    self.y2 = float(vals[3])

                    break
                except Exception as e:
                    print(str(e) + ' : Re-enter values..')
                # end try
            # wend
            self._calibrated = True
            self.clear()
            self._state = self._states['d']

            # initialize interpolators
            self._xio = interp1d([self._calibration_coords[0, 0],
                                  self._calibration_coords[1, 0]],
                                 [self.x1, self.x2],
                                 fill_value='extrapolate')
            self._yio = interp1d([self._calibration_coords[2, 1],
                                  self._calibration_coords[3, 1]],
                                 [self.y1, self.y2],
                                 fill_value='extrapolate')
        # end if

        return self._calibrated
    # end func

    @staticmethod
    def mouse_callback(event, x, y, flags, self):
        if(self._state == self._states['c']):
            if (event == (cv2.EVENT_LBUTTONDOWN) and flags != (cv2.EVENT_LBUTTONDOWN + cv2.EVENT_FLAG_CTRLKEY)):
                self._add_calibration_coord(x, y)

                color = None
                if(len(self._calibration_coords)<=2):color = (200, 50, 0)
                else:color = (50, 200, 0)
                self._gui_circle(x, y, color=color, radius=7)
                self._gui_message(self._calibration_prompt())
            # end if
        elif(self._state == self._states['d']):
            if (event == (cv2.EVENT_LBUTTONDOWN) and flags == (cv2.EVENT_LBUTTONDOWN + cv2.EVENT_FLAG_CTRLKEY)):
                self._gui_circle(x, y, add_to_master=False)
                self._digitization_coords.append([x, y,
                                                  float(self._xio(x)),
                                                  float(self._yio(y))])
            # end if
        # end if
    # end func

    def clear(self):
        self._work_img = self._master_img.copy()

        if(not self._calibrated):
            self._gui_message(self._calibration_prompt())
        # end if
        cv2.imshow(self._window_name, self._work_img)
    # end func

    def _gui_circle(self, x, y, color=(0,0,0), radius=5, add_to_master=True):
        if(add_to_master):
            cv2.circle(self._master_img, (x, y), radius, color, -1)
        # end if

        cv2.circle(self._work_img, (x, y), radius, color, -1)
        cv2.imshow(self._window_name, self._work_img)
    # end func

    def _gui_message(self, text, x=None, y=None, color=(0, 0, 0), scale=5, clean=True):
        if(clean): self._work_img = self._master_img.copy()

        if(x is None and y is None):
            x, y = self._work_img.shape[1]//5, self._work_img.shape[0]//2
        # end if

        font = cv2.FONT_HERSHEY_SIMPLEX
        print(self._work_img.shape, x,y,text)
        cv2.putText(self._work_img, text, (x, y), font,
                    scale, color, 10)
        cv2.imshow(self._window_name, self._work_img)
    # end func

    def draw_digitization(self):
        self._work_img = self._master_img.copy()
        for c in self._digitization_coords:
            self._gui_circle(c[0], c[1], add_to_master=False)
        # end for
        cv2.imshow(self._window_name, self._work_img)
    # end func

    def show(self):
        cv2.imshow(self._window_name, self._work_img)
        cv2.resizeWindow(self._window_name, 1200, 800)
    # end func

    def undo_digitization_step(self):
        self._digitization_coords = self._digitization_coords[:-1]
        self.draw_digitization()
    # end func
# end class

def app():
    #fn = 'OA.BK31.-OA.BW31..png'
    fn = 'OA.BK29.-OA.BW29..png'
    state = State(fn)
    state.show()

    while(1):
        k = cv2.waitKey(7)

        if(k == ord('r')): # reset
            cv2.destroyAllWindows()
            state = State(fn)
            state.show()
        elif(k == ord('q')): # quit
            cv2.destroyAllWindows()
            break
        elif(k == ord('c')): # clear digitization
            state._digitization_coords = []
            state.clear()
        elif(k == ord('u')): # undo digitization step
            state.undo_digitization_step()
        elif(k == ord('p')): # output digitization to screen
            if(state._calibrated): state._translator.output()
            else: print('Nothing to print..')
        elif(k == ord('s')): # output digitization to file
            if(state._calibrated): state._translator.output(file=open('a.txt', 'w'))
            else: print('Nothing to save..')
        # end if

        if(not state._calibrated):
            state._complete_calibration()
        # end if
    # wend
    cv2.destroyAllWindows()
# end func

if __name__=="__main__":
    app()
# end if
