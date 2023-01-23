#!/usr/bin/env python
import os, sys, re
import cv2
import numpy as np
from scipy.interpolate import interp1d
from PIL.PngImagePlugin import PngImageFile, PngInfo
import json
from pyproj import Geod
import click

class State:
    def __init__(self, image_fn, profile_fn=None):
        self._image_fn = image_fn
        self._profile_fn = profile_fn
        self._master_img = cv2.imread(self._image_fn, 1)
        self._work_img = self._master_img.copy()
        self._window_name = 'image'
        self._calibrated = False
        self._px1 = None
        self._py1 = None
        self._px2 = None
        self._py2 = None

        self._xio = None
        self._yio = None
        self._digitization_coords = []
        self._geod = None
        self._az = None
        self._baz = None

        cv2.namedWindow(self._window_name, cv2.WINDOW_NORMAL)
        cv2.setMouseCallback(self._window_name, State.mouse_callback, self)
    # end func

    def _initialize_translator(self):
        self._geod = Geod(ellps="WGS84")

        # Compute azimuth and backazimuth
        self._az, self._baz, _ = self._geod.inv(self._meta['lon1'],
                                                self._meta['lat1'],
                                                self._meta['lon2'],
                                                self._meta['lat2'])
    # end func

    def _calibrate(self):
        print('Calibrating..')
        if(not self._calibrated):
            # load metadata
            try:
                img = PngImageFile(self._image_fn)
                self._meta = json.loads(img.text['profile_meta'])
                img.close()
            except Exception as e:
                print(str(e))
                assert 0, 'Failed to load meta-data from {}. Aborting..'.format(self._image_fn)
            # end try

            # find min/max pixel coordinates
            gray = cv2.cvtColor(self._master_img.copy(), cv2.COLOR_BGR2GRAY)
            blur = cv2.GaussianBlur(gray, (5, 5), 0)
            edges = cv2.Canny(blur, threshold1=200, threshold2=240, L2gradient=True)
            contours, _ = cv2.findContours(edges, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
            sorted_contours = sorted(contours, key=cv2.contourArea, reverse=True)
            approx = cv2.approxPolyDP(sorted_contours[0],
                                      0.01*cv2.arcLength(sorted_contours[0], True), True)
            if(len(approx)!=4): assert 0, 'Calibration failed. Aborting..'

            coords = np.squeeze(np.array(approx))
            self._px1 = np.max(np.sort(coords[:, 0])[:2])
            self._py1 = np.max(np.sort(coords[:, 1])[:2])
            self._px2 = np.min(np.sort(coords[:, 0])[-2:])
            self._py2 = np.min(np.sort(coords[:, 1])[-2:])

            print('Calibration success: pixel_min({},{})->coord_min({:.2f},{:.2f}), pixel_max({},{})->coord_max({:.2f},{:.2f})'.format(
                   self._px1, self._py1, self._meta['x1'], self._meta['y1'],
                   self._px2, self._py2, self._meta['x2'], self._meta['y2']))

            # initialize interpolators
            self._xio = interp1d([self._px1, self._px2],
                                 [self._meta['x1'], self._meta['x2']],
                                 fill_value='extrapolate')
            self._yio = interp1d([self._py1, self._py2],
                                 [self._meta['y1'], self._meta['y2']],
                                 fill_value='extrapolate')

            # draw calibration markers
            self._gui_cross_hairs(self._px1, self._py1, color=(255,0,255), radius=25,
                                  thickness=2, add_to_master=True)
            self._gui_cross_hairs(self._px2, self._py2, color=(255,0,255), radius=25,
                                  thickness=2, add_to_master=True)

            self._calibrated = True
            if(self._profile_fn): self._load_profile()
            self._initialize_translator()
        # end if
    # end func

    def _load_profile(self):
        print('Loading profile {} ..'.format(self._profile_fn))
        lines = open(self._profile_fn, 'r').readlines()
        assert lines[0][0] == '#', 'Invalid profile file. Aborting..'

        pmeta = json.loads(lines[0][1:])
        # sanity check: ensure calibration matches
        msg = 'Calibration mismatch detected between image {} and profile {}. Aborting..'.format(self._image_fn, self._profile_fn)
        assert self._px1 == pmeta['px1'], msg
        assert self._py1 == pmeta['py1'], msg
        assert self._px2 == pmeta['px2'], msg
        assert self._py2 == pmeta['py2'], msg
        assert self._meta['x1'] == pmeta['x1'], msg
        assert self._meta['y1'] == pmeta['y1'], msg
        assert self._meta['x2'] == pmeta['x2'], msg
        assert self._meta['y2'] == pmeta['y2'], msg

        self._digitization_coords = []
        for i in np.arange(2, len(lines)): # two line header
            line = lines[i]
            print(line)
            px, py, lon, lat, depth = map(float, list(filter(len,re.split('\s+', line))))
            px = int(px)
            py = int(py)
            self._digitization_coords.append([px, py, lon, lat, depth])
        # end for

        self.draw_digitization()
    # end func

    @staticmethod
    def mouse_callback(event, x, y, flags, self):
        if(not self._calibrated): return
        if (event == (cv2.EVENT_LBUTTONDOWN) and flags == (cv2.EVENT_LBUTTONDOWN + cv2.EVENT_FLAG_CTRLKEY)):
            self._gui_circle(x, y, add_to_master=False)
            self._digitization_coords.append([x, y,
                                              float(self._xio(x)),
                                              float(self._yio(y))])
        # end if
    # end func

    def clear(self):
        self._work_img = self._master_img.copy()
        cv2.imshow(self._window_name, self._work_img)
    # end func

    def _gui_circle(self, x, y, color=(0,0,0), radius=5, add_to_master=True):
        if(add_to_master):
            cv2.circle(self._master_img, (x, y), radius, color, -1)
        # end if

        cv2.circle(self._work_img, (x, y), radius, color, -1)
        cv2.imshow(self._window_name, self._work_img)
    # end func

    def _gui_cross_hairs(self, x, y, color=(0,0,0), radius=10, thickness=1, add_to_master=True):
        if(add_to_master):
            cv2.line(self._master_img, (x-radius-thickness, y), (x+radius+thickness, y), color, thickness=thickness)
            cv2.line(self._master_img, (x, y-radius-thickness), (x, y+radius+thickness), color, thickness=thickness)
            cv2.circle(self._master_img, (x, y), radius, color, thickness=thickness)
        # end if

        cv2.line(self._work_img, (x-radius-thickness, y), (x+radius+thickness, y), color, thickness=thickness)
        cv2.line(self._work_img, (x, y-radius-thickness), (x, y+radius+thickness), color, thickness=thickness)
        cv2.circle(self._work_img, (x, y), radius, color, thickness=thickness)
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
        cv2.resizeWindow(self._window_name, 1600, 900)
    # end func

    def undo_digitization_step(self):
        self._digitization_coords = self._digitization_coords[:-1]
        self.draw_digitization()
    # end func

    def output(self, file=sys.stdout):
        d = '{{"px1":{}, "py1":{},' \
            '"px2":{}, "py2":{},' \
            '"x1":{}, "y1":{},' \
            '"x2":{}, "y2":{}}}'.format(self._px1,
                                        self._py1,
                                        self._px2,
                                        self._py2,
                                        self._meta['x1'],
                                        self._meta['y1'],
                                        self._meta['x2'],
                                        self._meta['y2'])
        ds = '#' + json.dumps(json.loads(d))
        print(ds, file=file)
        print('#pixel-x pixel-y longitude    latitude    distance    depth', file=file)
        for dc in self._digitization_coords:
            px = dc[0]
            py = dc[1]
            dist = dc[2]
            depth = dc[3]
            elon = None
            elat = None

            if(dist > 0):
                elon, elat, _ = self._geod.fwd(self._meta['lon1'],
                                               self._meta['lat1'],
                                               self._az, dist * 1e3) # dist in m
            else:
                elon, elat, _ = self._geod.fwd(self._meta['lon1'], self._meta['lat1'],
                                               self._baz, -dist * 1e3) # dist in m
            # end if

            print('{:7d} {:7d} {:3.7f} {:3.7f} {:5.7f} {:5.7f}'.format(px, py, elon, elat, dist, depth), file=file)
        # end for
    # end func
# end class

gui_help = \
    """
    =========================================================================
        DIGITIZER GUI HELP
            Mouse: 
                CTRL+LEFT-click : add digitization point
                
            Key-stroke:
                h : prints this help message
                p : prints current digitization points to terminal
                s : saves current digitization points to output file
                c : clears all digitization points
                u : undo previous digitization step
                q : quits session, after saving digitization to output file
                
    =========================================================================
    """
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('image-fn', type=click.Path(exists=True, dir_okay=False))
@click.option('--output-file-name', default=None, type=click.Path(dir_okay=False),
              help='Output file-name; default is <image-fn>.txt')
@click.option('--load-profile', default=None, type=click.Path(exists=True, dir_okay=False),
              help='Load existing digitization for image')
def app(image_fn, output_file_name, load_profile):
    """
        --------CCP Vertical-Profile Digitizer--------\n
        IMAGE_FN: Image file-name (output of plot_ccp.py)

        Input images are automatically calibrated. Upon successful calibration,
        markers are placed at the origin and at (maxX, maxY).

        Press 'h' on GUI for further help
    """

    print(gui_help)
    state = State(image_fn, profile_fn=load_profile)
    state.show()

    ofn = None
    if(output_file_name): ofn = output_file_name
    else: ofn = os.path.join(os.path.dirname(image_fn),
                             os.path.splitext(os.path.basename(image_fn))[0] + '.txt')

    while(1):
        k = cv2.waitKey(7)

        if(k == ord('h')):
            print(gui_help)
        elif(k == ord('c')): # clear digitization
            print('Clearing all digitization points..')
            state._digitization_coords = []
            state.clear()
        elif(k == ord('u')): # undo digitization step
            state.undo_digitization_step()
        elif(k == ord('p')): # output digitization to screen
            state.output()
        elif(k == ord('s')): # output digitization to file
            print('Saving output to {} ..'.format(ofn))
            state.output(file=open(ofn, 'w'))
        elif(k == ord('q')): # quit
            print('Saving output to {} ..'.format(ofn))
            state.output(file=open(ofn, 'w'))
            cv2.destroyAllWindows()
            break
        # end if

        if(not state._calibrated): state._calibrate()
    # wend

    cv2.destroyAllWindows()
# end func

if __name__=="__main__":
    app()
# end if
