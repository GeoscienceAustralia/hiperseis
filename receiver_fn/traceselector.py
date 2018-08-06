from matplotlib.widgets import AxesWidget

class TraceSelector(AxesWidget):
    def __init__(self, ax, txt):
        AxesWidget.__init__(self, ax)

        self.t = txt
        self.connect_event('button_press_event', self._clicked)
        self._prev_color = ''
        self.cnt = 0
        self.observers = {}
        ax.set_navigate(False)

    def toggle_color():
        bbox_p = self.t.get_bbox_patch()
        if not self._prev_color:
            self.t.set_bbox(dict(facecolor='white'))
            
        self.t.set_bbox(dict(facecolor=self._prev_color))
        self._prev_color = bbox_p.get_facecolor()

    def on_clicked(self, func):
        cid = self.cnt
        self.observers[cid] = func
        self.cnt += 1
        return cid

    def disconnect(self, cid):
        try:
            del self.observers[cid]
        except KeyError:
            pass

    def _clicked(self, event):
        print('Hello there!')
        if self.ignore(event):
            return
        if event.inaxes != self.ax:
            return

        if self.t.contains(event)[0]:
            self.toggle_color()
