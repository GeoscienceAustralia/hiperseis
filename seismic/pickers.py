from obspy.signal.trigger import ar_pick, pk_baer
from phasepapy.phasepicker.aicdpicker import AICDPicker
from phasepapy.phasepicker.ktpicker import KTPicker
from phasepapy.phasepicker.fbpicker import FBPicker


# write custom picker classes here
class PKBaer:
    """
    our wrapper class for obspy pk_baer picker class.

    For inputs and outputs refer to obspy.signal.trigge.pk_bear class
    """

    def __init__(self, tdownmax, tupevent, thr1, thr2, preset_len, p_dur):
        self.tdownmax = tdownmax
        self.tupevent = tupevent
        self.thr1 = thr1
        self.thr2 = thr2
        self.preset_len = preset_len
        self.p_dur = p_dur

    def picks(self, tr):
        """
        Parameters
        ----------
        tr: obspy.core.trace.Trace
            obspy trace instance
        samp_int: int
            number of samples per second

        Returns
        -------
        pptime: int
            pptime sample number of parrival
        pptime sample number of parrival
        pfm: str
            pfm direction of first motion (U or D)
        """
        reltrc = tr.data
        samp_int = tr.stats.sampling_rate
        pptime, pfm = pk_baer(reltrc=tr,
                              samp_int=samp_int,
                              tdownmax=self.tdownmax,
                              tupevent=self.tupevent,
                              thr1=self.thr1,
                              thr2=self.thr2,
                              preset_len=self.preset_len,
                              p_dur=self.p_dur)

        return pptime / samp_int, pfm

pickermaps = {
    'aicdpicker': AICDPicker,
    'fbpicker': FBPicker,
    'ktpicker': KTPicker,
    'pkbaer': PKBaer,
}

