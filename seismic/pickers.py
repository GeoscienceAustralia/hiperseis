from phasepapy.phasepicker.aicdpicker import AICDPicker
from phasepapy.phasepicker.ktpicker import KTPicker
from phasepapy.phasepicker.fbpicker import FBPicker


# write custom picker classes here
all_pickers = {
    'aicdpicker': AICDPicker,
    'fbpicker': FBPicker,
    'ktpicker': KTPicker,
}
