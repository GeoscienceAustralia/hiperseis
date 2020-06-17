python ~/dev/hiperseis/seismic/plot_network_event_dataset.py --channel-order ZRT synth_events_3L.h5 synth_events_3L.pdf

python ~/dev/hiperseis/seismic/receiver_fn/generate_rf.py --deconv-domain freq synth_events_3L.h5 synth_rf_3L_ZRT_fd_rev1.h5
python ~/dev/hiperseis/seismic/receiver_fn/generate_rf.py --deconv-domain time synth_events_3L.h5 synth_rf_3L_ZRT_td_rev1.h5
python ~/dev/hiperseis/seismic/receiver_fn/generate_rf.py --deconv-domain iter synth_events_3L.h5 synth_rf_3L_ZRT_it_rev1.h5

python ~/dev/hiperseis/seismic/receiver_fn/rf_quality_filter.py synth_rf_3L_ZRT_fd_rev1.h5 synth_rf_3L_ZRT_fd_rev1_qual.h5
python ~/dev/hiperseis/seismic/receiver_fn/rf_quality_filter.py synth_rf_3L_ZRT_td_rev1.h5 synth_rf_3L_ZRT_td_rev1_qual.h5
python ~/dev/hiperseis/seismic/receiver_fn/rf_quality_filter.py synth_rf_3L_ZRT_it_rev1.h5 synth_rf_3L_ZRT_it_rev1_qual.h5

python ~/dev/hiperseis/seismic/receiver_fn/bulk_rf_report.py --hk-weights  0.34 0.33 0.33 --hk-hpf-freq 0.2 synth_rf_3L_ZRT_fd_rev1_qual.h5 synth_rf_3L_ZRT_fd_rev1_qual_nofilt.pdf
python ~/dev/hiperseis/seismic/receiver_fn/bulk_rf_report.py --hk-weights  0.34 0.33 0.33 --hk-hpf-freq 0.2 synth_rf_3L_ZRT_td_rev1_qual.h5 synth_rf_3L_ZRT_td_rev1_qual_nofilt.pdf
python ~/dev/hiperseis/seismic/receiver_fn/bulk_rf_report.py --hk-weights  0.34 0.33 0.33 --hk-hpf-freq 0.2 synth_rf_3L_ZRT_it_rev1_qual.h5 synth_rf_3L_ZRT_it_rev1_qual_nofilt.pdf

python ~/dev/hiperseis/seismic/inversion/wavefield_decomp/runners.py single-job --waveform-file synth_events_3L.h5 --network SY --station OAA --output-file synth_events_3L_soln.h5 config_3L.json

python ~/dev/hiperseis/seismic/inversion/wavefield_decomp/plot_nd_batch.py --output-file synth_events_3L_soln.pdf synth_events_3L_soln.h5

