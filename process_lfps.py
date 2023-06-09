##

import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt, find_peaks, hilbert, windows

def butter_bandpass(lowcut, highcut, fs, order=3):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=3):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)

    filtered = filtfilt(b, a, data) #filtfilt

    analytic_signal = hilbert(filtered)
    phase = np.arctan2(filtered,analytic_signal.imag)

    amp = np.abs(analytic_signal)

    return filtered, phase, amp, analytic_signal


# multitaper spectral estimates

def create_data_segments(data,seg_len,seg_step,df_col_names):
    """data: array, of values to be broken up into distinct segments, typically one trial
       seg_len: scalar, number of samples per segment, e.g. 500
       seg_step: scalar, number of samples to advance the window e.g. 250
       df_col_names: list, of string column name for data, and for segment label col
       TODO: figure out how to adaptively find the seg_len and seg_overlap that
       yield integer number of main and overlap windows
        """

    full_len = len(data)
    # seg_start = np.arange(0, full_len - seg_len, seg_step)
    seg_start = np.arange(0, full_len, seg_step)
    seg_end = np.arange(seg_len, full_len+1, seg_step) #add one bc non-inclusive otherwise

    seg_df = []
    for i in range(len(seg_end)):

        #remember indexing is NOT right-inclusive
        d = {df_col_names[0]: data[seg_start[i]:seg_end[i]],
             df_col_names[1]: np.repeat(i,seg_len)}

        tmp_df = pd.DataFrame(d)
        seg_df.append(tmp_df)

    segmented_data_df = pd.concat(seg_df)

    return segmented_data_df

def compute_dpss_tapers(seg_len,time_half_bw_pdt,k_windows,sym):
    """This wrapper uses scipy.signal.windows.dpss to compute tapers and saves
    them to a pandas df with the taper label
    sym: bool, True(default) generates symmetric window for use in filter design;
    False, generates periodic window for use in spectral analysis"""

    tapers = windows.dpss(seg_len,time_half_bw_pdt,k_windows,sym)

    tapers_df = pd.DataFrame([])
    df = []
    for t in range(tapers.shape[0]):

        d = {'taper': tapers[t],
             'taper_label': np.repeat(t,seg_len)}
        df.append(pd.DataFrame(d))

    tapers_df = pd.concat(df)

    return tapers_df

def compute_tapered_timeseries(segmented_data_df,data_col_name,tapers_df):

    seg_labels = list(set(segmented_data_df[seglabel_col_name].to_list()))
    taper_labels = list(set(tapers_df[taperlabel_col_name].to_list()))

    df = []
    for seg in seg_labels:

        subseg = segmented_data_df[segmented_data_df['segment_label'] == seg]

        dataseg = subseg[data_col_name].values

        for taper in taper_labels:

            subtaper = tapers_df[tapers_df['taper_label'] == taper]
            tapervals = subtaper['taper'].values

            taperdata = dataseg * tapervals

            d = {'tapered_data': taperdata,
                 'segment_label': np.repeat(seg,len(taperdata)),
                 'taper_label': np.repeat(taper,len(taperdata))
            }

            df.append(pd.DataFrame(d))

    tapereddata_df = pd.concat(df)

    return tapereddata_df

def compute_butterHilbert_phase_per_tapereddata(tapereddata_df,seglabel_col_name,taperlabel_col_name,data_col_name,freq_band,fs):
    """freq_band, list, of low and high frequency cutoffs"""

    seg_labels = list(set(tapereddata_df[seglabel_col_name].to_list()))
    taper_labels = list(set(tapereddata_df[taperlabel_col_name].to_list()))

    low = freq_band[0]
    high = freq_band[1]

    df = []
    for seg in seg_labels:

        subseg = tapereddata_df[tapereddata_df[seglabel_col_name] == seg]

        for taper in taper_labels:

            subtaper = subseg[subseg[taperlabel_col_name] == taper]
            dataseg = subtaper[data_col_name].values

            _, phase, _, _ = butter_bandpass_filter(dataseg,low,high,fs)

            d = {'phase': phase,
                 'segment_label': np.repeat(seg,len(phase)),
                 'taper_label': np.repeat(taper,len(phase))
            }
            df.append(pd.DataFrame(d))

    phase_per_seg_df = pd.concat(df)

    return phase_per_seg_df

def circular_mean(phase_angles):

    tmpcos = [np.cos(phi) for phi in phase_angles]
    tmpsin = [np.sin(phi) for phi in phase_angles]

    meancos = np.mean(tmpcos)
    meansin = np.mean(tmpsin)

    circ_mean = np.arctan2(meansin,meancos)
    return circ_mean

def label_phase_cycles(df, lowcut, highcut, fs, lfp_trial_ind):
    """This function specifically: filters a raw LFP, finds troughs in the
    resulting phase time Series, and adds the trough-based cycle labels to the
    dataframe the lfp came in label
    df: dataframe, need to make sure that if this is a subset of another DataFrame
    you pass in the subdf.reset_index() version, so that the cycle labels can be
    assigned for the corresponding indices"""

    _, phase, _, _ = butter_bandpass_filter(df['lfp'].values,lowcut,highcut,fs)

    # find the troughs
    troughs = find_peaks(-phase,
                         height=None,
                         threshold=None,
                         distance=None,
                         prominence=None,
                         width=None,
                         wlen=None,
                         rel_height=0.5,
                         plateau_size=None)

    all_inds = np.hstack((0,troughs[0],len(phase)-1))
    periods = np.diff(all_inds)
    trial_adjusted_inds = [i+(lfp_trial_ind*df.shape[0]) for i in all_inds]

    cycle_labels = np.arange(0,len(periods))

    df['cycle_labels'] = pd.cut(df.index.to_list(),
                                  bins=trial_adjusted_inds,
                                  labels=cycle_labels,
                                  include_lowest=True
                                 )

    return df

def label_filtered_cycles(df, lowcut, highcut, fs, add_filtered_to_df, add_phase_to_df):
    """This function specifically: filters a raw LFP, finds troughs in the
    resulting phase time Series, and adds the trough-based cycle labels to the
    dataframe the lfp came in label
    df: dataframe, need to make sure that if this is a subset of another DataFrame
    you pass in the subdf.reset_index() version, so that the cycle labels can be
    assigned for the corresponding indices"""

    filtered, phase, _, _ = butter_bandpass_filter(df['lfp'].values,lowcut,highcut,fs)

    # find the troughs
    troughs = find_peaks(-filtered,
                         height=None,
                         threshold=None,
                         distance=None,
                         prominence=None,
                         width=None,
                         wlen=None,
                         rel_height=0.5,
                         plateau_size=None)

    all_inds = np.hstack((0,troughs[0],len(filtered)-1))
    periods = np.diff(all_inds)
    trial_adjusted_inds = [i+df.index[0] for i in all_inds]

    cycle_labels = np.arange(0,len(periods))

    df['cycle_labels'] = pd.cut(df.index.to_list(),
                                  bins=trial_adjusted_inds,
                                  labels=cycle_labels,
                                  include_lowest=True
                                 )

    if add_filtered_to_df:

         df['filtered_' + str(lowcut) + '_' + str(highcut)] = filtered

    if add_phase_to_df:

         df['phase_'+str(lowcut)+'_'+str(highcut)] = phase

    return df

def make_rhythm_on_off_timeseries(df,low_freq,high_freq,fs):

    #define the range of durations of a given rhythm's period
    cycle_boundaries=[(1/high_freq)*fs, (1/low_freq)*fs]

    #grab the list of cycle labels from the dataframe
    cycle_labels = list(set(df['cycle_labels'].to_list()))

    rhythm_on = []
    for cycle in cycle_labels:

        subcycle = df[df['cycle_labels']==cycle]
        dur = subcycle.shape[0]

        #if the period is within the duration boundaries, label all the timesamples
        #corresponding to that cycle with a 1 (rhythm is ON during this cycle)
        if (dur >= cycle_boundaries[0]) & (dur <= cycle_boundaries[1]):

            rhythm_on.append(np.repeat(1,dur))

        else:

            rhythm_on.append(np.repeat(0,dur))

    rhythm_on = np.concatenate(rhythm_on)

    df['rhythm_on'] = rhythm_on

    return df
