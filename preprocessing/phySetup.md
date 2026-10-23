```
# Installing phy2
cd to the folder you have this file: phy2_environment.yml

conda env create -f phy2_environment.yml

# Running phy2
1. open your conda terminal
2. conda activate phy2
3. cd to the desired Mouse_session_imec0_shankx/sorter_output folder

If it is the first time you're running phy for that probe/shank, you may need to extract waveforms. See how to check below.
* run this next in your conda terminal: phy extract-waveforms params.py
* then do the next step

If you have previously extracted waveforms, you can go ahead and open the gui for spike sorting
* run this in your conda terminal: phy template-gui params.py


If you don't know if you already extracted waveforms, don't guess. Accidentally re-running it makes your life exponentially worse. Instead, just check inside your Mouse_session_imec0_shankx/sorter_output folder. If you see the following three files, you already ran extract-waveforms params.py:

_phy_spikes_subset.channels.npy
_phy_spikes_subset.spikes.npy
_phy_spikes_subset.waveforms.npy


```
