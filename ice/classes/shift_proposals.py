import numpy as np
from scipy.signal import find_peaks
import ruptures as rpt
import logging


class ShiftProposals:
    filter_len=200
    def __init__(self,control_peaks,edited_peaks,mapping):
        self.base_order='ATCG'
        self.mapping=mapping
        self.edit_array=self._remap(edited_peaks,list(mapping.values()),self.base_order)
        self.control_array=self._remap(control_peaks,list(mapping.keys()),self.base_order)
        self.offset=self._compute_offset(list(mapping.values()),list(mapping.keys()))

    @staticmethod
    def _compute_offset(ctrl,edit):
        ctrl_idx = np.array(list(ctrl)).astype(np.float)
        exp_idx= np.array(list(edit)).astype(np.float)
        return np.nanmedian(ctrl_idx-exp_idx)

    @staticmethod
    def _remap(peaks,mapping,base_order):

        peak_array = np.array([peaks[base] for base in base_order]).astype(np.float)

        peak_array = np.expand_dims(peak_array, axis=0)

        # None values will be mapped to 0
        peak_array[:, :, -1] = [0, 0, 0, 0]

        indxs = np.array(list(mapping)).astype(np.float)
        indxs[np.isnan(indxs)] = -1
        indxs = indxs.astype(int)

        # get peaks from mapping
        aligned_array = np.take(peak_array, indxs, axis=2)

        #normalize peaks
        aligned_array = (aligned_array / np.sum(aligned_array, axis=1))
        aligned_array[np.isnan(aligned_array)]=0


        return aligned_array

    def _compute_discordance(self):
        discord = np.sum(abs(self.control_array - self.edit_array), axis=1).ravel()
        self.discordance=discord

        algo = rpt.Pelt(model="rbf").fit(discord.ravel())
        changepoints = algo.predict(pen=10)

        changepoints.pop(len(self.edit_array))
        if len(changepoints)==0:
            logging.warning('no discordance change points detected')
        else:
            self.changpoint=min(changepoints)
        return self.changpoint


    def _create_filters(self):

        # create filterstack from the right hand side all the way to the discordance point
        #TODO filter is protected namespace// plz fix
        filter_stack={}
        binarized_control=np.expand_dims((self.control_array == np.max(self.control_array, axis=1)).all(axis=0).astype(int),0)
        len_control=binarized_control.shape[-1]
        for i in np.arange(len_control-self.changpoint):

            temp_filter=np.zeros((1,4,self.filter_len))
            bc=binarized_control[:, :, -(self.filter_len + i):len_control - i]
            temp_filter[:,:,-bc.shape[-1]:]=bc
            filter_stack[i]=temp_filter

        self.filter_stack=filter_stack


    def _compute_single_convolution(self):
        convolved = np.zeros((self.edit_array.shape[-1]))
        filter=self.filter_stack[0]
        for idx in range(self.edit_array.shape[-1]):
            edit_slice = self.edit_array[:, :, idx:self.filter_len + idx]
            if edit_slice.shape[-1] < self.filter_len:
                padded = np.zeros((4, self.filter_len))
                padded[:, -edit_slice.shape[-1]:] = edit_slice

                convolved[idx] = np.sum(filter * padded)
            else:
                convolved[idx] = np.sum(filter * edit_slice)
        self.single_convolution = convolved

    def _compute_convolution_stack(self):

        convolution_stack = {}
        for key, filter in self.filter_stack.items():

            convolved = np.zeros((self.edit_array.shape[-1]))

            for idx in range(self.edit_array.shape[-1]):
                edit_slice = self.edit_array[:, :, idx:self.filter_len + idx]
                if edit_slice.shape[-1] < self.filter_len:
                    padded = np.zeros((4, self.filter_len))
                    padded[:, -edit_slice.shape[-1]:] = edit_slice

                    convolved[idx] = np.sum(filter * padded)
                else:
                    convolved[idx] = np.sum(filter * edit_slice)
            convolution_stack[key] = convolved


        convolution_stack_array = []
        for key in convolution_stack.keys():
            convolution_stack_array.append(np.roll(convolution_stack[key], key))
        self.convolution_stack = np.asarray(convolution_stack_array)



    def find_deletions(self):
        self._compute_discordance()
        self._create_filters()
        self._compute_single_convolution()

        peak_locs=find_peaks(self.single_convolution,prominence=25)[0]
        edit_len=self.edit_array.shape[-1]

        deletion_sizes=[]
        for peak_loc in peak_locs:
            deletion_sizes.append(edit_len - peak_loc - self.filter_len)

        self.deletion_sizes=deletion_sizes
        return deletion_sizes

    def compute_change_ranges(self):
        starts=self._compute_discordance()
        self._create_filters()
        self._compute_convolution()
        deletions=self.find_deletions()
        combinations = []

        for start in starts:
            for deletion in deletions:
                combinations.append([start-self.offset, start-self.offset + deletion])

        return combinations

