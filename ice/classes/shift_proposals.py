import numpy as np
import ruptures as rpt
import logging
from itertools import combinations
np.seterr(divide='ignore')
import pickle

class ShiftProposals:

    filter_len=50
    sweep_range= np.arange(-30, 31)
    sg_sweep_range= np.arange(-50, 51)

    def __init__(self,
                 control_peaks:dict,
                 edited_peaks:dict,
                 mapping:dict,
                 epc=None,
                 guide_targets=None,
                 changpoint=None):

        self.base_order='ATCG'
        self.mapping=mapping
        self.edit_array=self._remap(edited_peaks,list(mapping.values()),self.base_order)
        self.control_array=self._remap(control_peaks,list(mapping.keys()),self.base_order)
        self.epc=epc
        self.guide_targets=guide_targets
        self.changpoint=changpoint

    @staticmethod
    def dump_dictionary(mapping,control,edit,output):
        shift_dict = {}
        shift_dict['mapping']=mapping
        shift_dict['control_peaks']=control
        shift_dict['edited_peaks']=edit

        with open(output, 'wb') as handle:
            pickle.dump(shift_dict, handle)



    def convert_primary_bc(self):
        pass


    @staticmethod
    def _remap(peaks,mapping,base_order) -> np.ndarray:

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
        with np.errstate(divide='ignore', invalid='ignore'):
            aligned_array = (aligned_array / np.sum(aligned_array, axis=1))
        aligned_array[np.isnan(aligned_array)]=0


        return aligned_array

    def _compute_discordance(self) -> None:
        discord = np.sum(abs(self.control_array - self.edit_array), axis=1).ravel()
        self.discordance=discord

        algo = rpt.Pelt(model="rbf").fit(discord.ravel())
        changepoints = algo.predict(pen=10)


        if len(changepoints)==0:
            logging.warning('no discordance change points detected')
            self.changpoint=None
        else:
            if len(self.edit_array) in changepoints:
                changepoints.pop(len(self.edit_array))
            self.changpoint=min(changepoints)


    def _create_filters(self) -> None:

        # create filterstack from the right hand side all the way to the discordance point
        filter_stack={}
        binarized_control=np.expand_dims((self.control_array == np.max(self.control_array, axis=1)).all(axis=0).astype(int),0)
        len_control=binarized_control.shape[-1]
        for i in np.arange(len_control):

            temp_filter=np.zeros((1,4,self.filter_len))
            bc=binarized_control[:, :, -(self.filter_len + i):len_control - i]
            temp_filter[:,:,-bc.shape[-1]:]=bc
            filter_stack[i]=temp_filter

        self.filter_stack=filter_stack





    def _compute_single_convolution(self) -> None:
        convolved = np.zeros((self.edit_array.shape[-1]))
        filter=self.filter_stack[0]
        for idx in range(self.edit_array.shape[-1]):
            edit_slice = self.edit_array[:, :, idx:self.filter_len + idx]
            if edit_slice.shape[-1] < self.filter_len:
                padded = np.zeros((4, self.filter_len))
                padded[:, -edit_slice.shape[-1]:] = edit_slice

                convolved[idx] = np.sum(filter * padded)/self.filter_len
            else:
                convolved[idx] = np.sum(filter * edit_slice)/self.filter_len
        self.single_convolution = convolved

    def _compute_convolution_stack(self) -> None:

        convolution_stack = {}
        for key, filter in self.filter_stack.items():

            convolved = np.zeros((self.edit_array.shape[-1]))

            for idx in range(self.edit_array.shape[-1]):
                edit_slice = self.edit_array[:, :, idx:self.filter_len + idx]
                if edit_slice.shape[-1] < self.filter_len:
                    padded = np.zeros((4, self.filter_len))
                    padded[:, -edit_slice.shape[-1]:] = edit_slice

                    convolved[idx] = np.sum(filter * padded)/self.filter_len
                else:
                    convolved[idx] = np.sum(filter * edit_slice)/self.filter_len
            convolution_stack[key] = convolved


        convolution_stack_array = []
        for key in convolution_stack.keys():
            convolution_stack_array.append(np.roll(convolution_stack[key], key))

        convolution_stack_array_unrolled = []
        for key in convolution_stack.keys():
            convolution_stack_array_unrolled.append(convolution_stack[key])

        self.convolution_stack = np.asarray(convolution_stack_array)
        self.convolution_stack_array_unrolled=np.asarray(convolution_stack_array_unrolled)


    def _find_deletions_from_stack(self) -> np.ndarray:
        if self.changpoint is None:
            self._compute_discordance()

        if self.changpoint is not None:
            self._create_filters()
            self._compute_convolution_stack()
            rolled_stack = np.roll(self.convolution_stack, self.filter_len)

            variance_stack = np.max(rolled_stack, axis=0)
            edit_len = self.edit_array.shape[-1]
            return edit_len-variance_stack.argsort()[-30:][::-1]
        else:
            return []






    def get_multiguide_proposals(self) -> list:
        '''
        This method generates a list of proposals based around MG dropout locations and the
        :return:
        '''

        deletions = np.asarray(self._find_deletions_from_stack())
        deletions=deletions[deletions<200]

        print(f'deletions found : {deletions}')

        dropout_dict = {}

        for combo in combinations(self.guide_targets, 2):
            dropout_size = combo[1].cutsite - combo[0].cutsite
            dropout_dict[dropout_size] = combo

        deletion_sizes = list(dropout_dict.keys())
        mg_proposals=[]

        for deletion in deletions:
            # find the closet guide to both

            nearest_index = (np.abs(np.asarray(deletion_sizes) - deletion)).argmin()
            g1, g2 = dropout_dict[deletion_sizes[nearest_index]]

            deletion_delta = deletion - deletion_sizes[nearest_index]

            default_dels = [0, deletion_delta]



            left_offset = self.sweep_range - default_dels[0]
            right_offset = default_dels[1] - self.sweep_range

            for r in np.arange(len(self.sweep_range)):
                '''
                r is how many shifts to the left  this is

                How do we deal with the shifting indels:
                we see there are two values for each cutsite
                cut1=(additional_deletions,0)
                cut2=(0,additional_deletions)

                for the additional deletions values, those can either be positive or negative. If they're positive, that
                means you're deleting, if they're negative they're "insertion like" in the sense that you're moving
                the cutsite away


                '''
                cut1_del = (int(left_offset[r]), 0)
                cut2_del = (0, int(right_offset[r]))
                dropout = self.epc.multiplex_proposal(
                    g1.cutsite,
                    g2.cutsite,
                    g1.label,
                    g2.label,
                    cut1_del=cut1_del, cut2_del=cut2_del,
                    dropout=True
                )
                mg_proposals.append(dropout)



        ## creating dropout + change

        return mg_proposals



    def get_singleguide_proposals(self) -> list:
        '''
        This method generates a list of proposals based around MG dropout locations and the
        :return:
        '''

        deletions = np.asarray(self._find_deletions_from_stack())

        # for singleplex we only tolerate deletions greater less than 60
        deletions=deletions[deletions<60]
        print(f'deletions found : {deletions}')


        sg_proposals=[]

        for deletion in deletions:
            # find the closet guide to both

            default_dels = [0, deletion]

            left_offset = self.sg_sweep_range - default_dels[0]
            right_offset = default_dels[1] - self.sg_sweep_range

            for r in np.arange(len(self.sg_sweep_range)):
                '''
                This is 

                '''
                cut1_del = (int(left_offset[r]), int(right_offset[r]))

                shift_cut=self.epc.single_cut_edit_proposal(self.guide_targets[0].cutsite,
                                             self.guide_targets[0].label,
                                             del_before=cut1_del[0], del_after=cut1_del[1])


                sg_proposals.append(shift_cut)

        return sg_proposals



