#%%
import logging
from os import initgroups
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
# import spectrum_utils.plot as sup
# import spectrum_utils.spectrum as sus
import matplotlib.ticker as mtick
import copy
from adjustText import adjust_text

ANNOTATION_COLORS = {'a': '#388E3C', 'b': '#1976D2', 'c': '#00796B',
          'x': '#7B1FA2', 'y': '#D32F2F', 'z': '#F57C00',
          '?': '#212121', 'f': '#212121', None: '#212121'}

class SpectraPlot:
    """
    Mass Spec mirror plots with optional manually annotated spectra
    """
    def __init__(self, mass_charge_ratios:list=[], intensities:list=[]) -> None:
        """[summary]

        Args:
            mass_charge_ratios (list, optional): [description]. Defaults to [].
            intensities (list, optional): [description]. Defaults to [].
        """
        self._mass_charge_ratios = np.array([])
        self._intensities = np.array([])
        self._labelled_spectra = {}
        self.intensities = intensities
        self.mass_charge_ratios = mass_charge_ratios
        if len(self.intensities) != len(self.mass_charge_ratios):
            raise ValueError("size of intensities and mass_charge_ratios must be the same")
    
    @property
    def intensities(self):
        """Raw intensities

        Returns:
            list: normalized intensity values
        """
        return self._intensities
    @intensities.setter
    def intensities(self, intensities):
        """Sets raw intensity values
        Args:
            intensities ([type]): [description]
        """
        logging.info("setting intensities...")
        self._intensities = np.array(intensities)

    @property
    def normalized_intensities(self):
        """Intensitiees normalized to largest value

        Returns:
            list : normalized intensities
        """
        return self._intensities  / self._intensities.max() * 100
    @property 
    def mass_charge_ratios(self):
        return self._mass_charge_ratios
    @mass_charge_ratios.setter
    def mass_charge_ratios(self, mass_charge_ratios):
        logging.info("setting mass charge ratios...")
        self._mass_charge_ratios = mass_charge_ratios
    
    @property
    def labelled_spectra(self):
        return self._labelled_spectra
        

    @labelled_spectra.setter
    def labelled_spectra(self, labelled_spectra):
        
        is_valid_label =  type(labelled_spectra) == dict
        for label, loc in labelled_spectra.items():
            is_valid_label = is_valid_label and  type(loc) == dict
            is_valid_label = is_valid_label and   'm/z' in loc.keys()
            is_valid_label = is_valid_label and   'intensity' in loc.keys()
        if is_valid_label:
            self._labelled_spectra = labelled_spectra
            
        else:
            raise ValueError('labelled spectra must have format {label:{"m/z": m/z, "intensity":intensity}}')

    @property
    def normalized_labelled_spectra(self):
        tmp_ls = copy.deepcopy(self._labelled_spectra)
        for _, item in tmp_ls.items():
            item['intensity'] = item['intensity'] / self._intensities.max() * 100
        return tmp_ls

    def filter(self, percent_cutoff = 5):
        intensities = []
        mass_charge_ratios = []
        for mz, i in zip(self.mass_charge_ratios, self.intensities):
            if i / self.intensities.max() * 100 > percent_cutoff:
                intensities.append(i)
                mass_charge_ratios.append(mz)
        self.intensities = intensities
        self.mass_charge_ratios = mass_charge_ratios
        
    def _mirror_plot(self, ax=None, intensities=None, mass_charge_ratios=None, grid = True, colors = None):
        if ax is None:
            fig, ax = plt.subplots()
        if colors is None:
            colors = ANNOTATION_COLORS[None]
        if intensities is None:
            intensities = self.normalized_intensities
        if mass_charge_ratios is None:
            mass_charge_ratios = self.mass_charge_ratios
        if len(mass_charge_ratios) != len(intensities):
            raise ValueError("size of mass_charge_ratios and intensities must match")
        if type(colors) is list and len(colors) != len(mass_charge_ratios):
            raise ValueError("size of mass_charge_ratios and colors must match")
        vlines = []
        if type(colors) is list:
            for mz, iz, col in  zip(mass_charge_ratios, intensities, colors):
                vlines.append(ax.vlines(x=mz, ymin=0, ymax=iz, color = col))
        else:
            for mz, iz in  zip(mass_charge_ratios, intensities):
                vlines.append(ax.vlines(x=mz, ymin=0, ymax=iz, color = colors))
        if grid:
            ax.grid(which='minor',  alpha =0.2)
            ax.grid(which='major', alpha=0.5)
            ax.minorticks_on()
        ax.set_xlabel('m/z')
        ax.set_ylabel('Intensity')
        ax.set_ylim(0,110)
        ax.yaxis.set_major_formatter(mtick.PercentFormatter())
        return ax, vlines

    def _get_label_color(self, label):
        if label is None:
            return ANNOTATION_COLORS[None]
        for annot, annot_color in ANNOTATION_COLORS.items():
            if label.startswith(annot):
                return annot_color
        return ANNOTATION_COLORS[None]
    def _labelled_mirrorplot(self, ax = None, vlines = [], adjust_annotation = False):
        if ax is None:
            fig, ax = plt.subplots()
        labels = list(self.labelled_spectra.keys())
        intensities = [self.normalized_labelled_spectra[label]['intensity'] for label in labels]
        mz = [self._labelled_spectra[label]['m/z'] for label in labels]
        colors = [self._get_label_color(label) for label in labels]

        ax, labelled_vlines = self._mirror_plot(ax, intensities, mz, colors=colors)
        vlines = vlines + labelled_vlines
        texts = []
        for m, i, l, c in zip(mz, intensities, labels, colors):
            texts.append(ax.text(x=m,y=i+1,s=l,color=c, rotation='vertical', ha='center', va='bottom'))
        # texts = ax.text(x = mz, y = intensities, s = labels, color = colors)
        if adjust_annotation:
            adjust_text(texts, add_objects=vlines,
                        autoalign=False, only_move={'points':'y', 'text':'y', 'objects':'y'},
                        ha='center', va='bottom')
        return ax

        

    def mirror_plot(self, ax=None, label=True, grid=True, adjust_annotation = False):
        if ax is None:
            fig, ax = plt.subplots()
        ax, vlines = self._mirror_plot(ax, grid=grid)
        if label:
            self._labelled_mirrorplot(ax, vlines, adjust_annotation)
        return ax

# # %%
# label = {
#     'b+1':{'m/z': 5, 'intensity': 44},
#     'c':{'m/z': 25, 'intensity': 21}}
# spectraplot = SpectraPlot(mass_charge_ratios=[4,5,6, 8, 19], intensities=[15,44,3, 34, 21])
# spectraplot.labelled_spectra = label
# ax, vlines = spectraplot.mirror_plot(label = True)
# plt.show()
