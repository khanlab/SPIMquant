"""
This file creates an overlay comparison between the image and its mask using Napari.
A set of pairs of image/mask files are assumed to exist when this program is executed, and their paths
will be passed in as snakemake parameters to this program. Then this program will
run Napari to display the overlay of each pair of image files
"""


import napari
import numpy as np
from magicgui import magicgui
import os


PATH_PREFIX = snakemake.input.pred_result_dir
BLOCK_WIDTH = 4096  # take as much as possible
iim = -1
while True:
    iim += 1
    if not os.path.exists(f'{PATH_PREFIX}/inference_im_{iim}.npy'):
        print('All inferences have been viewed. Exiting...')
        break

    im = np.load(f'{PATH_PREFIX}/inference_im_{iim}.npy')[:BLOCK_WIDTH, :BLOCK_WIDTH, :BLOCK_WIDTH]
    seg = np.load(f'{PATH_PREFIX}/inference_{iim}.npy')[:, :BLOCK_WIDTH, :BLOCK_WIDTH, :BLOCK_WIDTH]
    seg_all = seg.argmax(axis=0)
    seg_bin = np.ones_like(im, np.int64)
    viewer = napari.Viewer()
    viewer.add_image(im, name='im', colormap='bop blue')
    viewer.add_image(seg_all, name='seg', colormap='tab10')
    for i in range(10):
        seg_bin[0, 0, i] = i
    seg_bin_layer = viewer.add_image(seg_bin, name='seg_bin', colormap='tab10')

    DEFAULT_THRES = 1.01
    SLIDER_DICT = {"widget_type": "FloatSlider", 'max': 1.01, 'min': -1.01}
    NCLASSES = 10

    INIT_STR = "seg_bin = seg[0] > 0"
    @magicgui(
        expr={"widget_type": "TextEdit"},
        save_button={"widget_type": "PushButton", "text": "save thresholds"},
        load_button={"widget_type": "PushButton", "text": "load thresholds"},
    )
    def widget_seg_reweight(
        expr=INIT_STR,
        save_button=False,
        load_button=False,
    ):
        global seg_bin
        locals = {'seg_bin': seg_bin, 'seg': seg}
        exec(expr, globals(), locals)
        seg_bin = locals['seg_bin']
        seg_bin_layer.data = np.array(seg_bin, dtype=np.int64)
        seg_bin_layer.refresh()


    def get_class_slider_floats() -> np.array:
        sf = np.zeros((NCLASSES, ), dtype=np.float32)
        widget_dict = widget_seg_reweight.asdict()
        for i in range(NCLASSES):
            key = f'slider_float{i + 1}'
            sf[i] = widget_dict[key]
        return sf


    viewer.window.add_dock_widget(widget_seg_reweight)
    napari.run()
