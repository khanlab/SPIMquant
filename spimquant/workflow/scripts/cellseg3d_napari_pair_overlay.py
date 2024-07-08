import napari
import numpy as np
from magicgui import magicgui
from magicgui.experimental import guiclass, button
from magicgui import widgets
from skimage.morphology import binary_erosion, binary_dilation
from cvpl_tools.np_algs import watershed, instance_to_binary, round_object_detection
import glob
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


    # @widget_seg_reweight.save_button.changed.connect
    # def save_button_callback():
    #     with open('./model_config/thres_conf.ini', 'w') as outfile:
    #         outfile.write(widget_seg_reweight.expr.value)
    #
    #
    # @widget_seg_reweight.load_button.changed.connect
    # def load_button_callback():
    #     with open('./model_config/thres_conf.ini', 'r') as infile:
    #         expr = infile.read().strip()
    #         widget_seg_reweight.expr.value = expr


    viewer.window.add_dock_widget(widget_seg_reweight)
    napari.run()
