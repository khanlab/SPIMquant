from spimquant.cvpl_tools.fs import ImFileType, ImReadSetting, ImIO


def reader_from_cmd(cmd: list, true_im_dim, uint8_only=False):
    fm = ImFileType.FORMAT_UINT8
    stack_axis = None
    concat_instead = False
    path_pattern = ''
    trim_front = 0
    trim_end = 0
    trim_common = False
    for i in range(0, len(cmd)):
        term = cmd[i]
        if i == 0:
            path_pattern = term
        elif term == '--probs':
            if uint8_only:
                raise ValueError('ERROR: Attempt to specify "--probs" option (implicit float32) when loading a uint8 '
                                 'image!')
            fm = ImFileType.FORMAT_FLOAT32
        elif term == '--stack':
            stack_axis = int(cmd[i + 1])
        elif term == '--concat':
            stack_axis = int(cmd[i + 1])
            concat_instead = True
        elif term == '--ord':
            if uint8_only:
                raise ValueError('ERROR: Attempt to specify "--ord" option (implicit uint16) when loading a uint8 '
                                 'image!')
            fm = ImFileType.FORMAT_UINT16
        elif term == '--trim_front':
            trim_front = int(cmd[i + 1])
        elif term == '--trim_end':
            trim_end = int(cmd[i + 1])
        elif term == '--trim_common':
            trim_common = True
    im_read_setting = ImReadSetting(
        im_format=fm,
        true_im_dim=true_im_dim,
        stack_axis=stack_axis,
        concat_instead=concat_instead,
        allow_pickle=True)
    imid_to_paths = ImIO.read_filenames(im_read_setting, path_pattern)
    if trim_front != 0 or trim_end != 0:
        imid_to_paths.trim_key_prefix_and_suffix(trim_front, trim_end)
    if trim_common:
        imid_to_paths.trim_key_common_prefix_and_suffix()

    return path_pattern, im_read_setting, imid_to_paths