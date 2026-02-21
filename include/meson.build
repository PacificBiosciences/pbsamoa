inc_dir = include_directories('.')

if meson.is_subproject()
    subdir_done()
endif

install_subdir('pbsamoa', install_dir: get_option('includedir'))
