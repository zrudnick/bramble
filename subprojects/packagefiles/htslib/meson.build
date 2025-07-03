project('htslib', 'c')

# Get absolute paths
build_dir = meson.current_build_dir()
install_dir = join_paths(build_dir, 'install')
include_dir = join_paths(install_dir, 'include')
lib_dir = join_paths(install_dir, 'lib')

# Create the directory structure during setup
run_command('mkdir', '-p', include_dir, check: true)
run_command('mkdir', '-p', lib_dir, check: true)

# All-in-one target that does configure, build, and install
htslib_build_and_install = custom_target('htslib-build-and-install',
    command: [
        'sh', '-c', '''
        cd @CURRENT_SOURCE_DIR@/../htslib-1.22 && 
        echo "=== Configuring ===" &&
        ./configure --prefix="''' + install_dir + '''" --disable-gcs --disable-s3 &&
        echo "=== Building ===" &&
        make -j8 &&
        echo "=== Installing ===" &&
        make install &&
        echo "=== Checking installation ===" &&
        ls -la "''' + install_dir + '''/lib/" &&
        ls -la "''' + install_dir + '''/include/" &&
        echo "=== Copying library ===" &&
        cp "''' + install_dir + '''/lib/libhts.a" "''' + install_dir + '''/../libhts.a"
        '''
    ],
    output: 'libhts.a',
    build_by_default: true,
    console: true,
)

## Step 1: Configure htslib
#htslib_configure = custom_target('htslib-configure',
#    command: [
#        'sh', '-c', '''
#        cd @CURRENT_SOURCE_DIR@/../htslib-1.22 && 
#        echo "Configuring in: $(pwd)" &&
#        echo "Install prefix will be: ''' + install_dir + '''" &&
#        ./configure --prefix="''' + install_dir + '''" --disable-shared --enable-static --disable-gcs --disable-s3
#        '''
#    ],
#    output: 'configure.done',
#    build_by_default: true,
#    console: true,
#)
#
## Step 2: Build htslib
#htslib_make = custom_target('htslib-make',
#    command: [
#        'sh', '-c', '''
#        cd @CURRENT_SOURCE_DIR@/../htslib-1.22 && 
#        echo "Building in: $(pwd)" &&
#        make -j8 && 
#        touch @OUTPUT@
#        '''
#    ],
#    depends: htslib_configure,
#    output: 'make.done',
#    build_by_default: true,
#    console: true,
#)
#
## Step 3: Install htslib
#htslib_install = custom_target('htslib-install',
#    command: [
#        'sh', '-c', '''
#        cd @CURRENT_SOURCE_DIR@/../htslib-1.22 && 
#        echo "Installing from: $(pwd)" &&
#        echo "Install directory should be: ''' + install_dir + '''" &&
#        make install &&
#        echo "Contents of install directory:" &&
#        ls -la "''' + install_dir + '''" || echo "Install directory not found!" &&
#        ls -la "''' + install_dir + '''/lib" || echo "Lib directory not found!" &&
#        ls -la "''' + install_dir + '''/include" || echo "Include directory not found!" &&
#        touch @OUTPUT@
#        '''
#    ],
#    depends: htslib_make,
#    output: 'install.done',
#    build_by_default: true,
#    console: true,
#)
#
## Step 4: Create the actual library artifact that Meson can link against
#htslib_lib = custom_target('htslib-lib',
#    command: [
#        'sh', '-c', '''
#        echo "Looking for library at: ''' + install_dir + '''/lib/libhts.a" &&
#        ls -la "''' + install_dir + '''/lib/" &&
#        cp "''' + install_dir + '''/lib/libhts.a" @OUTPUT@
#        '''
#    ],
#    depends: htslib_install,
#    output: 'libhts.a',
#    build_by_default: true,
#    console: true,
#)

# Create the dependency
htslib_inc = include_directories(join_paths('install', 'include'))
htslib_dep = declare_dependency(
    link_with: htslib_build_and_install,
    include_directories: include_directories(join_paths('install', 'include')),
)

# This is what makes it available to `get_variable()`
meson.override_dependency('htslib', htslib_dep)
htslib_lib = htslib_build_and_install
htslib_inc = include_directories(join_paths('install', 'include'))

