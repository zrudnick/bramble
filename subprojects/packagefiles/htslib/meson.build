project('htslib', 'c')

# Get absolute paths
source_dir = meson.current_source_dir()
build_dir = meson.current_build_dir()
install_dir = join_paths(build_dir, 'install')
include_dir = join_paths(install_dir, 'include')
lib_dir = join_paths(install_dir, 'lib')

# Create the directory structure during setup
run_command('mkdir', '-p', include_dir, check: true)
run_command('mkdir', '-p', lib_dir, check: true)

# Get source files to track --- uglyyy
source_files = files([
    '../htslib-1.22/configure',
    '../htslib-1.22/configure.ac',
    '../htslib-1.22/version.sh',
])

# Get all C source files (this will make rebuilds more accurate)
c_source_files = run_command('find', source_dir + '/../htslib-1.22', '-name', '*.c', '-o', '-name', '*.h', '-i', '-name', '*.in', check: false)
if c_source_files.returncode() == 0
    # Add the source files if find command succeeded
    source_files += files(c_source_files.stdout().strip().split('\n'))
endif

# Step 1: Configure htslib
htslib_configure = custom_target('htslib-configure',
    input: source_files,  # Use all source files
    command: [
        'sh', '-c', '''
        cd ''' + source_dir + '''/../htslib-1.22 &&
        echo "=== CONFIGURE DEBUG INFO ===" &&
        echo "Configuring in: $(pwd)" &&
        echo "Install prefix will be: ''' + install_dir + '''" &&
        echo "Output file will be: ''' + build_dir + '''/configure.done" &&
        echo "Checking if configure.done exists: $(ls -la ''' + build_dir + '''/configure.done 2>/dev/null || echo 'NOT FOUND')" &&
        echo "Checking if configure script exists: $(ls -la configure 2>/dev/null || echo 'NOT FOUND')" &&
        echo "==============================" &&
        
        # Only configure if config.status doesn't exist or is older than configure
        if [ ! -f config.status ] || [ configure -nt config.status ]; then
            echo "Running configure..."
            ./configure --prefix="''' + install_dir + '''" --disable-shared --enable-static --disable-gcs --disable-s3
        else
            echo "Configure is up to date, skipping..."
        fi
        
        # Create our stamp file with a checksum of key files
        echo "$(stat -c %Y configure 2>/dev/null || stat -f %m configure)" > "''' + build_dir + '''/configure.done"
        echo "Created stamp file: ''' + build_dir + '''/configure.done"
        '''
    ],
    output: 'configure.done',
    console: true,
)

# Step 2: Build htslib  
htslib_make = custom_target('htslib-make',
    input: [htslib_configure],  # Depends on configure step
    command: [
        'sh', '-c', '''
        cd ''' + source_dir + '''/../htslib-1.22 && 
        make -j8 && 
        touch "''' + build_dir + '''/make.done"
        '''
    ],
    output: 'make.done',
    console: true,
)

# Step 3: Install htslib
htslib_install = custom_target('htslib-install',
    input: [htslib_make],  # Depends on make step
    command: [
        'sh', '-c', '''
        cd ''' + source_dir + '''/../htslib-1.22 && 
        make install &&
        touch "''' + build_dir + '''/install.done"
        '''
    ],
    output: 'install.done',
    console: true,
)

# Step 4: Create the library artifact
htslib_lib = custom_target('htslib-lib',
    input: [htslib_install],  # Depends on install step
    command: [
        'cp', join_paths(install_dir, 'lib', 'libhts.a'), '@OUTPUT@'
    ],
    output: 'libhts.a',
    console: true,
)

# Create the dependency
htslib_inc = include_directories(join_paths('install', 'include'))
htslib_dep = declare_dependency(
    link_with: htslib_lib,
    include_directories: include_directories(join_paths('install', 'include')),
)

# This is what makes it available to `get_variable()`
meson.override_dependency('htslib', htslib_dep)
